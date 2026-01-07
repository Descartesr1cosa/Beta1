#include "6_boundary/Boundary.h"
#include <iostream>
#include <cstdlib>

// ------------------------------------------------------------
// Setup
// ------------------------------------------------------------
void BoundaryCore::SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par)
{
    grd_ = grd;
    fld_ = fld;
    topo_ = topo;
    par_ = par;

    if (!grd_ || !fld_ || !topo_ || !par_)
        throw std::runtime_error("[BoundaryCore] SetUp got null pointer");

    BuildPhysicalPatchCache();

    // 如果你不注册任何 handler，也要能跑：
    // 这里不强制插 registry，而是 Resolve 失败时直接回退到 Default*Copy。
}

// ------------------------------------------------------------
// Registry: Physical
// ------------------------------------------------------------
void BoundaryCore::RegisterPhysical(StaggerLocation loc,
                                    const std::string &field_name,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    BOUND::PhysicalKey k{loc, field_name, bc_name};
    phy_reg_[k] = std::move(h);
}

void BoundaryCore::RegisterPhysical(StaggerLocation loc,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    RegisterPhysical(loc, "", bc_name, std::move(h));
}

void BoundaryCore::SetDefaultPhysical(StaggerLocation loc, BOUND::PhysicalHandler h)
{
    // default 的约定：bc_name="" 且 field_name=""
    RegisterPhysical(loc, "", "", std::move(h));
}

// ------------------------------------------------------------
// Apply Physical
// ------------------------------------------------------------
void BoundaryCore::ApplyPhysical(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    const StaggerLocation loc = desc.location;
    const int nghost = desc.nghost;

    for (const auto &patch : phy_patches_)
    {
        const int ib = patch.this_block;
        FieldBlock &U = fld_->field(fid, ib);

        // 1) 把 topo 的 node patch box 转成该 loc 下的 ghost slab box
        const Block &blk = grd_->grids(ib);
        const Box3 ghost_box = MakeFaceGhostSlabBox(blk, loc, patch.base_box, patch.direction, nghost);

        // 2) 构造一个临时 region：把 base_box 覆盖成 ghost_box，传给 handler
        BOUND::PhysicalRegion r = patch;
        r.base_box = ghost_box;

        // 3) Resolve handler（不绑定，每次查 registry）
        auto h = ResolvePhysical(loc, field_name, patch.bc_name);
        if (h)
            h(U, fld_, r, nghost);
        else
            DefaultPhysicalCopy(U, fld_, r, nghost);
    }
}

void BoundaryCore::ApplyPhysical(const std::vector<std::string> &field_names)
{
    for (const auto &fn : field_names)
        ApplyPhysical(fn);
}

// ------------------------------------------------------------
// Resolve handlers (priority search)
// ------------------------------------------------------------
BOUND::PhysicalHandler BoundaryCore::ResolvePhysical(StaggerLocation loc,
                                                     const std::string &field_name,
                                                     const std::string &bc_name) const
{
    // 优先级：
    // 1) (loc, field, bc)
    // 2) (loc, "",    bc)
    // 3) (loc, field, "")
    // 4) (loc, "",    "")
    auto find_one = [&](const std::string &f, const std::string &b) -> BOUND::PhysicalHandler
    {
        BOUND::PhysicalKey k{loc, f, b};
        auto it = phy_reg_.find(k);
        if (it != phy_reg_.end())
            return it->second;
        return nullptr;
    };

    if (auto h = find_one(field_name, bc_name))
        return h;
    if (auto h = find_one("", bc_name))
        return h;
    if (auto h = find_one(field_name, ""))
        return h;
    if (auto h = find_one("", ""))
        return h;
    return nullptr;
}

// ------------------------------------------------------------
// Physical patch cache
// ------------------------------------------------------------
void BoundaryCore::BuildPhysicalPatchCache()
{
    phy_patches_.clear();
    phy_patches_.reserve(topo_->physical_patches.size());

    for (const auto &p : topo_->physical_patches)
    {
        BOUND::PhysicalRegion r;
        r.this_block = p.this_block;
        r.this_block_name = p.this_block_name;
        r.bc_id = p.bc_id;
        r.bc_name = p.bc_name;
        r.direction = p.direction;
        r.raw = p.raw;

        // cycle：如果 raw 提供就用 raw；否则也可从 direction 推导
        if (p.raw)
        {
            r.cycle.i = p.raw->cycle[0];
            r.cycle.j = p.raw->cycle[1];
            r.cycle.k = p.raw->cycle[2];
        }
        else
        {
            r.cycle = {0, 0, 0};
            const int ax = std::abs(p.direction);
            const int sgn = (p.direction > 0) ? +1 : -1;
            if (ax == 1)
                r.cycle.i = sgn;
            if (ax == 2)
                r.cycle.j = sgn;
            if (ax == 3)
                r.cycle.k = sgn;
        }

        // base_box：这里暂存 topo 的 node patch box
        r.base_box = p.this_box_node;

        phy_patches_.push_back(r);
    }
}

// ------------------------------------------------------------
// Default handlers
// ------------------------------------------------------------
void BoundaryCore::DefaultPhysicalCopy(FieldBlock &U, Field * /*fld*/,
                                       const BOUND::PhysicalRegion &r, int /*nghost*/)
{
    // 约定：r.base_box 在 ApplyPhysical 里已被覆盖成“ghost slab box”
    const Box3 &g = r.base_box;

    // 计算 interior 参考面索引：靠近边界的第一个 interior layer
    // 使用 U.block + U.desc.location 推出 interior hi
    const Block &blk = U.get_block();
    const Int3 hi_in = LocInnerHi(blk, U.descriptor().location);

    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;

    int i_ref = 0, j_ref = 0, k_ref = 0;
    // 默认用“相同 (j,k)”复制，normal 方向取 interior 的贴边层
    // 例如 X-: i_ref=0; X+: i_ref=hi_in.i-1
    auto ref_index = [&](int i, int j, int k, int &ir, int &jr, int &kr)
    {
        ir = i;
        jr = j;
        kr = k;
        if (ax == 1)
            ir = (sgn < 0) ? 0 : (hi_in.i - 1);
        if (ax == 2)
            jr = (sgn < 0) ? 0 : (hi_in.j - 1);
        if (ax == 3)
            kr = (sgn < 0) ? 0 : (hi_in.k - 1);
    };

    const int ncomp = U.descriptor().ncomp;

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                ref_index(i, j, k, i_ref, j_ref, k_ref);
                for (int m = 0; m < ncomp; ++m)
                    U(i, j, k, m) = U(i_ref, j_ref, k_ref, m);
            }
}
