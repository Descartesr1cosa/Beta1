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

    BuildPhysicalPatternsCachedInnerSlabs();

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

    auto pit = phy_patterns_.find(loc);
    if (pit == phy_patterns_.end())
        return;

    for (const auto &cached : pit->second.regions)
    {
        FieldBlock &U = fld_->field(fid, cached.this_block);

        // 运行时仅 O(1)：inner_slab -> ghost_slab
        BOUND::PhysicalRegion work = cached;
        work.box = MakeGhostSlabFromInner(cached.inner_slab, cached.direction, nghost);

        auto h = ResolvePhysical(loc, field_name, cached.bc_name);
        if (h)
            h(U, fld_, work, nghost);
        else
            DefaultPhysicalCopy(U, fld_, work, nghost);
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
void BoundaryCore::BuildPhysicalPatternsCachedInnerSlabs()
{
    phy_patterns_.clear();

    // 只构建你当前会用到的 location；这里示例构建常用 7 类
    const std::vector<StaggerLocation> locs = {
        StaggerLocation::Cell,
        StaggerLocation::FaceXi, StaggerLocation::FaceEt, StaggerLocation::FaceZe,
        StaggerLocation::EdgeXi, StaggerLocation::EdgeEt, StaggerLocation::EdgeZe};

    for (auto loc : locs)
    {
        BOUND::PhysicalPattern pat;
        pat.location = loc;
        phy_patterns_[loc] = std::move(pat);
    }

    for (const auto &p : topo_->physical_patches)
    {
        const int ib = p.this_block;
        const Block &blk = grd_->grids(ib); // 你工程里取 block 的接口按实际改

        const Box3 face_node_box = p.this_box_node; // topo 给的 node patch box（半开区间）
        const int dir = p.direction;

        for (auto loc : locs)
        {
            BOUND::PhysicalRegion r;
            r.this_block = p.this_block;
            r.this_block_name = p.this_block_name;
            r.bc_id = p.bc_id;
            r.bc_name = p.bc_name;
            r.direction = dir;
            r.raw = p.raw;

            // cycle 可选：也可完全由 direction 推导
            if (p.raw)
            {
                r.cycle.i = p.raw->cycle[0];
                r.cycle.j = p.raw->cycle[1];
                r.cycle.k = p.raw->cycle[2];
            }

            // 关键：Build 阶段推导“域内贴边1层”的 inner_slab（loc 坐标）
            r.inner_slab = MakeInnerSlabBox_OneLayer(blk, loc, face_node_box, dir);

            phy_patterns_[loc].regions.push_back(std::move(r));
        }
    }
}

// ------------------------------------------------------------
// Default handlers
// ------------------------------------------------------------
void BoundaryCore::DefaultPhysicalCopy(FieldBlock &U, Field * /*fld*/,
                                       const BOUND::PhysicalRegion &r, int /*nghost*/)
{
    const Box3 &g = r.box;            // ghost slab：需要写入
    const Box3 &inner = r.inner_slab; // 域内贴边一层：参考

    const int ax = std::abs(r.direction); // 1/2/3
    const int sgn = (r.direction > 0) ? +1 : -1;

    // 法向参考索引：inner slab 的那一层（厚度=1）
    const int i_ref = (ax == 1) ? ((sgn < 0) ? inner.lo.i : (inner.hi.i - 1)) : 0;
    const int j_ref = (ax == 2) ? ((sgn < 0) ? inner.lo.j : (inner.hi.j - 1)) : 0;
    const int k_ref = (ax == 3) ? ((sgn < 0) ? inner.lo.k : (inner.hi.k - 1)) : 0;

    const int ncomp = U.descriptor().ncomp;

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                const int ii = (ax == 1) ? i_ref : i;
                const int jj = (ax == 2) ? j_ref : j;
                const int kk = (ax == 3) ? k_ref : k;

                for (int m = 0; m < ncomp; ++m)
                    U(i, j, k, m) = U(ii, jj, kk, m);
            }
}