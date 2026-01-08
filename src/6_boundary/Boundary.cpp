#include "6_boundary/Boundary.h"
#include <iostream>
#include <cstdlib>

// ------------------------------------------------------------
// Setup
// ------------------------------------------------------------
void BoundaryCore::SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par, const std::vector<std::string> &field_names)
{
    grd_ = grd;
    fld_ = fld;
    topo_ = topo;
    par_ = par;

    if (!grd_ || !fld_ || !topo_ || !par_)
        throw std::runtime_error("[BoundaryCore] SetUp got null pointer");

    // 1) 从 field_ids 推导需要的 locations
    enabled_locs_.clear();
    for (auto field : field_names)
    {
        int fid = fld_->field_id(field);
        const auto &desc = fld_->descriptor(fid);
        enabled_locs_.insert(desc.location);
    }

    // 2) 只为这些 locations 构建 Pattern（inner_slab，不扩 ghost）
    BuildPhysicalPatterns();
}

// ------------------------------------------------------------
// Registry: Physical
// ------------------------------------------------------------
void BoundaryCore::RegisterPhysical(StaggerLocation loc,
                                    const std::string &field_name,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    phy_reg_[BOUND::PhysicalKey{loc, field_name, bc_name}] = std::move(h);
}

void BoundaryCore::RegisterPhysical(const std::string &field_name,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    if (!fld_)
        throw std::runtime_error("[BoundaryCore] RegisterPhysical called before SetUp");

    const int fid = fld_->field_id(field_name);
    const auto &desc = fld_->descriptor(fid);
    const StaggerLocation loc = desc.location;

    // 强制要求：SetUp(field_ids) 已经包含该 field 的 location（更严格：包含该 fid）
    if (enabled_locs_.find(loc) == enabled_locs_.end())
        throw std::runtime_error("[BoundaryCore] RegisterPhysical: location not enabled by SetUp(field_ids).");

    RegisterPhysical(loc, field_name, bc_name, std::move(h));
}

void BoundaryCore::CheckPhysicalHandlers(const std::vector<std::string> &field_names) const
{
    // 0) 基本一致性检查
    if (!fld_ || !topo_ || !grd_)
        throw std::runtime_error("[BoundaryCore] AssertPhysicalHandlers called before SetUp");

    // 1) 收集：loc -> (bc_name -> sample_region*)
    std::map<StaggerLocation, std::map<std::string, const BOUND::PhysicalRegion *>> bcset;
    for (const auto &[loc, pat] : phy_patterns_)
    {
        for (const auto &r : pat.regions)
        {
            if (!bcset[loc].count(r.bc_name))
                bcset[loc][r.bc_name] = &r; // 存一个样例用于报错，如果没找到，可用该记录输出信息
        }
    }

    struct Missing
    {
        std::string field;
        StaggerLocation loc{};
        std::string bc;
        int example_block = -1;
        int bc_id = -1;
        int direction = 0;
    };

    std::vector<Missing> missing;

    // 2) 对每个 field 检查覆盖情况
    for (const auto &field_name : field_names)
    {
        const int fid = fld_->field_id(field_name);
        const auto &desc = fld_->descriptor(fid);
        const StaggerLocation loc = desc.location;

        // 2.1) pattern 必须存在，否则 SetUp(field_ids) 或 BuildPhysicalPatterns 有问题
        auto pit = phy_patterns_.find(loc);
        if (pit == phy_patterns_.end())
        {
            throw std::runtime_error("[BoundaryCore] No physical pattern built for location of field: " + field_name);
        }

        // 2.2) 该 loc 下出现过的每个 bc_name，都必须能找到 handler
        auto bit = bcset.find(loc);
        if (bit == bcset.end())
        {
            // 该 location 没有任何 patch（可能是全周期、或没有外边界），这不算错
            continue;
        }

        for (const auto &[bc_name, sample] : bit->second)
        {
            BOUND::PhysicalKey key{loc, field_name, bc_name};
            if (phy_reg_.find(key) == phy_reg_.end())
            {
                missing.push_back(Missing{
                    field_name, loc, bc_name,
                    sample ? sample->this_block : -1,
                    sample ? sample->bc_id : -1,
                    sample ? sample->direction : 0});
            }
        }
    }

    // 3) 缺失则严格停止
    if (!missing.empty())
    {
        std::string msg;
        msg += "[BoundaryCore] Missing Physical BC handlers:\n";
        for (const auto &m : missing)
        {
            msg += "  field=" + m.field + " loc=" + std::to_string((int)m.loc) + " bc_name=" + m.bc + " (example: block=" + std::to_string(m.example_block) + " bc_id=" + std::to_string(m.bc_id) + " dir=" + std::to_string(m.direction) + ")\n";
        }
        throw std::runtime_error(msg);
        std::abort();
        // 或者你想“直接停止”也可以：std::cerr<<msg; std::abort();
    }
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
        throw std::runtime_error("[BoundaryCore] ApplyPhysical: pattern not built for this field/location: " + field_name);

    for (const auto &cached : pit->second.regions)
    {
        FieldBlock &U = fld_->field(fid, cached.this_block);

        // 运行时仅 O(1)：inner_slab -> ghost_slab
        BOUND::PhysicalRegion work = cached;
        work.box = MakeGhostSlabFromInner(cached.inner_slab, cached.direction, nghost);

        auto h = ResolvePhysical(loc, field_name, cached.bc_name);
        if (!h)
            throw std::runtime_error("[BoundaryCore] ApplyPhysical: missing handler for field=" + field_name + " bc=" + cached.bc_name);
        h(U, fld_, work, nghost);
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
    // if (auto h = find_one("", bc_name))
    //     return h;
    // if (auto h = find_one(field_name, ""))
    //     return h;
    // if (auto h = find_one("", ""))
    //     return h;
    return nullptr;
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