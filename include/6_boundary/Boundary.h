#pragma once

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdio>

#include "6_boundary/Boundary_Type.h"

#include "2_topology/2_MPCNS_Topology.h"
#include "1_grid/1_MPCNS_Grid.h"

// 维护两个 registry：
// PhysicalRegistry phy_reg_：(location, field_name, bc_name) → PhysicalHandler
// CouplingRegistry cpl_reg_：(src, dst, location, channel_tag, dst_field) → CouplingHandler
// 提供 Apply：
// ApplyPhysical(field_name)：对这个 field 的所有物理 patch 施加 BC
// ApplyCouplingPair(src,dst)：把 coupling buffer 写入 dst ghost（并允许处理）
class BoundaryCore
{
public:
    BoundaryCore() = default;
    ~BoundaryCore() = default;

    // ------------------------------------------------------------
    // Setup
    // ------------------------------------------------------------
    void SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par);

    // ------------------------------------------------------------
    // Registry: Physical BC
    // ------------------------------------------------------------
    // 注册：对某个 location + 某个 field + 某个 bc_name 的处理器
    void RegisterPhysical(StaggerLocation loc,
                          const std::string &field_name,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);

    // 注册：对某个 location + 所有 field + 某个 bc_name 的通用处理器
    void RegisterPhysical(StaggerLocation loc,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);

    // 设置：对某个 location 的默认处理器（fallback）
    void SetDefaultPhysical(StaggerLocation loc, BOUND::PhysicalHandler h);

    // ------------------------------------------------------------
    // Registry: Coupling BC
    // ------------------------------------------------------------
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          const std::string &channel_tag,
                          const std::string &dst_field_name,
                          BOUND::CouplingHandler h);

    // src/dst/location 下所有 channel 的默认 handler
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          BOUND::CouplingHandler h);

    // src/dst 的全默认 handler（可选）
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          BOUND::CouplingHandler h);

    // ------------------------------------------------------------
    // Apply
    // ------------------------------------------------------------
    void ApplyPhysical(const std::string &field_name);
    void ApplyPhysical(const std::vector<std::string> &field_names);

    // 对一个 coupling pair 执行：把 buffer 写入 dst ghost
    void ApplyCouplingPair(const std::string &src, const std::string &dst);

protected:
    // ------------------------------------------------------------
    // You may override these defaults if needed
    // ------------------------------------------------------------
    static void DefaultPhysicalCopy(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int nghost);
    static void DefaultCouplingCopy(FieldBlock &Udst, Field *fld,
                                    CouplingBufferBlock &buf,
                                    const std::string &src,
                                    const std::string &dst,
                                    const std::string &channel_tag);

private:
    // ------------------------------------------------------------
    // Pointers
    // ------------------------------------------------------------
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;

    // ------------------------------------------------------------
    // Cached physical patch list (from topo_->physical_patches)
    //
    // 这里的 PhysicalRegion.base_box 我们存“node patch box”（p.this_box_node）
    // Apply 时按 field.location/nghost 动态算 ghost slab box，再临时覆盖 base_box 传给 handler。
    // ------------------------------------------------------------
    std::vector<BOUND::PhysicalRegion> phy_patches_;

    // ------------------------------------------------------------
    // Registries
    // ------------------------------------------------------------
    BOUND::PhysicalRegistry phy_reg_;
    BOUND::CouplingRegistry cpl_reg_;

    // ------------------------------------------------------------
    // Internal helpers: resolve handlers
    // ------------------------------------------------------------
    BOUND::PhysicalHandler ResolvePhysical(StaggerLocation loc,
                                           const std::string &field_name,
                                           const std::string &bc_name) const;

    BOUND::CouplingHandler ResolveCoupling(const std::string &src,
                                           const std::string &dst,
                                           StaggerLocation loc,
                                           const std::string &channel_tag,
                                           const std::string &dst_field_name) const;

    // ------------------------------------------------------------
    // Internal helpers: geometry for ghost slab on a boundary face
    // This copies the logic style used in your coupling buffer builder.
    // ------------------------------------------------------------
    static Int3 LocDelta(StaggerLocation loc);
    static Int3 LocInnerHi(const Block &blk, StaggerLocation loc);

    static void ConvertTangent(int lo_n, int hi_n, int delta, int &lo, int &hi);

    // 给定：block、field location、边界 face 的 node-box、dir_code(±1/±2/±3)、nghost
    // 返回：该 location 下的 ghost slab box（half-open [lo,hi)）
    static Box3 MakeFaceGhostSlabBox(const Block &blk,
                                     StaggerLocation loc,
                                     const Box3 &face_node_box,
                                     int dir_code,
                                     int nghost);

    // Build cached physical patch list
    void BuildPhysicalPatchCache();
};