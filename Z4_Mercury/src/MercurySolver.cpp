// Core
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

// Z4_Mercury
#include "MercurySolver.h"

MercurySolver::MercurySolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo,
                             Param *par)
    : grd_(grd),
      topo_(topo),
      fld_(fld),
      halo_(halo),
      par_(par)
{
  // ---- Cache field ids ----
  fid_.Init(fld_);

  // ---- Build IO Module ----
  constexpr int NRES = 13; // 只统计H Na 守恒变量和感应磁场
  io_.Setup(par_, grd_, fld_, NRES);

  {
    std::vector<std::string> bin_name = {"U_H", "U_Na", "U_b"};
    io_.ClearRestartFields();
    io_.SetRestartFields(bin_name);

    io_.SetTecplotMode(IOModule::TecplotMode::CellAsNode);
    std::vector<std::string> tec_block_name = {}; // 全部物理块输出
    io_.SetTecplotBlock(tec_block_name);

    std::vector<std::string> plt_name = {"PV_H", "PV_Na", "B_cell"};
    io_.SetTecplotFields(plt_name);

    std::string fld_name = "PV_H";
    std::vector<std::string> var_name = {"u_H", "v_H", "w_H", "p_H", "T_H"};
    io_.SetTecplotFieldComponentNames(fld_name, var_name);

    fld_name = "PV_Na";
    var_name = {"u_Na", "v_Na", "w_Na", "p_Na", "T_Na"};
    io_.SetTecplotFieldComponentNames(fld_name, var_name);

    fld_name = "B_cell";
    var_name = {"Bx", "By", "Bz"};
    io_.SetTecplotFieldComponentNames(fld_name, var_name);
  }

  // ---- Calc Constants ----
  calc_physical_constant(par_);

  // ---- components ----
  {
    // 0) 需要添加边界的物理场
    std::vector<std::string> bnd_fields = {"U_H", "U_Na", "U_b"};

    // 1) 初始化 BoundaryCore：Build 阶段会按 location 缓存每个 patch 的 inner_slab（法向1层）
    bound_.SetUp(grd_, fld_, topo_, par_, bnd_fields);

    // 2) 注册“默认物理边界处理”：至少把你会用到的 location 都设 default
    //    做法：location + ("","") 作为最终 fallback
    // for (const auto &fn : bnd_fields)
    // {
    //   bound_.RegisterPhysical(fn, "",
    //                           [](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    //                           {
    //                             BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    //                           });
    // }

    // 3) 注册“特定边界类型”的处理（示例：Cell + U_ + Solid_Surface）
    bound_.RegisterPhysical("U_H", "Outflow",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_H", "Pole",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_H", "Farfield",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              this->BC_UH_Farfield_H(U, fld, r, ngh);
                            });
    bound_.RegisterPhysical("U_Na", "Outflow",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_Na", "Pole",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_Na", "Farfield",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              this->BC_UH_Farfield_Na(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_b", "Outflow",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_b", "Pole",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
                            });

    bound_.RegisterPhysical("U_b", "Farfield",
                            [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                            {
                              this->BC_UH_Farfield_b(U, fld, r, ngh);
                            });

    // 4) coupling handler 注册（可选：你也可以先不注册，让默认 copy 工作）
    // 注意：coupling buffers 的 build 发生在 Field 注册 coupling channels 并 build buffers 之后
    // 这里仅注册 handler，不需要绑定
    bound_.RegisterCoupling("Solid", "Fluid", StaggerLocation::Cell, "U_b", "U_b",
                            [](FieldBlock &Udst, Field *fld, CouplingBufferBlock &buf,
                               const std::string &src, const std::string &dst, const std::string &tag)
                            {
                              BoundaryCore::DefaultCouplingCopy(Udst, fld, buf, src, dst, tag);
                            });
    bound_.RegisterCoupling("Fluid", "Solid", StaggerLocation::Cell, "U_b", "U_b",
                            [](FieldBlock &Udst, Field *fld, CouplingBufferBlock &buf,
                               const std::string &src, const std::string &dst, const std::string &tag)
                            {
                              BoundaryCore::DefaultCouplingCopy(Udst, fld, buf, src, dst, tag);
                            });

    // 5) 严格检查：缺失就 throw（你想要的模式）
    bound_.CheckPhysicalHandlers(bnd_fields);
  }

  // ---- Initialization ----
  if (par_->GetBoo("continue_calc"))
  {
    io_.ReadRestartBinFile();
    io_.ReadRunDataFile();
  }
  initial_.Initialization(fld_, fid_);

  // ---- components ----
  control_.Setup(par_);
}