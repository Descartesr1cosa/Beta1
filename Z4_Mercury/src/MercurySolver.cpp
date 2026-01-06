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
  io_.Setup(par_, grd_, fld_, 13);

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
  calc_physical_constant(par);

  // ---- components ----
  // bound_.SetUp(grd_, fld_, topo_, par_);

  // ---- Initialization ----
  initial_.Initialization(fld_, fid_);

  // ---- components ----
  // control_.SetUp(par_, 8);
}