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

  // ---- Calc Constants ----
  calc_physical_constant(par);

  // ---- components ----
  // output_.SetUp(par_, fld_);
  // bound_.SetUp(grd_, fld_, topo_, par_);

  // ---- Initialization ----
  initial_.Initialization(fld_, fid_);

  // ---- components ----
  // control_.SetUp(par_, 8);
}