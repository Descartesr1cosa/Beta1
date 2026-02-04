
#include "MercurySolver.h"

void MercurySolver::Build_E_explicit_edge_()
{
    AddIdealEdgeEMF_();

    AddHallEdgeEMF_();

    // AddAmbipolarEdgeEMF_();

    // 后续CT只会用到inner的电场，不会使用nghost区域，因此只需对Pole处理即可
    // bound_.add_Edge_pole_boundary("E_xi"); // pole边界处理
    // bound_.add_Edge_pole_boundary("E_eta");
}
