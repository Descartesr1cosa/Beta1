
#include "MercurySolver.h"

void MercurySolver::Build_E_explicit_edge_()
{
    AddIdealEdgeEMF_();

    Calc_J_Edge();
    AddHallEdgeEMF_();
    AddResistiveEdgeEMF_(); // Add magnetic diffusion in solid (and optionally fluid)

    // AddAmbipolarEdgeEMF_();

    // 后续CT只会用到inner的电场，不会使用nghost区域，因此只需对Pole处理即可
    // bound_.add_Edge_pole_boundary("E_xi"); // pole边界处理
    // bound_.add_Edge_pole_boundary("E_eta");
}

void MercurySolver::Calc_J_Edge()
{
    //  ComputeJ_AtEdges_Inner_();
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);
    }

    mercury_bound_.Sync("Jedge"); // ApplyBC_EdgeJ_();
}
