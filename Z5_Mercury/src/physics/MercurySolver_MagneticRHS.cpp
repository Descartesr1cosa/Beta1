#include "MercurySolver.h"

void MercurySolver::AssembleRHS_Induction_CT_()
{
    const int nb = fld_->num_blocks();

    // 0) E_edge 清零
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        auto zero_electric = [&](FieldBlock &E)
        {
            Int3 lo = E.get_lo(), hi = E.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };
        zero_electric(Exi);
        zero_electric(Eeta);
        zero_electric(Eze);
    }

    // Explicit Advance
    {
        // 1) Calculate the E fields
        Build_E_explicit_edge_();

        // 3) E_edge 也要进边界/halo（否则 curl 更新 B_face 时边界附近会乱）
        // mercury_bound_.Sync("Eedge"); // 你需要加一个 group：fields={E_xi,E_eta,E_zeta}

        // 4) curl(E_edge) -> RHS_Bface，然后 Bface += dt*RHS
        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Exi = fld_->field(fid_.fid_E.xi, ib);
            auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
            auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

            auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
            auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
            auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
            if (!Exi.is_allocated())
                continue;

            // multiper 的符号你要和你的取向一致：
            // 一般 CT 是 B^{n+1} = B^n - dt * curl(E)
            CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
        }
    }
}