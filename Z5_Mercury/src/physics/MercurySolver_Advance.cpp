
#include "MercurySolver.h"

void MercurySolver::Time_Advance()
{

    // ---------- Hall 子循环数 ----------
    const double safety = 0.8;
    int nsub = 1;
    // if (dt_hall > 0.0)
    //     nsub = std::max(1, (int)std::ceil(dt / (safety * dt_hall)));
    // nsub = std::min(nsub, 200);
    dt_sub = dt / nsub;
    if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    {
        std::printf("[HallSub] step=%d dt=%.3e dt_hall=%.3e nsub=%d dt_sub=%.3e\n\n",
                    run_data_->step, dt, dt_hall, nsub, dt / double(nsub));
    }

    // 1) dq/db set to ZERO
    ZeroRHS_();

    // 2) RHS for Flow
    AssembleRHS_Fluid_(); // = Scheme_U_ + AddSourceToRHS_Fluid()

    // 3) RHS for Magnetic fields：CT (Edge EMF -> Curl -> Face RHS)
    AssembleRHS_Induction_CT_();

    // 4) Euler 1st order Advance
    ApplyUpdate_Euler_(); // U += dt*RHS, B_face += dt*RHS_Bface

    // 5) 低密度/负压修复（尽量按 Fortran：邻域平均 + 重建 E）
    // RepairNonPhysical_();

    // 2) Hall 子步：只更新 Bface
    mercury_bound_.Sync("Ucell");
    mercury_bound_.Sync("Bface");

    calc_PV();
    calc_Uplus();

    {
        const int nb2 = fld_->num_blocks();
        // --------------------- persistent buffers (allocate once) ---------------------
        struct Buf
        {
            bool ok{false};
            Int3 lo, hi;
            std::vector<double> v; // inner only
        };

        auto make_buf_like = [&](FieldBlock &F) -> Buf
        {
            Buf b;
            if (!F.is_allocated())
                return b;
            b.ok = true;
            b.lo = F.inner_lo();
            b.hi = F.inner_hi();
            const int ni = b.hi.i - b.lo.i;
            const int nj = b.hi.j - b.lo.j;
            const int nk = b.hi.k - b.lo.k;
            b.v.assign((size_t)ni * (size_t)nj * (size_t)nk, 0.0);
            return b;
        };

        auto pack0 = [&](FieldBlock &F, Buf &b)
        {
            if (!b.ok)
                return;
            size_t idx = 0;
            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k, ++idx)
                        b.v[idx] = F(i, j, k, 0);
        };

        auto pack_rhs = [&](FieldBlock &RHS, Buf &b)
        {
            if (!b.ok)
                return;
            size_t idx = 0;
            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k, ++idx)
                        b.v[idx] = RHS(i, j, k, 0);
        };

        // write: B = B0 + (fac * dt_step) * K
        auto write_B = [&](FieldBlock &B, const Buf &B0, const Buf &K, double dt_step, double fac)
        {
            if (!B0.ok)
                return;
            const double a = fac * dt_step;
            size_t idx = 0;
            for (int i = B0.lo.i; i < B0.hi.i; ++i)
                for (int j = B0.lo.j; j < B0.hi.j; ++j)
                    for (int k = B0.lo.k; k < B0.hi.k; ++k, ++idx)
                        B(i, j, k, 0) = B0.v[idx] + a * K.v[idx];
        };

        // S += w * K
        auto axpy = [&](Buf &S, const Buf &K, double w)
        {
            if (!S.ok)
                return;
            for (size_t t = 0; t < S.v.size(); ++t)
                S.v[t] += w * K.v[t];
        };

        auto zero_buf = [&](Buf &b)
        {
            if (!b.ok)
                return;
            std::fill(b.v.begin(), b.v.end(), 0.0);
        };

        std::vector<Buf> B0_xi(nb2), B0_eta(nb2), B0_ze(nb2);
        std::vector<Buf> K_xi(nb2), K_eta(nb2), K_ze(nb2);
        std::vector<Buf> S_xi(nb2), S_eta(nb2), S_ze(nb2);

        // allocate once based on Bface inner extents
        for (int ib = 0; ib < nb2; ++ib)
        {
            auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
            auto &Beta = fld_->field(fid_.fid_B.eta, ib);
            auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

            B0_xi[ib] = make_buf_like(Bxi);
            B0_eta[ib] = make_buf_like(Beta);
            B0_ze[ib] = make_buf_like(Bze);

            K_xi[ib] = make_buf_like(Bxi);
            K_eta[ib] = make_buf_like(Beta);
            K_ze[ib] = make_buf_like(Bze);

            S_xi[ib] = make_buf_like(Bxi);
            S_eta[ib] = make_buf_like(Beta);
            S_ze[ib] = make_buf_like(Bze);
        }

        auto rk4_hall_bface_only = [&](double dt_step)
        {
            // snapshot B0, reset S
            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
                auto &Beta = fld_->field(fid_.fid_B.eta, ib);
                auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

                if (B0_xi[ib].ok)
                    pack0(Bxi, B0_xi[ib]);
                if (B0_eta[ib].ok)
                    pack0(Beta, B0_eta[ib]);
                if (B0_ze[ib].ok)
                    pack0(Bze, B0_ze[ib]);

                zero_buf(S_xi[ib]);
                zero_buf(S_eta[ib]);
                zero_buf(S_ze[ib]);
            }

            // ---------------- Stage 1: K1 = f(B0); S = K1; B = B0 + 0.5 dt K1
            mercury_bound_.Sync("Bface");
            AssembleRHS_Induction_CT_HallOnly_();

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
                auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
                auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

                if (K_xi[ib].ok)
                {
                    pack_rhs(RHSBxi, K_xi[ib]);
                    S_xi[ib].v = K_xi[ib].v;
                }
                if (K_eta[ib].ok)
                {
                    pack_rhs(RHSBeta, K_eta[ib]);
                    S_eta[ib].v = K_eta[ib].v;
                }
                if (K_ze[ib].ok)
                {
                    pack_rhs(RHSBze, K_ze[ib]);
                    S_ze[ib].v = K_ze[ib].v;
                }
            }

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
                auto &Beta = fld_->field(fid_.fid_B.eta, ib);
                auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

                if (B0_xi[ib].ok)
                    write_B(Bxi, B0_xi[ib], K_xi[ib], dt_step, 0.5);
                if (B0_eta[ib].ok)
                    write_B(Beta, B0_eta[ib], K_eta[ib], dt_step, 0.5);
                if (B0_ze[ib].ok)
                    write_B(Bze, B0_ze[ib], K_ze[ib], dt_step, 0.5);
            }

            // ---------------- Stage 2: K2 = f(B0+0.5dtK1); S += 2K2; B = B0 + 0.5 dt K2
            mercury_bound_.Sync("Bface");
            AssembleRHS_Induction_CT_HallOnly_();

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
                auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
                auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

                if (K_xi[ib].ok)
                {
                    pack_rhs(RHSBxi, K_xi[ib]);
                    axpy(S_xi[ib], K_xi[ib], 2.0);
                }
                if (K_eta[ib].ok)
                {
                    pack_rhs(RHSBeta, K_eta[ib]);
                    axpy(S_eta[ib], K_eta[ib], 2.0);
                }
                if (K_ze[ib].ok)
                {
                    pack_rhs(RHSBze, K_ze[ib]);
                    axpy(S_ze[ib], K_ze[ib], 2.0);
                }
            }

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
                auto &Beta = fld_->field(fid_.fid_B.eta, ib);
                auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

                if (B0_xi[ib].ok)
                    write_B(Bxi, B0_xi[ib], K_xi[ib], dt_step, 0.5);
                if (B0_eta[ib].ok)
                    write_B(Beta, B0_eta[ib], K_eta[ib], dt_step, 0.5);
                if (B0_ze[ib].ok)
                    write_B(Bze, B0_ze[ib], K_ze[ib], dt_step, 0.5);
            }

            // ---------------- Stage 3: K3 = f(B0+0.5dtK2); S += 2K3; B = B0 + 1.0 dt K3
            mercury_bound_.Sync("Bface");
            AssembleRHS_Induction_CT_HallOnly_();

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
                auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
                auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

                if (K_xi[ib].ok)
                {
                    pack_rhs(RHSBxi, K_xi[ib]);
                    axpy(S_xi[ib], K_xi[ib], 2.0);
                }
                if (K_eta[ib].ok)
                {
                    pack_rhs(RHSBeta, K_eta[ib]);
                    axpy(S_eta[ib], K_eta[ib], 2.0);
                }
                if (K_ze[ib].ok)
                {
                    pack_rhs(RHSBze, K_ze[ib]);
                    axpy(S_ze[ib], K_ze[ib], 2.0);
                }
            }

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
                auto &Beta = fld_->field(fid_.fid_B.eta, ib);
                auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

                if (B0_xi[ib].ok)
                    write_B(Bxi, B0_xi[ib], K_xi[ib], dt_step, 1.0);
                if (B0_eta[ib].ok)
                    write_B(Beta, B0_eta[ib], K_eta[ib], dt_step, 1.0);
                if (B0_ze[ib].ok)
                    write_B(Bze, B0_ze[ib], K_ze[ib], dt_step, 1.0);
            }

            // ---------------- Stage 4: K4 = f(B0+dtK3); S += K4; final B = B0 + (dt/6) S
            mercury_bound_.Sync("Bface");
            AssembleRHS_Induction_CT_HallOnly_();

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
                auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
                auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

                if (K_xi[ib].ok)
                {
                    pack_rhs(RHSBxi, K_xi[ib]);
                    axpy(S_xi[ib], K_xi[ib], 1.0);
                }
                if (K_eta[ib].ok)
                {
                    pack_rhs(RHSBeta, K_eta[ib]);
                    axpy(S_eta[ib], K_eta[ib], 1.0);
                }
                if (K_ze[ib].ok)
                {
                    pack_rhs(RHSBze, K_ze[ib]);
                    axpy(S_ze[ib], K_ze[ib], 1.0);
                }
            }

            for (int ib = 0; ib < nb2; ++ib)
            {
                auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
                auto &Beta = fld_->field(fid_.fid_B.eta, ib);
                auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

                if (B0_xi[ib].ok)
                    write_B(Bxi, B0_xi[ib], S_xi[ib], dt_step, 1.0 / 6.0);
                if (B0_eta[ib].ok)
                    write_B(Beta, B0_eta[ib], S_eta[ib], dt_step, 1.0 / 6.0);
                if (B0_ze[ib].ok)
                    write_B(Bze, B0_ze[ib], S_ze[ib], dt_step, 1.0 / 6.0);
            }

            mercury_bound_.Sync("Bface"); // for next substep J=curl(B)
        };

        for (int s = 0; s < nsub; ++s)
        {
            rk4_hall_bface_only(dt_sub);
        }
    }

    // for (int s = 0; s < nsub; ++s)
    // {
    //     // 只组装 Hall 的 RHS_b（不动 U 的 RHS）
    //     AssembleRHS_Induction_CT_HallOnly_();

    //     // 只更新 Bface: Bface += dt_sub * RHS_b
    //     ApplyUpdate_Euler_BfaceOnly_(dt_sub);

    //     // 更新后做一次 Bface 同步，供下一个子步算 J=curl(B)
    //     mercury_bound_.Sync("Bface");
    // }
}

void MercurySolver::ZeroRHS_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_et = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_ze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; // = Jac.inner_lo();
        Int3 hi; // = Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = RHSH.inner_lo();
            hi = RHSH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSH(i, j, k, m) = 0.0;
        }

        if (RHSN.is_allocated())
        {
            lo = RHSN.inner_lo();
            hi = RHSN.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSN(i, j, k, m) = 0.0;
        }

        if (RHSB_xi.is_allocated())
        {
            lo = RHSB_xi.inner_lo();
            hi = RHSB_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_xi(i, j, k, 0) = 0.0;
        }

        if (RHSB_et.is_allocated())
        {
            lo = RHSB_et.inner_lo();
            hi = RHSB_et.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_et(i, j, k, 0) = 0.0;
        }

        if (RHSB_ze.is_allocated())
        {
            lo = RHSB_ze.inner_lo();
            hi = RHSB_ze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_ze(i, j, k, 0) = 0.0;
        }
    }
}

void MercurySolver::ApplyUpdate_Euler_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        FieldBlock &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        FieldBlock &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; //= Jac.inner_lo();
        Int3 hi; //= Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = UH.inner_lo();
            hi = UH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UH(i, j, k, m) += dt * RHSH(i, j, k, m);
        }

        if (RHSN.is_allocated())
        {
            lo = UN.inner_lo();
            hi = UN.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UN(i, j, k, m) += dt * RHSN(i, j, k, m);
        }

        if (RHSB_xi.is_allocated())
        {
            lo = Ub_xi.inner_lo();
            hi = Ub_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt * RHSB_xi(i, j, k, 0);
        }

        if (RHSB_eta.is_allocated())
        {
            lo = Ub_eta.inner_lo();
            hi = Ub_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt * RHSB_eta(i, j, k, 0);
        }

        if (RHSB_zeta.is_allocated())
        {
            lo = Ub_zeta.inner_lo();
            hi = Ub_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt * RHSB_zeta(i, j, k, 0);
        }
    }
}

// Hall-only CT RHS assembly: builds RHSB_* from Ehall only (E is cleared first).
// Assumptions:
//   - Calc_J_Edge() computes Jxi/Jeta/Jzeta from current Bface and does needed Sync("Jedge") internally.
//   - AddHallEdgeEMF_() adds Ehall to E (so E must be zeroed first).
//   - mercury_bound_.Sync("Eedge") exists and syncs E_xi/E_eta/E_zeta (halo + physical BC/pole if you configured it).
//   - CTOperators::CurlEdgeToFace(..., multiper=-1.0) produces RHS consistent with B^{n+1}=B^n+dt*RHS.
// void MercurySolver::AssembleRHS_Induction_CT_HallOnly_()
// {
//     const int nb = fld_->num_blocks();

//     // 0) Clear E_edge
//     for (int ib = 0; ib < nb; ++ib)
//     {
//         auto &Exi = fld_->field(fid_.fid_E.xi, ib);
//         auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
//         auto &Eze = fld_->field(fid_.fid_E.zeta, ib);
//         if (!Exi.is_allocated())
//             continue;

//         auto zero_electric = [&](FieldBlock &E)
//         {
//             Int3 lo = E.get_lo(), hi = E.get_hi();
//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         E(i, j, k, 0) = 0.0;
//         };

//         zero_electric(Exi);
//         zero_electric(Eeta);
//         zero_electric(Eze);
//     }

//     // 1) Ensure Bface is ready for J = curl(B)
//     // mercury_bound_.Sync("Bface");

//     // 2) J_edge from current Bface
//     Calc_J_Edge();
//     // mercury_bound_.DebugCheckSurfaceTangentialEdgeField(
//     //     "J-inner-solid", "J_xi", "J_eta", "J_zeta", "Coupled-Solid");

//     // mercury_bound_.DebugCheckSurfaceTangentialEdgeFieldGhost(
//     //     "J-ghost-solid", "J_xi", "J_eta", "J_zeta", "Coupled-Solid", 1);

//     // 3) Add Hall EMF into E (E is currently zero => E = Ehall)
//     AddHallEdgeEMF_();

//     // 4) Sync Eedge before taking curl(E)
//     mercury_bound_.Sync("Eedge");

//     // mercury_bound_.DebugCheckSurfaceTangentialEdgeField(
//     //     "E-inner-solid", "E_xi", "E_eta", "E_zeta", "Coupled-Solid");

//     // mercury_bound_.DebugCheckSurfaceTangentialEdgeFieldGhost(
//     //     "E-ghost-solid", "E_xi", "E_eta", "E_zeta", "Coupled-Solid", 1);

//     // 5) Clear RHSB_* (important if CurlEdgeToFace writes += instead of =)
//     for (int ib = 0; ib < nb; ++ib)
//     {
//         auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
//         auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
//         auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
//         if (!RHSBxi.is_allocated())
//             continue;

//         auto zero_rhs = [&](FieldBlock &F)
//         {
//             Int3 lo = F.inner_lo(), hi = F.inner_hi();
//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         F(i, j, k, 0) = 0.0;
//         };

//         zero_rhs(RHSBxi);
//         zero_rhs(RHSBeta);
//         zero_rhs(RHSBze);
//     }

//     // 6) curl(Ehall) -> RHS_Bface
//     for (int ib = 0; ib < nb; ++ib)
//     {
//         auto &Exi = fld_->field(fid_.fid_E.xi, ib);
//         auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
//         auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

//         auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
//         auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
//         auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

//         if (!Exi.is_allocated())
//             continue;

//         // B^{n+1} = B^n - dt * curl(E)  ==> RHS = -curl(E)
//         CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
//     }
// }

void MercurySolver::AssembleRHS_Induction_CT_HallOnly_()
{
    const int nb = fld_->num_blocks();

    double dt_stage = dt_sub;

    // 0) Clear E_edge
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

    // 1) J_edge from current Bface
    Calc_J_Edge();

    // 2) E = Ehall
    AddHallEdgeEMF_();

    // 3) Sync Eedge
    mercury_bound_.Sync("Eedge");

    // 4) Clear RHSB
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
        if (!RHSBxi.is_allocated())
            continue;

        auto zero_rhs = [&](FieldBlock &F)
        {
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        F(i, j, k, 0) = 0.0;
        };

        zero_rhs(RHSBxi);
        zero_rhs(RHSBeta);
        zero_rhs(RHSBze);
    }

    // ---- dB limiter parameters ----
    const double C_dB_far = 0.01;  // 先保守，建议先 0.005~0.01
    const double B_floor = 1.0e-2; // 按你的无量纲量级可再调
    const double tiny = 1.0e-300;

    double clip_ratio_max_l = 0.0;
    double clip_count_l = 0;

    // 5) curl(Ehall) -> RHS_Bface, then limit dB directly
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        // RHS = -curl(Ehall)
        CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, -1.0);

        auto limit_rhs_face = [&](FieldBlock &Bf, FieldBlock &RHSf)
        {
            Int3 lo = RHSf.inner_lo(), hi = RHSf.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // const double Bref = std::max(std::abs(Bf(i, j, k, 0)), B_floor);
                        const double Bref = std::max(std::abs(B0f_ref(i, j, k)), B_floor);
                        const double dB = std::abs(dt_stage * RHSf(i, j, k, 0));
                        const double dBmax = C_dB_far * Bref;

                        if (dB > dBmax)
                        {
                            const double fac = dBmax / (dB + tiny);
                            RHSf(i, j, k, 0) *= fac;

                            clip_count_l += 1.0;
                            clip_ratio_max_l = std::max(clip_ratio_max_l, dB / (dBmax + tiny));
                        }
                    }
        };

        limit_rhs_face(Bxi, RHSBxi);
        limit_rhs_face(Beta, RHSBeta);
        limit_rhs_face(Bze, RHSBze);
    }

    double clip_ratio_max_g = clip_ratio_max_l;
    PARALLEL::mpi_max(&clip_ratio_max_l, &clip_ratio_max_g, 1);

    double clip_count_g = clip_count_l;
    PARALLEL::mpi_sum(&clip_count_l, &clip_count_g, 1);

    if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    {
        std::printf("[Hall-dB-Limiter] step=%d  clips=%lf  maxRatio=%.3e  C_dB=%.3e  dt_stage=%.3e\n",
                    run_data_->step, clip_count_g, clip_ratio_max_g, C_dB_far, dt_stage);
        std::fflush(stdout);
    }
}
// Update only Bface with a supplied dt_step (used for Hall subcycling).
void MercurySolver::ApplyUpdate_Euler_BfaceOnly_(double dt_step)
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        FieldBlock &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        FieldBlock &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo, hi;

        if (RHSB_xi.is_allocated())
        {
            lo = Ub_xi.inner_lo();
            hi = Ub_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt_step * RHSB_xi(i, j, k, 0);
        }

        if (RHSB_eta.is_allocated())
        {
            lo = Ub_eta.inner_lo();
            hi = Ub_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt_step * RHSB_eta(i, j, k, 0);
        }

        if (RHSB_zeta.is_allocated())
        {
            lo = Ub_zeta.inner_lo();
            hi = Ub_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt_step * RHSB_zeta(i, j, k, 0);
        }
    }
}