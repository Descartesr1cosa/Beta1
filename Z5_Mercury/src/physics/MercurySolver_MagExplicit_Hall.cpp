#include "MercurySolver.h"

void MercurySolver::AddHallEdgeEMF_()
{
    ComputeJ_AtEdges_Inner_();
    // ApplyBC_EdgeJ_();
    AddHallE_AtEdges_EnergyPreserving_(); // 只填 Ehall_xi/eta/zeta（线积分量）加入E
}

void MercurySolver::ComputeJ_AtEdges_Inner_()
{
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
}

void MercurySolver::AddHallE_AtEdges_EnergyPreserving_()
{
    // ------------------------------------------------------------
    // Hall coefficient:
    //   E_hall = alpha * (J x B)
    //   alpha typically = +1/(n e)  (sign can be absorbed here)
    // ------------------------------------------------------------
    const double hall_coeff = hall_coef;
    const double N_floor = 1e-300; // 防止除零

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &UH = fld_->field(fid_.fid_U_H, iblk);
        auto &UNa = fld_->field(fid_.fid_U_Na, iblk);

        if (!UH.is_allocated() || !UNa.is_allocated())
            continue;

        auto &Bind_xi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Bind_eta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bind_zeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Badd_xi = fld_->field(fid_.fid_Badd.xi, iblk);
        auto &Badd_eta = fld_->field(fid_.fid_Badd.eta, iblk);
        auto &Badd_zeta = fld_->field(fid_.fid_Badd.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

        auto &E_xi = fld_->field(fid_.fid_E.xi, iblk);
        auto &E_eta = fld_->field(fid_.fid_E.eta, iblk);
        auto &E_zeta = fld_->field(fid_.fid_E.zeta, iblk);

        // edge cache: 9 comps (row-major 3x3)
        auto &pinvGT_xi = fld_->field(fid_.fid_pinvGT.xi, iblk);
        auto &pinvGT_eta = fld_->field(fid_.fid_pinvGT.eta, iblk);
        auto &pinvGT_zeta = fld_->field(fid_.fid_pinvGT.zeta, iblk);

        auto &pinvAT_xi = fld_->field(fid_.fid_pinvAT.xi, iblk);
        auto &pinvAT_eta = fld_->field(fid_.fid_pinvAT.eta, iblk);
        auto &pinvAT_zeta = fld_->field(fid_.fid_pinvAT.zeta, iblk);

        auto &x = grd_->grids(iblk).x;
        auto &y = grd_->grids(iblk).y;
        auto &z = grd_->grids(iblk).z;

        // small helpers
        auto matvec3 = [&](FieldBlock &M9, int i, int j, int k,
                           double c0, double c1, double c2) -> Vec3
        {
            Vec3 v;
            v.vec[0] = M9(i, j, k, 0) * c0 + M9(i, j, k, 1) * c1 + M9(i, j, k, 2) * c2;
            v.vec[1] = M9(i, j, k, 3) * c0 + M9(i, j, k, 4) * c1 + M9(i, j, k, 5) * c2;
            v.vec[2] = M9(i, j, k, 6) * c0 + M9(i, j, k, 7) * c1 + M9(i, j, k, 8) * c2;
            return v;
        };

        auto cache_is_valid = [&](FieldBlock &M9, int i, int j, int k) -> bool
        {
            double s = 0.0;
            for (int m = 0; m < 9; ++m)
                s += std::abs(M9(i, j, k, m));
            return s > 0.0;
        };

        auto Bxi = [&](int i, int j, int k, int m) -> double
        {
            return Bind_xi(i, j, k, m) + Badd_xi(i, j, k, m);
        };
        auto Beta = [&](int i, int j, int k, int m) -> double
        {
            return Bind_eta(i, j, k, m) + Badd_eta(i, j, k, m);
        };
        auto Bzeta = [&](int i, int j, int k, int m) -> double
        {
            return Bind_zeta(i, j, k, m) + Badd_zeta(i, j, k, m);
        };

        auto NUM = [&](int i, int j, int k, int m) -> double
        {
            return UH(i, j, k, 0) * rho_ref / M_H + UNa(i, j, k, 0) * rho_ref / M_Na;
        };

        // ============================================================
        // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
        // ============================================================
        {
            Int3 lo = E_xi.inner_lo();
            Int3 hi = E_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_xi, i, j, k) || !cache_is_valid(pinvAT_xi, i, j, k))
                        {
                            E_xi(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at xi-edge: average 4 surrounding cells (j- and k- directions)
                        double Num = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i, j - 1, k, 0) +
                                             NUM(i, j, k - 1, 0) +
                                             NUM(i, j - 1, k - 1, 0));
                        double alpha = hall_coeff / (Num + N_floor);

                        // Phi (2-form) co-located at xi-edge center
                        double Phi_xi = 0.0;
                        for (int di : {0, 1})
                            for (int dj : {0, -1})
                                for (int dk : {0, -1})
                                    Phi_xi += (Bxi(i + di, j + dj, k + dk, 0));
                        Phi_xi *= 0.125;

                        double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i, j, k - 1, 0));
                        double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i, j - 1, k, 0));

                        // j (1-form) co-located at same xi-edge center
                        double j_xi = Jxi(i, j, k, 0);
                        double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i + 1, j, k, 0) + Jeta(i, j - 1, k, 0) + Jeta(i + 1, j - 1, k, 0));
                        double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i + 1, j, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i + 1, j, k - 1, 0));

                        // map to physical vectors
                        Vec3 Jvec = matvec3(pinvGT_xi, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_xi, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i + 1, j, k) - x(i, j, k);
                        dr.vec[1] = y(i + 1, j, k) - y(i, j, k);
                        dr.vec[2] = z(i + 1, j, k) - z(i, j, k);

                        E_xi(i, j, k, 0) += Evec * dr; // line integral
                    }
        }

        // ============================================================
        // 2) EdgeEt
        // ============================================================
        {
            Int3 lo = E_eta.inner_lo();
            Int3 hi = E_eta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_eta, i, j, k) || !cache_is_valid(pinvAT_eta, i, j, k))
                        {
                            E_eta(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at eta-edge: average 4 surrounding cells (i- and k- directions)
                        double rho = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i - 1, j, k, 0) +
                                             NUM(i, j, k - 1, 0) +
                                             NUM(i - 1, j, k - 1, 0));
                        double alpha = hall_coeff / (rho + N_floor);

                        // Phi co-located at eta-edge center
                        double Phi_eta = 0.0;
                        for (int di : {0, -1})
                            for (int dj : {0, 1})
                                for (int dk : {0, -1})
                                    Phi_eta += Beta(i + di, j + dj, k + dk, 0);
                        Phi_eta *= 0.125;

                        double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j, k - 1, 0));
                        double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i - 1, j, k, 0));

                        // j co-located at eta-edge center
                        double j_eta = Jeta(i, j, k, 0);
                        double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j + 1, k, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j + 1, k, 0));
                        double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i, j + 1, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i, j + 1, k - 1, 0));

                        Vec3 Jvec = matvec3(pinvGT_eta, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_eta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i, j + 1, k) - x(i, j, k);
                        dr.vec[1] = y(i, j + 1, k) - y(i, j, k);
                        dr.vec[2] = z(i, j + 1, k) - z(i, j, k);

                        E_eta(i, j, k, 0) += Evec * dr;
                    }
        }

        // ============================================================
        // 3) EdgeZe
        // ============================================================
        {
            Int3 lo = E_zeta.inner_lo();
            Int3 hi = E_zeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_zeta, i, j, k) || !cache_is_valid(pinvAT_zeta, i, j, k))
                        {
                            E_zeta(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at zeta-edge: average 4 surrounding cells (i- and j- directions)
                        double rho = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i - 1, j, k, 0) +
                                             NUM(i, j - 1, k, 0) +
                                             NUM(i - 1, j - 1, k, 0));
                        double alpha = hall_coeff / (rho + N_floor);

                        // Phi co-located at zeta-edge center
                        double Phi_zeta = 0.0;
                        for (int di : {0, -1})
                            for (int dj : {0, -1})
                                for (int dk : {0, 1})
                                    Phi_zeta += Bzeta(i + di, j + dj, k + dk, 0);
                        Phi_zeta *= 0.125;

                        double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j - 1, k, 0));
                        double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i - 1, j, k, 0));

                        // j co-located at zeta-edge center
                        double j_zeta = Jzeta(i, j, k, 0);
                        double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j, k + 1, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j, k + 1, 0));
                        double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i, j, k + 1, 0) + Jeta(i, j - 1, k, 0) + Jeta(i, j - 1, k + 1, 0));

                        Vec3 Jvec = matvec3(pinvGT_zeta, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_zeta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i, j, k + 1) - x(i, j, k);
                        dr.vec[1] = y(i, j, k + 1) - y(i, j, k);
                        dr.vec[2] = z(i, j, k + 1) - z(i, j, k);

                        E_zeta(i, j, k, 0) += Evec * dr;
                    }
        }
    }
}