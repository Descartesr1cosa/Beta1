#include "MercurySolver.h"

void MercurySolver::AddHallEdgeEMF_()
{
    // Add Ehall_xi/eta/zeta（integration along Edge）to E
    // ------------------------------------------------------------
    // Hall coefficient:
    //   E_hall = alpha * (J x B)
    //   alpha typically = +1/(n e)  (sign can be absorbed here)
    // ------------------------------------------------------------
    const double hall_coeff = hall_coef;
    // const double N_floor = ne_hall_floor; // 1e-300; // 防止除零
    constexpr double C_eta = 0.0;

    const double C_alpha = 0.05; // 0.5; // 先试 0.3 ~ 0.8
    const double limit_C_alpha_eps = 1e-30;

    auto limit_alpha = [&](double alpha_phy, double Babs, double h2) -> double
    {
        if (h2 <= 0.0)
            return 0.0;

        // B 很小时，不限制
        if (Babs <= 1e-14)
            return alpha_phy;

        const double alpha_max = C_alpha * h2 / (Babs * dt_sub + limit_C_alpha_eps);

        const double amag = std::min(std::abs(alpha_phy), alpha_max);
        return std::copysign(amag, alpha_phy);
    };

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Ehall_et = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Ehall_ze = fld_->field(fid_.fid_Ehall.zeta, iblk);
        if (!Ehall_xi.is_allocated() || !Ehall_et.is_allocated() || !Ehall_ze.is_allocated())
            continue;

        {
            Int3 lo = Ehall_xi.inner_lo();
            Int3 hi = Ehall_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ehall_xi(i, j, k, 0) = 0.0;
        }
        {
            Int3 lo = Ehall_et.inner_lo();
            Int3 hi = Ehall_et.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ehall_et(i, j, k, 0) = 0.0;
        }
        {
            Int3 lo = Ehall_ze.inner_lo();
            Int3 hi = Ehall_ze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ehall_ze(i, j, k, 0) = 0.0;
        }
    }

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

        auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);

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

        auto &dlx = fld_->field("dl_xi", iblk);
        auto &dle = fld_->field("dl_eta", iblk);
        auto &dlz = fld_->field("dl_zeta", iblk);

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
            // return UH(i, j, k, 0) * rho_ref / M_H + UNa(i, j, k, 0) * rho_ref / M_Na;
            return UH(i, j, k, 0) / M_H + UNa(i, j, k, 0) / M_Na;
        };

        const double h_eps = 1e-12;
        auto hmin2 = [&](int i, int j, int k) -> double
        {
            double hx = dlx.is_allocated() ? dlx(i, j, k, 0) : 1e100;
            double he = dle.is_allocated() ? dle(i, j, k, 0) : 1e100;
            double hz = dlz.is_allocated() ? dlz(i, j, k, 0) : 1e100;
            double h = std::min(hx, std::min(he, hz));
            if (h <= h_eps)
                return 0.0;
            return h * h;
        };
        // ============================================================
        // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
        // ============================================================
        {
            Int3 lo = Ehall_xi.inner_lo();
            Int3 hi = Ehall_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_xi, i, j, k) || !cache_is_valid(pinvAT_xi, i, j, k))
                        {
                            Ehall_xi(i, j, k, 0) = 0.0;
                            continue;
                        }

                        const double h2 = hmin2(i, j, k);
                        if (h2 <= 0.0)
                            continue;

                        // rho at xi-edge: average 4 surrounding cells (j- and k- directions)
                        double Num = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i, j - 1, k, 0) +
                                             NUM(i, j, k - 1, 0) +
                                             NUM(i, j - 1, k - 1, 0));
                        // if (Num < ne_cut_hall)
                        // {
                        //     Ehall_xi(i, j, k, 0) = 0.0;
                        //     continue;
                        // } // EdgeXi

                        const double ne_true = std::max(Num, 0.0);
                        // 1) 平滑 floor（避免 max 的硬拐点）
                        const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                        // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                        const double s = ne_true / (ne_true + ne_hall_cut);
                        // 3) 统一有效 alpha（Hall 强度）
                        double alpha_phy = hall_coeff * s / ne_eff;

                        // // 可选：如果你仍想“极低密度完全关掉”，把 s 再乘一个 smoothstep 或者加个很小阈值
                        // double alpha = hall_coeff / (Num + N_floor);

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

                        // -------- limiter  --------
                        const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] +
                                                      Bvec.vec[1] * Bvec.vec[1] +
                                                      Bvec.vec[2] * Bvec.vec[2]);
                        const double alpha = limit_alpha(alpha_phy, Babs, h2);

                        Evec *= alpha;

                        // NEW: Hall stabilization (Ohmic / hyper-resistive in spirit)
                        // const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] + Bvec.vec[1] * Bvec.vec[1] + Bvec.vec[2] * Bvec.vec[2]);
                        const double eta_h = C_eta * std::abs(alpha) * Babs; // C_eta ~ 0.05 ~ 0.5
                        Evec.vec[0] += eta_h * Jvec.vec[0];
                        Evec.vec[1] += eta_h * Jvec.vec[1];
                        Evec.vec[2] += eta_h * Jvec.vec[2];

                        Vec3 dr;
                        dr.vec[0] = x(i + 1, j, k) - x(i, j, k);
                        dr.vec[1] = y(i + 1, j, k) - y(i, j, k);
                        dr.vec[2] = z(i + 1, j, k) - z(i, j, k);

                        Ehall_xi(i, j, k, 0) = Evec * dr; // line integral
                    }
        }

        // ============================================================
        // 2) EdgeEt
        // ============================================================
        {
            Int3 lo = Ehall_eta.inner_lo();
            Int3 hi = Ehall_eta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_eta, i, j, k) || !cache_is_valid(pinvAT_eta, i, j, k))
                        {
                            Ehall_eta(i, j, k, 0) = 0.0;
                            continue;
                        }

                        const double h2 = hmin2(i, j, k);
                        if (h2 <= 0.0)
                            continue;

                        // rho at eta-edge: average 4 surrounding cells (i- and k- directions)
                        double Num = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i - 1, j, k, 0) +
                                             NUM(i, j, k - 1, 0) +
                                             NUM(i - 1, j, k - 1, 0));
                        // if (Num < ne_cut_hall)
                        // {
                        //     Ehall_eta(i, j, k, 0) = 0.0;
                        //     continue;
                        // } // EdgeEt
                        // double alpha = hall_coeff / (Num + N_floor);

                        const double ne_true = std::max(Num, 0.0);
                        // 1) 平滑 floor（避免 max 的硬拐点）
                        const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                        // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                        const double s = ne_true / (ne_true + ne_hall_cut);
                        // 3) 统一有效 alpha（Hall 强度）
                        double alpha_phy = hall_coeff * s / ne_eff;

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

                        // -------- limiter  --------
                        const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] +
                                                      Bvec.vec[1] * Bvec.vec[1] +
                                                      Bvec.vec[2] * Bvec.vec[2]);
                        const double alpha = limit_alpha(alpha_phy, Babs, h2);

                        Evec *= alpha;

                        // NEW: Hall stabilization (Ohmic / hyper-resistive in spirit)
                        // const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] + Bvec.vec[1] * Bvec.vec[1] + Bvec.vec[2] * Bvec.vec[2]);
                        const double eta_h = C_eta * std::abs(alpha) * Babs; // C_eta ~ 0.05 ~ 0.5
                        Evec.vec[0] += eta_h * Jvec.vec[0];
                        Evec.vec[1] += eta_h * Jvec.vec[1];
                        Evec.vec[2] += eta_h * Jvec.vec[2];

                        Vec3 dr;
                        dr.vec[0] = x(i, j + 1, k) - x(i, j, k);
                        dr.vec[1] = y(i, j + 1, k) - y(i, j, k);
                        dr.vec[2] = z(i, j + 1, k) - z(i, j, k);

                        Ehall_eta(i, j, k, 0) = Evec * dr;
                    }
        }

        // ============================================================
        // 3) EdgeZe
        // ============================================================
        {
            Int3 lo = Ehall_zeta.inner_lo();
            Int3 hi = Ehall_zeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_zeta, i, j, k) || !cache_is_valid(pinvAT_zeta, i, j, k))
                        {
                            Ehall_zeta(i, j, k, 0) += 0.0;
                            continue;
                        }

                        const double h2 = hmin2(i, j, k);
                        if (h2 <= 0.0)
                            continue;

                        // rho at zeta-edge: average 4 surrounding cells (i- and j- directions)
                        double Num = 0.25 * (NUM(i, j, k, 0) +
                                             NUM(i - 1, j, k, 0) +
                                             NUM(i, j - 1, k, 0) +
                                             NUM(i - 1, j - 1, k, 0));
                        // if (Num < ne_cut_hall)
                        // {
                        //     Ehall_zeta(i, j, k, 0) = 0.0;
                        //     continue;
                        // } // EdgeZe
                        // double alpha = hall_coeff / (Num + N_floor);

                        const double ne_true = std::max(Num, 0.0);
                        // 1) 平滑 floor（避免 max 的硬拐点）
                        const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                        // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                        const double s = ne_true / (ne_true + ne_hall_cut);
                        // 3) 统一有效 alpha（Hall 强度）
                        double alpha_phy = hall_coeff * s / ne_eff;

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

                        // -------- limiter  --------
                        const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] +
                                                      Bvec.vec[1] * Bvec.vec[1] +
                                                      Bvec.vec[2] * Bvec.vec[2]);
                        const double alpha = limit_alpha(alpha_phy, Babs, h2);

                        Evec *= alpha;

                        // NEW: Hall stabilization (Ohmic / hyper-resistive in spirit)
                        // const double Babs = std::sqrt(Bvec.vec[0] * Bvec.vec[0] + Bvec.vec[1] * Bvec.vec[1] + Bvec.vec[2] * Bvec.vec[2]);
                        const double eta_h = C_eta * std::abs(alpha) * Babs; // C_eta ~ 0.05 ~ 0.5
                        Evec.vec[0] += eta_h * Jvec.vec[0];
                        Evec.vec[1] += eta_h * Jvec.vec[1];
                        Evec.vec[2] += eta_h * Jvec.vec[2];

                        Vec3 dr;
                        dr.vec[0] = x(i, j, k + 1) - x(i, j, k);
                        dr.vec[1] = y(i, j, k + 1) - y(i, j, k);
                        dr.vec[2] = z(i, j, k + 1) - z(i, j, k);

                        Ehall_zeta(i, j, k, 0) = Evec * dr;
                    }
        }
    }

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);

        auto &x = grd_->grids(iblk).x;
        auto &y = grd_->grids(iblk).y;
        auto &z = grd_->grids(iblk).z;

        auto radius = [&](int i, int j, int k) -> double
        {
            double xx = x(i, j, k);
            double yy = y(i, j, k);
            double zz = z(i, j, k);
            return std::sqrt(xx * xx + yy * yy + zz * zz);
        };

        auto hall_factor_s = [&](double r) -> double
        {
            const double r0 = 1.1; // 内边界
            const double r1 = 1.2; // taper 外边界

            if (r <= r0)
                return 0.0;
            if (r >= r1)
                return 1.0;

            double xi = (r - r0) / (r1 - r0);  // 映射到 [0,1]
            return xi * xi * (3.0 - 2.0 * xi); // smoothstep
        };

        // ============================================================
        // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
        // ============================================================
        {
            Int3 lo = Ehall_xi.inner_lo();
            Int3 hi = Ehall_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double fac = hall_factor_s(radius(i, j, k));
                        Ehall_xi(i, j, k, 0) *= fac;
                    }
        }

        // ============================================================
        // 2) EdgeEt
        // ============================================================
        {
            Int3 lo = Ehall_eta.inner_lo();
            Int3 hi = Ehall_eta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double fac = hall_factor_s(radius(i, j, k));
                        Ehall_eta(i, j, k, 0) *= fac;
                    }
        }

        // ============================================================
        // 3) EdgeZe
        // ============================================================
        {
            Int3 lo = Ehall_zeta.inner_lo();
            Int3 hi = Ehall_zeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double fac = hall_factor_s(radius(i, j, k));
                        Ehall_zeta(i, j, k, 0) *= fac;
                    }
        }
    }

    mercury_bound_.Sync("Ehall");

    // mercury_bound_.DebugCheckSurfaceTangentialEdgeField(
    //     "Ehall-inner-solid", "Ehall_xi", "Ehall_eta", "Ehall_zeta", "Coupled-Solid");

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);
        if (!Ehall_xi.is_allocated() || !Ehall_eta.is_allocated() || !Ehall_zeta.is_allocated())
            continue;

        auto &E_xi = fld_->field(fid_.fid_E.xi, iblk);
        auto &E_eta = fld_->field(fid_.fid_E.eta, iblk);
        auto &E_zeta = fld_->field(fid_.fid_E.zeta, iblk);

        if (!E_xi.is_allocated() || !E_eta.is_allocated() || !E_zeta.is_allocated())
            continue;

        // ============================================================
        // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
        // ============================================================
        {
            Int3 lo = E_xi.inner_lo();
            Int3 hi = E_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_xi(i, j, k, 0) += Ehall_xi(i, j, k, 0);
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
                        E_eta(i, j, k, 0) += Ehall_eta(i, j, k, 0);
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
                        E_zeta(i, j, k, 0) += Ehall_zeta(i, j, k, 0);
        }
    }
}
