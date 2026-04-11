#include "MercurySolver.h"

void MercurySolver::AddResistiveEdgeEMF_()
{
    const int nb = fld_->num_blocks();

    const double r_cut_in = 0.8;
    const double r_cut_out = 0.98;
    const double r0 = 0.8;
    const double r1 = 0.98;
    const double w = 0.01;

    auto yita0_of_r = [&](double r) -> double
    {
        if (r <= r_cut_in || r >= r_cut_out)
            return 0.0;
        return 0.5 * (std::tanh((r - r0) / w) - std::tanh((r - r1) / w));
    };

    // ----  add resistive E on solid blocks: E += invRem8 * yita0_edge * J ----
    for (int ib = 0; ib < nb; ++ib)
    {
        Block &blk = fld_->grd->grids(ib);
        if (blk.block_name != "Solid")
            continue;

        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        // --- EdgeXi: use the edge midpoint of node (i,j,k) -> (i+1,j,k) ---
        {
            Int3 lo = Exi.inner_lo();
            Int3 hi = Exi.inner_hi();

            // 如果你后续在扩散项里不需要 i±1/j±1/k±1，就不必“严格 inner 收缩”；
            // 这里我们只用本点的 Jxi 和本边的坐标，所以不用收缩。
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i + 1, j, k));
                        const double ym = 0.5 * (y(i, j, k) + y(i + 1, j, k));
                        const double zm = 0.5 * (z(i, j, k) + z(i + 1, j, k));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Exi(i, j, k, 0) += (inver_Rem * yita0) * Jxi(i, j, k, 0);
                    }
        }

        // --- EdgeEt: node (i,j,k) -> (i,j+1,k) ---
        {
            Int3 lo = Eeta.inner_lo();
            Int3 hi = Eeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i, j + 1, k));
                        const double ym = 0.5 * (y(i, j, k) + y(i, j + 1, k));
                        const double zm = 0.5 * (z(i, j, k) + z(i, j + 1, k));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Eeta(i, j, k, 0) += (inver_Rem * yita0) * Jeta(i, j, k, 0);
                    }
        }

        // --- EdgeZe: node (i,j,k) -> (i,j,k+1) ---
        {
            Int3 lo = Eze.inner_lo();
            Int3 hi = Eze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i, j, k + 1));
                        const double ym = 0.5 * (y(i, j, k) + y(i, j, k + 1));
                        const double zm = 0.5 * (z(i, j, k) + z(i, j, k + 1));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Eze(i, j, k, 0) += (inver_Rem * yita0) * Jzeta(i, j, k, 0);
                    }
        }
    }
}

void MercurySolver::AddSecondResistiveEdgeEMF_()
{
    constexpr double eps = 1e-14;

    // 只读一个参数；建议先从 0.02 ~ 0.06 试
    const double c2 = par_->GetDou("HallResist2");

    // 下面这些都先写死，避免参数过多
    // chi = h^2 |lap(J)| / (|J| + Jref + eps)
    // 当 chi 小时不开；chi 大时逐渐打开
    // constexpr double chi0 = 0.15;
    // constexpr double chi1 = 0.45;
    constexpr double chi0 = 0.05;
    constexpr double chi1 = 0.2;

    auto clamp01 = [](double x) -> double
    {
        return std::max(0.0, std::min(1.0, x));
    };

    auto smoothstep = [&](double x0, double x1, double x) -> double
    {
        if (x1 <= x0)
            return (x > x0 ? 1.0 : 0.0);
        const double t = clamp01((x - x0) / (x1 - x0));
        return t * t * (3.0 - 2.0 * t);
    };

    auto max4 = [](double a, double b, double c, double d) -> double
    {
        return std::max(std::max(a, b), std::max(c, d));
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eet = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jet = fld_->field(fid_.fid_J.eta, ib);
        auto &Jze = fld_->field(fid_.fid_J.zeta, ib);

        auto &dl_xi = fld_->field(fid_.Edge_dl.xi, ib);
        auto &dl_et = fld_->field(fid_.Edge_dl.eta, ib);
        auto &dl_ze = fld_->field(fid_.Edge_dl.zeta, ib);

        auto &beta = hall_face_scratch_[ib].beta; // 由 BuildHallFaceEMF_Rusanov_() 填好

        if (!Exi.is_allocated() || !Jxi.is_allocated())
            continue;

        // ============================================================
        // xi-edge
        // ============================================================
        {
            Int3 lo = Exi.inner_lo();
            Int3 hi = Exi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i, j - 1, k - 1),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double hx = std::max(dl_xi(i, j, k, 0), eps);
                        const double hy = std::max(0.5 * (dl_et(i, j, k, 0) + dl_et(i + 1, j, k, 0)), eps);
                        const double hz = std::max(0.5 * (dl_ze(i, j, k, 0) + dl_ze(i + 1, j, k, 0)), eps);

                        const double h = std::max(std::min({hx, hy, hz}), eps);
                        const double h2 = h * h;

                        const double Jc = Jxi(i, j, k, 0);

                        const double lapJ =
                            (Jxi(i + 1, j, k, 0) - 2.0 * Jc + Jxi(i - 1, j, k, 0)) / h2 +
                            (Jxi(i, j + 1, k, 0) - 2.0 * Jc + Jxi(i, j - 1, k, 0)) / h2 +
                            (Jxi(i, j, k + 1, 0) - 2.0 * Jc + Jxi(i, j, k - 1, 0)) / h2;

                        // 很便宜的局部参考尺度，避免 |J| 很小时误触发
                        // const double Jref = 0.05 * beta_e + 1e-12;
                        const double Jref = 0.2 * beta_e + 1e-12;

                        const double chi = h2 * std::abs(lapJ) / (std::abs(Jc) + Jref + eps);
                        const double S = smoothstep(chi0, chi1, chi);

                        // const double eta2 = c2 * beta_e * h * S;
                        const double S_floor = 0.05;
                        const double eta2 = c2 * beta_e * h * (S_floor + (1.0 - S_floor) * S);

                        // E += eta J  -->  CT 下对应二阶磁扩散
                        Exi(i, j, k, 0) += eta2 * Jc;
                    }
        }

        // ============================================================
        // eta-edge
        // ============================================================
        {
            Int3 lo = Eet.inner_lo();
            Int3 hi = Eet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i - 1, j, k - 1),
                                                   beta(i - 1, j, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double hx = std::max(0.5 * (dl_xi(i, j, k, 0) + dl_xi(i, j + 1, k, 0)), eps);
                        const double hy = std::max(dl_et(i, j, k, 0), eps);
                        const double hz = std::max(0.5 * (dl_ze(i, j, k, 0) + dl_ze(i, j + 1, k, 0)), eps);

                        const double h = std::max(std::min({hx, hy, hz}), eps);
                        const double h2 = h * h;

                        const double Jc = Jet(i, j, k, 0);

                        const double lapJ =
                            (Jet(i + 1, j, k, 0) - 2.0 * Jc + Jet(i - 1, j, k, 0)) / h2 +
                            (Jet(i, j + 1, k, 0) - 2.0 * Jc + Jet(i, j - 1, k, 0)) / h2 +
                            (Jet(i, j, k + 1, 0) - 2.0 * Jc + Jet(i, j, k - 1, 0)) / h2;

                        // const double Jref = 0.05 * beta_e + 1e-12;
                        const double Jref = 0.2 * beta_e + 1e-12;

                        const double chi = h2 * std::abs(lapJ) / (std::abs(Jc) + Jref + eps);
                        const double S = smoothstep(chi0, chi1, chi);

                        // const double eta2 = c2 * beta_e * h * S;
                        const double S_floor = 0.05;
                        const double eta2 = c2 * beta_e * h * (S_floor + (1.0 - S_floor) * S);

                        Eet(i, j, k, 0) += eta2 * Jc;
                    }
        }

        // ============================================================
        // zeta-edge
        // ============================================================
        {
            Int3 lo = Eze.inner_lo();
            Int3 hi = Eze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i - 1, j - 1, k),
                                                   beta(i - 1, j, k),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k));

                        const double hx = std::max(0.5 * (dl_xi(i, j, k, 0) + dl_xi(i, j, k + 1, 0)), eps);
                        const double hy = std::max(0.5 * (dl_et(i, j, k, 0) + dl_et(i, j, k + 1, 0)), eps);
                        const double hz = std::max(dl_ze(i, j, k, 0), eps);

                        const double h = std::max(std::min({hx, hy, hz}), eps);
                        const double h2 = h * h;

                        const double Jc = Jze(i, j, k, 0);

                        const double lapJ =
                            (Jze(i + 1, j, k, 0) - 2.0 * Jc + Jze(i - 1, j, k, 0)) / h2 +
                            (Jze(i, j + 1, k, 0) - 2.0 * Jc + Jze(i, j - 1, k, 0)) / h2 +
                            (Jze(i, j, k + 1, 0) - 2.0 * Jc + Jze(i, j, k - 1, 0)) / h2;

                        // const double Jref = 0.05 * beta_e + 1e-12;
                        const double Jref = 0.2 * beta_e + 1e-12;

                        const double chi = h2 * std::abs(lapJ) / (std::abs(Jc) + Jref + eps);
                        const double S = smoothstep(chi0, chi1, chi);

                        // const double eta2 = c2 * beta_e * h * S;
                        const double S_floor = 0.05;
                        const double eta2 = c2 * beta_e * h * (S_floor + (1.0 - S_floor) * S);

                        Eze(i, j, k, 0) += eta2 * Jc;
                    }
        }
    }
}

void MercurySolver::AddHyperResistiveEdgeEMF_()
{
    constexpr double eps = 1e-14;
    const double chyp = par_->GetDou("HyperResist"); // 0.02; // 先小一点，后面再调

    auto max4 = [](double a, double b, double c, double d) -> double
    {
        return std::max(std::max(a, b), std::max(c, d));
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eet = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jet = fld_->field(fid_.fid_J.eta, ib);
        auto &Jze = fld_->field(fid_.fid_J.zeta, ib);

        auto &dl_xi = fld_->field(fid_.Edge_dl.xi, ib);
        auto &dl_et = fld_->field(fid_.Edge_dl.eta, ib);
        auto &dl_ze = fld_->field(fid_.Edge_dl.zeta, ib);

        auto &beta = hall_face_scratch_[ib].beta; // 已由 BuildHallFaceEMF_Rusanov_() 填好

        if (!Exi.is_allocated() || !Jxi.is_allocated())
            continue;

        // ------------------------------------------------------------
        // xi-edge
        // ------------------------------------------------------------
        {
            Int3 lo = Exi.inner_lo();
            Int3 hi = Exi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i, j - 1, k - 1),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double hx = std::max(dl_xi(i, j, k, 0), eps);
                        const double hy = std::max(0.5 * (dl_et(i, j, k, 0) + dl_et(i + 1, j, k, 0)), eps);
                        const double hz = std::max(0.5 * (dl_ze(i, j, k, 0) + dl_ze(i + 1, j, k, 0)), eps);

                        const double h2 = std::max(std::min({hx * hx, hy * hy, hz * hz}), eps);

                        const double lap =
                            (Jxi(i + 1, j, k, 0) - 2.0 * Jxi(i, j, k, 0) + Jxi(i - 1, j, k, 0)) / h2 +
                            (Jxi(i, j + 1, k, 0) - 2.0 * Jxi(i, j, k, 0) + Jxi(i, j - 1, k, 0)) / h2 +
                            (Jxi(i, j, k + 1, 0) - 2.0 * Jxi(i, j, k, 0) + Jxi(i, j, k - 1, 0)) / h2;

                        const double nu4 = chyp * beta_e * h2;
                        Exi(i, j, k, 0) += -nu4 * lap;
                    }
        }

        // ------------------------------------------------------------
        // eta-edge
        // ------------------------------------------------------------
        {
            Int3 lo = Eet.inner_lo();
            Int3 hi = Eet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i - 1, j, k - 1),
                                                   beta(i - 1, j, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double hx = std::max(0.5 * (dl_xi(i, j, k, 0) + dl_xi(i, j + 1, k, 0)), eps);
                        const double hy = std::max(dl_et(i, j, k, 0), eps);
                        const double hz = std::max(0.5 * (dl_ze(i, j, k, 0) + dl_ze(i, j + 1, k, 0)), eps);

                        const double h2 = std::max(std::min({hx * hx, hy * hy, hz * hz}), eps);

                        const double lap =
                            (Jet(i + 1, j, k, 0) - 2.0 * Jet(i, j, k, 0) + Jet(i - 1, j, k, 0)) / h2 +
                            (Jet(i, j + 1, k, 0) - 2.0 * Jet(i, j, k, 0) + Jet(i, j - 1, k, 0)) / h2 +
                            (Jet(i, j, k + 1, 0) - 2.0 * Jet(i, j, k, 0) + Jet(i, j, k - 1, 0)) / h2;

                        const double nu4 = chyp * beta_e * h2;
                        Eet(i, j, k, 0) += -nu4 * lap;
                    }
        }

        // ------------------------------------------------------------
        // zeta-edge
        // ------------------------------------------------------------
        {
            Int3 lo = Eze.inner_lo();
            Int3 hi = Eze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double beta_e = max4(beta(i - 1, j - 1, k),
                                                   beta(i - 1, j, k),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k));

                        const double hx = std::max(0.5 * (dl_xi(i, j, k, 0) + dl_xi(i, j, k + 1, 0)), eps);
                        const double hy = std::max(0.5 * (dl_et(i, j, k, 0) + dl_et(i, j, k + 1, 0)), eps);
                        const double hz = std::max(dl_ze(i, j, k, 0), eps);

                        const double h2 = std::max(std::min({hx * hx, hy * hy, hz * hz}), eps);

                        const double lap =
                            (Jze(i + 1, j, k, 0) - 2.0 * Jze(i, j, k, 0) + Jze(i - 1, j, k, 0)) / h2 +
                            (Jze(i, j + 1, k, 0) - 2.0 * Jze(i, j, k, 0) + Jze(i, j - 1, k, 0)) / h2 +
                            (Jze(i, j, k + 1, 0) - 2.0 * Jze(i, j, k, 0) + Jze(i, j, k - 1, 0)) / h2;

                        const double nu4 = chyp * beta_e * h2;
                        Eze(i, j, k, 0) += -nu4 * lap;
                    }
        }
    }
}
