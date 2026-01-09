
#include "MercurySolver.h"

void MercurySolver::Scheme_B_()
{
    auto minmod = [](double a, double b) -> double
    {
        if (a * b <= 0.0)
            return 0.0;
        return (std::abs(a) < std::abs(b)) ? a : b;
    };

    auto flux_B = [](const std::array<double, 6> &u, double ex, double ey, double ez,
                     std::array<double, 3> &fs)
    {
        const double Ue = ex * u[3] + ey * u[4] + ez * u[5];
        const double Be = ex * u[0] + ey * u[1] + ez * u[2];
        fs[0] = Ue * u[0] - u[3] * Be;
        fs[1] = Ue * u[1] - u[4] * Be;
        fs[2] = Ue * u[2] - u[5] * Be;
    };

    auto muscl_recon_B = [&](const std::array<double, 6> &uL2,
                             const std::array<double, 6> &uL1,
                             const std::array<double, 6> &uR1,
                             const std::array<double, 6> &uR2,
                             std::array<double, 6> &UL,
                             std::array<double, 6> &UR)
    {
        // 只重构 B(1:3)，速度(4:6)取 cell-center（与 Fortran 一致）
        for (int m = 0; m < 3; ++m)
        {
            const double sL = minmod(uL1[m] - uL2[m], uR1[m] - uL1[m]);
            const double sR = minmod(uR1[m] - uL1[m], uR2[m] - uR1[m]);
            UL[m] = uL1[m] + 0.5 * sL;
            UR[m] = uR1[m] - 0.5 * sR;
        }
        UL[3] = uL1[3];
        UL[4] = uL1[4];
        UL[5] = uL1[5];
        UR[3] = uR1[3];
        UR[4] = uR1[4];
        UR[5] = uR1[5];
    };

    auto fluxRusanov_MUSCL_B = [&](const std::array<double, 6> &uL2,
                                   const std::array<double, 6> &uL1,
                                   const std::array<double, 6> &uR1,
                                   const std::array<double, 6> &uR2,
                                   double ex, double ey, double ez,
                                   std::array<double, 3> &fn)
    {
        std::array<double, 6> UL{}, UR{};
        muscl_recon_B(uL2, uL1, uR1, uR2, UL, UR);

        std::array<double, 3> fL{}, fR{};
        flux_B(UL, ex, ey, ez, fL);
        flux_B(UR, ex, ey, ez, fR);

        const double UeL = ex * UL[3] + ey * UL[4] + ez * UL[5];
        const double UeR = ex * UR[3] + ey * UR[4] + ez * UR[5];
        const double amax = std::max(std::abs(UeL), std::abs(UeR));

        for (int l = 0; l < 3; ++l)
            fn[l] = 0.5 * (fL[l] + fR[l]) - 0.5 * amax * (UR[l] - UL[l]);
    };

    auto scheme_B_point = [&](const std::array<double, 6> &u0,
                              const std::array<double, 6> &u1,
                              const std::array<double, 6> &u2,
                              const std::array<double, 6> &u3,
                              const std::array<double, 6> &u4,
                              double ex, double ey, double ez,
                              std::array<double, 3> &ffs)
    {
        // i-1/2: (u0,u1,u2,u3)
        std::array<double, 3> fn12{}, fn23{};
        fluxRusanov_MUSCL_B(u0, u1, u2, u3, ex, ey, ez, fn12);
        // i+1/2: (u1,u2,u3,u4)
        fluxRusanov_MUSCL_B(u1, u2, u3, u4, ex, ey, ez, fn23);

        // Fortran 的 BM 修正项（保持一致）
        const double Be1 = ex * u1[0] + ey * u1[1] + ez * u1[2];
        const double Be3 = ex * u3[0] + ey * u3[1] + ez * u3[2];
        const double Bee = 0.5 * (Be3 - Be1);

        const double BM0 = u2[3] * Bee;
        const double BM1 = u2[4] * Bee;
        const double BM2 = u2[5] * Bee;

        ffs[0] = -(fn23[0] - fn12[0]) - BM0;
        ffs[1] = -(fn23[1] - fn12[1]) - BM1;
        ffs[2] = -(fn23[2] - fn12[2]) - BM2;
    };

    const double tinyV = 1e-30;

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);

        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);   // JDxi (FaceXi,3)
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);  // JDet (FaceEt,3)
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib); // JDze (FaceZe,3)

        FieldBlock &B = fld_->field(fid_.fid_Bcell, ib);   // B_cell (Cell,3,ngg)
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib); // U_plus (Cell,3,ngg)
        FieldBlock &RHS = fld_->field(fid_.fid_RHS_b, ib); // RHS_B  (Cell,3,0)

        if (!Jac.is_allocated() || !Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated() ||
            !B.is_allocated() || !Up.is_allocated() || !RHS.is_allocated())
            continue;

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        auto load_state = [&](int i, int j, int k, std::array<double, 6> &u)
        {
            u[0] = B(i, j, k, 0);
            u[1] = B(i, j, k, 1);
            u[2] = B(i, j, k, 2);
            u[3] = Up(i, j, k, 0);
            u[4] = Up(i, j, k, 1);
            u[5] = Up(i, j, k, 2);
        };

        auto metric_grad = [&](int dir, int i, int j, int k, double &ex, double &ey, double &ez)
        {
            const double V = std::abs(Jac(i, j, k, 0));
            if (V < tinyV)
            {
                ex = ey = ez = 0.0;
                return;
            }

            if (dir == 0) // xi: avg(JDxi(i), JDxi(i+1))/V
            {
                ex = 0.5 * (Axi(i, j, k, 0) + Axi(i - 1, j, k, 0)) / V;
                ey = 0.5 * (Axi(i, j, k, 1) + Axi(i - 1, j, k, 1)) / V;
                ez = 0.5 * (Axi(i, j, k, 2) + Axi(i - 1, j, k, 2)) / V;
            }
            else if (dir == 1) // eta: avg(JDet(j), JDet(j+1))/V
            {
                ex = 0.5 * (Aet(i, j, k, 0) + Aet(i, j - 1, k, 0)) / V;
                ey = 0.5 * (Aet(i, j, k, 1) + Aet(i, j - 1, k, 1)) / V;
                ez = 0.5 * (Aet(i, j, k, 2) + Aet(i, j - 1, k, 2)) / V;
            }
            else // dir == 2, zeta: avg(JDze(k), JDze(k+1))/V
            {
                ex = 0.5 * (Aze(i, j, k, 0) + Aze(i, j, k - 1, 0)) / V;
                ey = 0.5 * (Aze(i, j, k, 1) + Aze(i, j, k - 1, 1)) / V;
                ez = 0.5 * (Aze(i, j, k, 2) + Aze(i, j, k - 1, 2)) / V;
            }
        };

        // 三个方向 sweep：对应 Fortran tvd_B -> solutionspace_B(te,le=1..3)
        for (int dir = 0; dir < 3; ++dir)
        {
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // 5 点 stencil（要求 ngg >= 2）
                        std::array<double, 6> u0{}, u1{}, u2{}, u3{}, u4{};
                        if (dir == 0)
                        {
                            load_state(i - 2, j, k, u0);
                            load_state(i - 1, j, k, u1);
                            load_state(i, j, k, u2);
                            load_state(i + 1, j, k, u3);
                            load_state(i + 2, j, k, u4);
                        }
                        else if (dir == 1)
                        {
                            load_state(i, j - 2, k, u0);
                            load_state(i, j - 1, k, u1);
                            load_state(i, j, k, u2);
                            load_state(i, j + 1, k, u3);
                            load_state(i, j + 2, k, u4);
                        }
                        else
                        {
                            load_state(i, j, k - 2, u0);
                            load_state(i, j, k - 1, u1);
                            load_state(i, j, k, u2);
                            load_state(i, j, k + 1, u3);
                            load_state(i, j, k + 2, u4);
                        }

                        double ex, ey, ez;
                        metric_grad(dir, i, j, k, ex, ey, ez);

                        std::array<double, 3> ffs{};
                        scheme_B_point(u0, u1, u2, u3, u4, ex, ey, ez, ffs);

                        RHS(i, j, k, 0) += ffs[0];
                        RHS(i, j, k, 1) += ffs[1];
                        RHS(i, j, k, 2) += ffs[2];
                    }
        }
    }
}