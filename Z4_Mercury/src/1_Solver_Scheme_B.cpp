
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

        FieldBlock &B = fld_->field(fid_.fid_U_b, ib);     // B_cell (Cell,3,ngg)
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

void MercurySolver::AddSourceToRHS_B()
{
    // ---------- coefficients (default 0 if not provided) ----------
    double a7 = 0.0, a8 = 0.0;
    // {
    //     auto &mp = par_->GetDou_List("constant").data;
    //     if (mp.count("a7"))
    //         a7 = mp["a7"];
    //     if (mp.count("a8"))
    //         a8 = mp["a8"];
    // }

    const double ne_floor = 1e-30;
    const double inv23 = 1.0 / 23.0; // Na 质量比（与 calc_Uplus 一致）

    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) + b * (fp_j - fm_j) + c * (fp_k - fm_k));
    };

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        FieldBlock &Ub = fld_->field(fid_.fid_U_b, ib);    // qb(1:3)
        FieldBlock &B0 = fld_->field(fid_.fid_Badd, ib);   // B0
        FieldBlock &Bt = fld_->field(fid_.fid_Bcell, ib);  // Bt=qb+B0
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib); // u+

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UNa = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib); // p in comp=3
        FieldBlock &PVN = fld_->field(fid_.fid_PV_Na, ib);

        FieldBlock &RHS = fld_->field(fid_.fid_RHS_b, ib); // db(1:3)

        if (!Jac.is_allocated() || !Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated())
            continue;
        if (!Ub.is_allocated() || !B0.is_allocated() || !Bt.is_allocated() || !Up.is_allocated())
            continue;
        if (!UH.is_allocated() || !UNa.is_allocated() || !PVH.is_allocated() || !PVN.is_allocated())
            continue;
        if (!RHS.is_allocated())
            continue;

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        auto derv_at = [&](int i, int j, int k,
                           double &ax, double &ay, double &az,
                           double &bx, double &by, double &bz,
                           double &cx, double &cy, double &cz)
        {
            const double V = std::abs(Jac(i, j, k, 0));
            if (V <= 0.0)
            {
                ax = ay = az = bx = by = bz = cx = cy = cz = 0.0;
                return;
            }

            // Axi(i,j,k) is i+1/2 face: cell i uses faces (i-1) and (i)
            ax = 0.5 * (Axi(i - 1, j, k, 0) + Axi(i, j, k, 0)) / V;
            ay = 0.5 * (Axi(i - 1, j, k, 1) + Axi(i, j, k, 1)) / V;
            az = 0.5 * (Axi(i - 1, j, k, 2) + Axi(i, j, k, 2)) / V;

            bx = 0.5 * (Aet(i, j - 1, k, 0) + Aet(i, j, k, 0)) / V;
            by = 0.5 * (Aet(i, j - 1, k, 1) + Aet(i, j, k, 1)) / V;
            bz = 0.5 * (Aet(i, j - 1, k, 2) + Aet(i, j, k, 2)) / V;

            cx = 0.5 * (Aze(i, j, k - 1, 0) + Aze(i, j, k, 0)) / V;
            cy = 0.5 * (Aze(i, j, k - 1, 1) + Aze(i, j, k, 1)) / V;
            cz = 0.5 * (Aze(i, j, k - 1, 2) + Aze(i, j, k, 2)) / V;
        };

        auto ne_at = [&](int i, int j, int k) -> double
        {
            const double rhoH = std::max(UH(i, j, k, 0), 0.0);
            const double rhoNa = std::max(UNa(i, j, k, 0), 0.0);
            // number density ~ rho / m : m_H=1, m_Na=23
            return rhoH + rhoNa * inv23;
        };

        auto pe_at = [&](int i, int j, int k) -> double
        {
            return PVH(i, j, k, 3) + PVN(i, j, k, 3);
        };

        // compute curl(qb) at a point -> J
        auto J_at = [&](int i, int j, int k, double J[3])
        {
            double ax, ay, az, bx, by, bz, cx, cy, cz;
            derv_at(i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

            const double Bx_p_i = Ub(i + 1, j, k, 0), Bx_m_i = Ub(i - 1, j, k, 0);
            const double Bx_p_j = Ub(i, j + 1, k, 0), Bx_m_j = Ub(i, j - 1, k, 0);
            const double Bx_p_k = Ub(i, j, k + 1, 0), Bx_m_k = Ub(i, j, k - 1, 0);

            const double By_p_i = Ub(i + 1, j, k, 1), By_m_i = Ub(i - 1, j, k, 1);
            const double By_p_j = Ub(i, j + 1, k, 1), By_m_j = Ub(i, j - 1, k, 1);
            const double By_p_k = Ub(i, j, k + 1, 1), By_m_k = Ub(i, j, k - 1, 1);

            const double Bz_p_i = Ub(i + 1, j, k, 2), Bz_m_i = Ub(i - 1, j, k, 2);
            const double Bz_p_j = Ub(i, j + 1, k, 2), Bz_m_j = Ub(i, j - 1, k, 2);
            const double Bz_p_k = Ub(i, j, k + 1, 2), Bz_m_k = Ub(i, j, k - 1, 2);

            // Fortran form:
            const double Bzy = dd(ay, by, cy, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);
            const double Byz = dd(az, bz, cz, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);

            const double Bxz = dd(az, bz, cz, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);
            const double Bzx = dd(ax, bx, cx, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);

            const double Byx = dd(ax, bx, cx, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);
            const double Bxy = dd(ay, by, cy, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);

            J[0] = Bzy - Byz;
            J[1] = Bxz - Bzx;
            J[2] = Byx - Bxy;
        };

        auto sjbne_at = [&](int i, int j, int k, double out[3])
        {
            const double ne = ne_at(i, j, k);
            if (ne < ne_floor)
            {
                out[0] = out[1] = out[2] = 0.0;
                return;
            }

            double J[3];
            J_at(i, j, k, J);

            const double Btx = Bt(i, j, k, 0);
            const double Bty = Bt(i, j, k, 1);
            const double Btz = Bt(i, j, k, 2);

            // sjb = J x Bt
            const double sjbx = J[1] * Btz - J[2] * Bty;
            const double sjby = J[2] * Btx - J[0] * Btz;
            const double sjbz = J[0] * Bty - J[1] * Btx;

            const double invne = 1.0 / ne;
            out[0] = sjbx * invne;
            out[1] = sjby * invne;
            out[2] = sjbz * invne;
        };

        auto upb0_comp = [&](int i, int j, int k, int comp) -> double
        {
            const double ux = Up(i, j, k, 0), uy = Up(i, j, k, 1), uz = Up(i, j, k, 2);
            const double Bx = B0(i, j, k, 0), By = B0(i, j, k, 1), Bz = B0(i, j, k, 2);
            // u x B0
            const double cx = uy * Bz - uz * By;
            const double cy = uz * Bx - ux * Bz;
            const double cz = ux * By - uy * Bx;
            return (comp == 0) ? cx : ((comp == 1) ? cy : cz);
        };

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    double ax, ay, az, bx, by, bz, cx, cy, cz;
                    derv_at(i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

                    // -------- term 1: curl(u+ x B0) --------
                    const double upb0_1_p_i = upb0_comp(i + 1, j, k, 0), upb0_1_m_i = upb0_comp(i - 1, j, k, 0);
                    const double upb0_1_p_j = upb0_comp(i, j + 1, k, 0), upb0_1_m_j = upb0_comp(i, j - 1, k, 0);
                    const double upb0_1_p_k = upb0_comp(i, j, k + 1, 0), upb0_1_m_k = upb0_comp(i, j, k - 1, 0);

                    const double upb0_2_p_i = upb0_comp(i + 1, j, k, 1), upb0_2_m_i = upb0_comp(i - 1, j, k, 1);
                    const double upb0_2_p_j = upb0_comp(i, j + 1, k, 1), upb0_2_m_j = upb0_comp(i, j - 1, k, 1);
                    const double upb0_2_p_k = upb0_comp(i, j, k + 1, 1), upb0_2_m_k = upb0_comp(i, j, k - 1, 1);

                    const double upb0_3_p_i = upb0_comp(i + 1, j, k, 2), upb0_3_m_i = upb0_comp(i - 1, j, k, 2);
                    const double upb0_3_p_j = upb0_comp(i, j + 1, k, 2), upb0_3_m_j = upb0_comp(i, j - 1, k, 2);
                    const double upb0_3_p_k = upb0_comp(i, j, k + 1, 2), upb0_3_m_k = upb0_comp(i, j, k - 1, 2);

                    const double ubzy = dd(ay, by, cy, upb0_3_p_i, upb0_3_m_i, upb0_3_p_j, upb0_3_m_j, upb0_3_p_k, upb0_3_m_k);
                    const double ubyz = dd(az, bz, cz, upb0_2_p_i, upb0_2_m_i, upb0_2_p_j, upb0_2_m_j, upb0_2_p_k, upb0_2_m_k);

                    const double ubxz = dd(az, bz, cz, upb0_1_p_i, upb0_1_m_i, upb0_1_p_j, upb0_1_m_j, upb0_1_p_k, upb0_1_m_k);
                    const double ubzx = dd(ax, bx, cx, upb0_3_p_i, upb0_3_m_i, upb0_3_p_j, upb0_3_m_j, upb0_3_p_k, upb0_3_m_k);

                    const double ubyx = dd(ax, bx, cx, upb0_2_p_i, upb0_2_m_i, upb0_2_p_j, upb0_2_m_j, upb0_2_p_k, upb0_2_m_k);
                    const double ubxy = dd(ay, by, cy, upb0_1_p_i, upb0_1_m_i, upb0_1_p_j, upb0_1_m_j, upb0_1_p_k, upb0_1_m_k);

                    double srcBx = (ubzy - ubyz);
                    double srcBy = (ubxz - ubzx);
                    double srcBz = (ubyx - ubxy);

                    // -------- term 2: -a7 * curl( (J x Bt)/ne ) --------
                    if (a7 != 0.0)
                    {
                        double s_p[3], s_m[3];

                        // ujbxy etc use center derv + neighbor sjbne values (same as Fortran)
                        double sjbne_1_p_i[3], sjbne_1_m_i[3], sjbne_1_p_j[3], sjbne_1_m_j[3], sjbne_1_p_k[3], sjbne_1_m_k[3];
                        double sjbne_2_p_i[3], sjbne_2_m_i[3], sjbne_2_p_j[3], sjbne_2_m_j[3], sjbne_2_p_k[3], sjbne_2_m_k[3];
                        double sjbne_3_p_i[3], sjbne_3_m_i[3], sjbne_3_p_j[3], sjbne_3_m_j[3], sjbne_3_p_k[3], sjbne_3_m_k[3];

                        sjbne_at(i + 1, j, k, sjbne_1_p_i);
                        sjbne_at(i - 1, j, k, sjbne_1_m_i);
                        sjbne_at(i, j + 1, k, sjbne_1_p_j);
                        sjbne_at(i, j - 1, k, sjbne_1_m_j);
                        sjbne_at(i, j, k + 1, sjbne_1_p_k);
                        sjbne_at(i, j, k - 1, sjbne_1_m_k);

                        sjbne_at(i + 1, j, k, sjbne_2_p_i);
                        sjbne_at(i - 1, j, k, sjbne_2_m_i);
                        sjbne_at(i, j + 1, k, sjbne_2_p_j);
                        sjbne_at(i, j - 1, k, sjbne_2_m_j);
                        sjbne_at(i, j, k + 1, sjbne_2_p_k);
                        sjbne_at(i, j, k - 1, sjbne_2_m_k);

                        sjbne_at(i + 1, j, k, sjbne_3_p_i);
                        sjbne_at(i - 1, j, k, sjbne_3_m_i);
                        sjbne_at(i, j + 1, k, sjbne_3_p_j);
                        sjbne_at(i, j - 1, k, sjbne_3_m_j);
                        sjbne_at(i, j, k + 1, sjbne_3_p_k);
                        sjbne_at(i, j, k - 1, sjbne_3_m_k);

                        // Derivatives of sjbne components using center derv (Fortran ujbxy...):
                        const double ujbxy = dd(ay, by, cy, sjbne_1_p_i[0], sjbne_1_m_i[0], sjbne_1_p_j[0], sjbne_1_m_j[0], sjbne_1_p_k[0], sjbne_1_m_k[0]);
                        const double ujbxz = dd(az, bz, cz, sjbne_1_p_i[0], sjbne_1_m_i[0], sjbne_1_p_j[0], sjbne_1_m_j[0], sjbne_1_p_k[0], sjbne_1_m_k[0]);

                        const double ujbyx = dd(ax, bx, cx, sjbne_2_p_i[1], sjbne_2_m_i[1], sjbne_2_p_j[1], sjbne_2_m_j[1], sjbne_2_p_k[1], sjbne_2_m_k[1]);
                        const double ujbyz = dd(az, bz, cz, sjbne_2_p_i[1], sjbne_2_m_i[1], sjbne_2_p_j[1], sjbne_2_m_j[1], sjbne_2_p_k[1], sjbne_2_m_k[1]);

                        const double ujbzx = dd(ax, bx, cx, sjbne_3_p_i[2], sjbne_3_m_i[2], sjbne_3_p_j[2], sjbne_3_m_j[2], sjbne_3_p_k[2], sjbne_3_m_k[2]);
                        const double ujbzy = dd(ay, by, cy, sjbne_3_p_i[2], sjbne_3_m_i[2], sjbne_3_p_j[2], sjbne_3_m_j[2], sjbne_3_p_k[2], sjbne_3_m_k[2]);

                        const double curl_sjbne_x = (ujbzy - ujbyz);
                        const double curl_sjbne_y = (ujbxz - ujbzx);
                        const double curl_sjbne_z = (ujbyx - ujbxy);

                        srcBx += -a7 * curl_sjbne_x;
                        srcBy += -a7 * curl_sjbne_y;
                        srcBz += -a7 * curl_sjbne_z;
                    }

                    // -------- term 3: +a8 * (grad(1/ne) x grad(pe)) --------
                    if (a8 != 0.0)
                    {
                        const double invne_p_i = 1.0 / std::max(ne_at(i + 1, j, k), ne_floor);
                        const double invne_m_i = 1.0 / std::max(ne_at(i - 1, j, k), ne_floor);
                        const double invne_p_j = 1.0 / std::max(ne_at(i, j + 1, k), ne_floor);
                        const double invne_m_j = 1.0 / std::max(ne_at(i, j - 1, k), ne_floor);
                        const double invne_p_k = 1.0 / std::max(ne_at(i, j, k + 1), ne_floor);
                        const double invne_m_k = 1.0 / std::max(ne_at(i, j, k - 1), ne_floor);

                        const double pe_p_i = pe_at(i + 1, j, k), pe_m_i = pe_at(i - 1, j, k);
                        const double pe_p_j = pe_at(i, j + 1, k), pe_m_j = pe_at(i, j - 1, k);
                        const double pe_p_k = pe_at(i, j, k + 1), pe_m_k = pe_at(i, j, k - 1);

                        const double dpex = dd(ax, bx, cx, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        const double dpey = dd(ay, by, cy, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        const double dpez = dd(az, bz, cz, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);

                        const double dnex = dd(ax, bx, cx, invne_p_i, invne_m_i, invne_p_j, invne_m_j, invne_p_k, invne_m_k);
                        const double dney = dd(ay, by, cy, invne_p_i, invne_m_i, invne_p_j, invne_m_j, invne_p_k, invne_m_k);
                        const double dnez = dd(az, bz, cz, invne_p_i, invne_m_i, invne_p_j, invne_m_j, invne_p_k, invne_m_k);

                        srcBx += a8 * (dney * dpez - dnez * dpey);
                        srcBy += a8 * (dnez * dpex - dnex * dpez);
                        srcBz += a8 * (dnex * dpey - dney * dpex);
                    }

                    RHS(i, j, k, 0) += srcBx;
                    RHS(i, j, k, 1) += srcBy;
                    RHS(i, j, k, 2) += srcBz;
                }
    }

    // === Resistive magnetic diffusion (Mercury inner/shell):  -curl( yita0 * curl(Ub) ) / Rem8 ===
    // Put this near the end of MercurySolver::AddSourceToRHS_B()

    auto cst = par_->GetDou_List("constant").data;
    auto ref = par_->GetDou_List("REF").data;

    // mu_mag should be μ0 (physical). yitamax: max resistivity in SI (Ohm·m) in the legacy code convention.
    const double mu0 = cst.at("mu_mag");
    const double yitamax = par_->GetDou("eta_max_mercury"); // default matches legacy 1.25e7
    const double L_ref = ref.at("L_ref");                   // meters
    const double U_ref = ref.at("U");                       // m/s

    // Rem8 = L_ref * U_ref * mu0 / yitamax  => 1/Rem8 = yitamax / (mu0*L_ref*U_ref)
    const double invRem8 = yitamax / (mu0 * L_ref * U_ref);

    // shell profile parameters (defaults match the legacy Fortran intent)
    const double r_cut_in = (cst.count("yita_r_cut_in") ? cst.at("yita_r_cut_in") : 0.8);
    const double r_cut_out = (cst.count("yita_r_cut_out") ? cst.at("yita_r_cut_out") : 1.0);
    const double r0 = (cst.count("yita_r0") ? cst.at("yita_r0") : 0.8);
    const double r1 = (cst.count("yita_r1") ? cst.at("yita_r1") : 1.00);
    const double w = (cst.count("yita_w") ? cst.at("yita_w") : 0.01);

    auto yita0_of_r = [&](double r) -> double
    {
        if (r <= r_cut_in || r >= r_cut_out)
            return 0.0;
        // 0.5*(tanh((r-r0)/w) - tanh((r-r1)/w)) : ~1 inside [r0,r1], smooth edges
        return 0.5 * (std::tanh((r - r0) / w) - std::tanh((r - r1) / w));
    };

    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        FieldBlock &Ub = fld_->field(fid_.fid_U_b, ib); // induced B (qb)
        FieldBlock &RHS = fld_->field(fid_.fid_RHS_b, ib);

        if (!Jac.is_allocated() || !Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated())
            continue;
        if (!Ub.is_allocated() || !RHS.is_allocated())
            continue;

        Block &blk = fld_->grd->grids(ib);

        if (blk.block_name != "Solid")
            continue;

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        // reuse your derv_at / dd / J_at by capturing required blocks
        auto derv_at = [&](int i, int j, int k,
                           double &ax, double &ay, double &az,
                           double &bx, double &by, double &bz,
                           double &cx, double &cy, double &cz)
        {
            const double V = std::abs(Jac(i, j, k, 0));
            if (V <= 0.0)
            {
                ax = ay = az = bx = by = bz = cx = cy = cz = 0.0;
                return;
            }

            ax = 0.5 * (Axi(i - 1, j, k, 0) + Axi(i, j, k, 0)) / V;
            ay = 0.5 * (Axi(i - 1, j, k, 1) + Axi(i, j, k, 1)) / V;
            az = 0.5 * (Axi(i - 1, j, k, 2) + Axi(i, j, k, 2)) / V;

            bx = 0.5 * (Aet(i, j - 1, k, 0) + Aet(i, j, k, 0)) / V;
            by = 0.5 * (Aet(i, j - 1, k, 1) + Aet(i, j, k, 1)) / V;
            bz = 0.5 * (Aet(i, j - 1, k, 2) + Aet(i, j, k, 2)) / V;

            cx = 0.5 * (Aze(i, j, k - 1, 0) + Aze(i, j, k, 0)) / V;
            cy = 0.5 * (Aze(i, j, k - 1, 1) + Aze(i, j, k, 1)) / V;
            cz = 0.5 * (Aze(i, j, k - 1, 2) + Aze(i, j, k, 2)) / V;
        };

        auto dd = [](double a, double b, double c,
                     double fp_i, double fm_i,
                     double fp_j, double fm_j,
                     double fp_k, double fm_k) -> double
        {
            return 0.5 * (a * (fp_i - fm_i) + b * (fp_j - fm_j) + c * (fp_k - fm_k));
        };

        auto J_at = [&](int i, int j, int k, double J[3])
        {
            double ax, ay, az, bx, by, bz, cx, cy, cz;
            derv_at(i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

            const double Bx_p_i = Ub(i + 1, j, k, 0), Bx_m_i = Ub(i - 1, j, k, 0);
            const double Bx_p_j = Ub(i, j + 1, k, 0), Bx_m_j = Ub(i, j - 1, k, 0);
            const double Bx_p_k = Ub(i, j, k + 1, 0), Bx_m_k = Ub(i, j, k - 1, 0);

            const double By_p_i = Ub(i + 1, j, k, 1), By_m_i = Ub(i - 1, j, k, 1);
            const double By_p_j = Ub(i, j + 1, k, 1), By_m_j = Ub(i, j - 1, k, 1);
            const double By_p_k = Ub(i, j, k + 1, 1), By_m_k = Ub(i, j, k - 1, 1);

            const double Bz_p_i = Ub(i + 1, j, k, 2), Bz_m_i = Ub(i - 1, j, k, 2);
            const double Bz_p_j = Ub(i, j + 1, k, 2), Bz_m_j = Ub(i, j - 1, k, 2);
            const double Bz_p_k = Ub(i, j, k + 1, 2), Bz_m_k = Ub(i, j, k - 1, 2);

            const double Bzy = dd(ay, by, cy, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);
            const double Byz = dd(az, bz, cz, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);

            const double Bxz = dd(az, bz, cz, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);
            const double Bzx = dd(ax, bx, cx, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);

            const double Byx = dd(ax, bx, cx, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);
            const double Bxy = dd(ay, by, cy, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);

            J[0] = (Bzy - Byz);
            J[1] = (Bxz - Bzx);
            J[2] = (Byx - Bxy);
        };

        auto F_at = [&](int i, int j, int k, double F[3])
        {
            const double x = blk.dual_x(i + 1, j + 1, k + 1);
            const double y = blk.dual_y(i + 1, j + 1, k + 1);
            const double z = blk.dual_z(i + 1, j + 1, k + 1);
            const double r = std::sqrt(x * x + y * y + z * z);

            const double yita0 = yita0_of_r(r);
            if (yita0 == 0.0)
            {
                F[0] = F[1] = F[2] = 0.0;
                return;
            }

            double J[3];
            J_at(i, j, k, J);

            // F = yita0 * curl(B)
            F[0] = yita0 * J[0];
            F[1] = yita0 * J[1];
            F[2] = yita0 * J[2];
        };

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // quick reject using center radius (optional)
                    const double x = blk.dual_x(i + 1, j + 1, k + 1);
                    const double y = blk.dual_y(i + 1, j + 1, k + 1);
                    const double z = blk.dual_z(i + 1, j + 1, k + 1);
                    const double r = std::sqrt(x * x + y * y + z * z);
                    if (yita0_of_r(r) == 0.0)
                        continue;

                    double ax, ay, az, bx, by, bz, cx, cy, cz;
                    derv_at(i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

                    double Fpi[3], Fmi[3], Fpj[3], Fmj[3], Fpk[3], Fmk[3];
                    F_at(i + 1, j, k, Fpi);
                    F_at(i - 1, j, k, Fmi);
                    F_at(i, j + 1, k, Fpj);
                    F_at(i, j - 1, k, Fmj);
                    F_at(i, j, k + 1, Fpk);
                    F_at(i, j, k - 1, Fmk);

                    // curl(F): (dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy)
                    const double Fzy = dd(ay, by, cy, Fpi[2], Fmi[2], Fpj[2], Fmj[2], Fpk[2], Fmk[2]);
                    const double Fyz = dd(az, bz, cz, Fpi[1], Fmi[1], Fpj[1], Fmj[1], Fpk[1], Fmk[1]);

                    const double Fxz = dd(az, bz, cz, Fpi[0], Fmi[0], Fpj[0], Fmj[0], Fpk[0], Fmk[0]);
                    const double Fzx = dd(ax, bx, cx, Fpi[2], Fmi[2], Fpj[2], Fmj[2], Fpk[2], Fmk[2]);

                    const double Fyx = dd(ax, bx, cx, Fpi[1], Fmi[1], Fpj[1], Fmj[1], Fpk[1], Fmk[1]);
                    const double Fxy = dd(ay, by, cy, Fpi[0], Fmi[0], Fpj[0], Fmj[0], Fpk[0], Fmk[0]);

                    const double curlFx = (Fzy - Fyz);
                    const double curlFy = (Fxz - Fzx);
                    const double curlFz = (Fyx - Fxy);

                    RHS(i, j, k, 0) += -invRem8 * curlFx;
                    RHS(i, j, k, 1) += -invRem8 * curlFy;
                    RHS(i, j, k, 2) += -invRem8 * curlFz;
                }
    }
}