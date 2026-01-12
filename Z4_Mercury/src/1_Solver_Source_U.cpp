#include "MercurySolver.h"

// 添加 Fortran source_species 中对流体 dq 的“生效项”
// 依赖字段：U_H/U_Na, PV_H/PV_Na, U_plus, B_cell, U_b(用于curl), Na(neutral), Photo_rate, Jac, metric(Axi/Aet/Aze)
void MercurySolver::AddSourceToRHS_Fluid()
{
    // ---------- constants / refs ----------
    const auto &C = par_->GetDou_List("constant").data;
    const auto &R = par_->GetDou_List("REF").data;

    const double NA = C.at("NA");
    const double q_e = C.at("q_e");    // e (Coulomb)
    const double mu0 = C.at("mu_mag"); // μ0

    const double L_ref = R.at("L_ref");         // m
    const double U_ref = R.at("U");             // m/s
    const double B_ref = R.at("B_ref");         // Tesla
    const double n_ref = R.at("n");             // 1/m^3
    const double Mref = R.at("Molecular_mass"); // kg/mol (molar mass)

    const double rho_ref = (Mref / NA) * n_ref; // kg/m^3
    const double m_H = (Mref / NA);             // kg per particle (assume ref mass = proton)
    const double m_Na = 23.0 * m_H;

    // coefficients matching Fortran structure (see explanation in my previous message)
    const double a2 = (1e6 * q_e * L_ref * B_ref) / (rho_ref * U_ref);
    const double a3 = (B_ref * B_ref) / (mu0 * rho_ref * U_ref * U_ref);
    const double a1_Na = (L_ref / U_ref) * (m_Na / rho_ref) * 1e6; // mass source coeff for Na+

    // a4 = L_ref / U_ref(秒)，Fortran : a4 = sl8 / v8
    const double a4 = (L_ref / U_ref);

    const double ne_floor = 1e-30;

    // dd: Fortran 0.5*(a*(f(i+1)-f(i-1)) + b*(f(j+1)-f(j-1)) + c*(f(k+1)-f(k-1)))
    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) + b * (fp_j - fm_j) + c * (fp_k - fm_k));
    };

    // derv_at: 从 (Axi/Aet/Aze)/Jac 还原 derv(1..3,1..3) 的 ax..cz（和 1_Solver_Scheme_B.cpp 一致）
    auto derv_at = [&](FieldBlock &Jac, FieldBlock &Axi, FieldBlock &Aet, FieldBlock &Aze,
                       int i, int j, int k,
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

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UNa = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib);
        FieldBlock &PVN = fld_->field(fid_.fid_PV_Na, ib);

        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib);
        FieldBlock &Ub = fld_->field(fid_.fid_U_b, ib);
        FieldBlock &Bt = fld_->field(fid_.fid_Bcell, ib);

        FieldBlock &NaNeu = fld_->field(fld_->field_id("Na"), ib);
        FieldBlock &Photo = fld_->field(fld_->field_id("Photo_rate"), ib);

        FieldBlock &RHS_H = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHS_Na = fld_->field(fid_.fid_RHS_Na, ib);

        if (!Jac.is_allocated() || !Axi.is_allocated())
            continue;
        if (!UH.is_allocated() || !UNa.is_allocated())
            continue;
        if (!PVH.is_allocated() || !PVN.is_allocated())
            continue;
        if (!Up.is_allocated() || !Ub.is_allocated() || !Bt.is_allocated())
            continue;
        if (!RHS_H.is_allocated() || !RHS_Na.is_allocated())
            continue;
        if (!NaNeu.is_allocated() || !Photo.is_allocated())
            continue;

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // ---------- basic primitives ----------
                    const double uH = PVH(i, j, k, 0), vH = PVH(i, j, k, 1), wH = PVH(i, j, k, 2);
                    const double uN = PVN(i, j, k, 0), vN = PVN(i, j, k, 1), wN = PVN(i, j, k, 2);

                    const double upx = Up(i, j, k, 0), upy = Up(i, j, k, 1), upz = Up(i, j, k, 2);

                    const double Bx = Bt(i, j, k, 0), By = Bt(i, j, k, 1), Bz = Bt(i, j, k, 2);

                    // ---------- number densities in cm^-3 (to match Fortran roenum) ----------
                    const double rhoH_nd = std::max(UH(i, j, k, 0), 0.0);
                    const double rhoNa_nd = std::max(UNa(i, j, k, 0), 0.0);

                    const double nH_cm = (rho_ref * rhoH_nd) / m_H * 1e-6;
                    const double nNa_cm = (rho_ref * rhoNa_nd) / m_Na * 1e-6;
                    const double ne_cm = std::max(nH_cm + nNa_cm, ne_floor);

                    // -------------------- metric derv(ax..cz) --------------------
                    double ax, ay, az, bx, by, bz, cx, cy, cz;
                    derv_at(Jac, Axi, Aet, Aze, i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

                    // ---------- J and JxB ----------
                    double J[3];
                    {
                        const double Bx_p_i = Ub(i + 1, j, k, 0), Bx_m_i = Ub(i - 1, j, k, 0);
                        const double By_p_i = Ub(i + 1, j, k, 1), By_m_i = Ub(i - 1, j, k, 1);
                        const double Bz_p_i = Ub(i + 1, j, k, 2), Bz_m_i = Ub(i - 1, j, k, 2);

                        double Bx_p_j = 0, Bx_m_j = 0, By_p_j = 0, By_m_j = 0, Bz_p_j = 0, Bz_m_j = 0;
                        double Bx_p_k = 0, Bx_m_k = 0, By_p_k = 0, By_m_k = 0, Bz_p_k = 0, Bz_m_k = 0;

                        Bx_p_j = Ub(i, j + 1, k, 0);
                        Bx_m_j = Ub(i, j - 1, k, 0);
                        By_p_j = Ub(i, j + 1, k, 1);
                        By_m_j = Ub(i, j - 1, k, 1);
                        Bz_p_j = Ub(i, j + 1, k, 2);
                        Bz_m_j = Ub(i, j - 1, k, 2);

                        Bx_p_k = Ub(i, j, k + 1, 0);
                        Bx_m_k = Ub(i, j, k - 1, 0);
                        By_p_k = Ub(i, j, k + 1, 1);
                        By_m_k = Ub(i, j, k - 1, 1);
                        Bz_p_k = Ub(i, j, k + 1, 2);
                        Bz_m_k = Ub(i, j, k - 1, 2);

                        const double Bzy = dd(ay, by, cy, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);
                        const double Byz = dd(az, bz, cz, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);

                        const double Bxz = dd(az, bz, cz, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);
                        const double Bzx = dd(ax, bx, cx, Bz_p_i, Bz_m_i, Bz_p_j, Bz_m_j, Bz_p_k, Bz_m_k);

                        const double Byx = dd(ax, bx, cx, By_p_i, By_m_i, By_p_j, By_m_j, By_p_k, By_m_k);
                        const double Bxy = dd(ay, by, cy, Bx_p_i, Bx_m_i, Bx_p_j, Bx_m_j, Bx_p_k, Bx_m_k);

                        J[0] = Bzy - Byz;
                        J[1] = Bxz - Bzx;
                        J[2] = Byx - Bxy;
                    }

                    // JxB (对应 Fortran sjb)
                    const double sjbx = J[1] * Bz - J[2] * By;
                    const double sjby = J[2] * Bx - J[0] * Bz;
                    const double sjbz = J[0] * By - J[1] * Bx;

                    // -------------------- grad(pe) --------------------
                    double dpex, dpey, dpez;
                    {
                        const double pe_p_i = PVH(i + 1, j, k, 3) + PVN(i + 1, j, k, 3);
                        const double pe_m_i = PVH(i - 1, j, k, 3) + PVN(i - 1, j, k, 3);

                        double pe_p_j = 0, pe_m_j = 0, pe_p_k = 0, pe_m_k = 0;

                        pe_p_j = PVH(i, j + 1, k, 3) + PVN(i, j + 1, k, 3);
                        pe_m_j = PVH(i, j - 1, k, 3) + PVN(i, j - 1, k, 3);

                        pe_p_k = PVH(i, j, k + 1, 3) + PVN(i, j, k + 1, 3);
                        pe_m_k = PVH(i, j, k - 1, 3) + PVN(i, j, k - 1, 3);

                        dpex = dd(ax, bx, cx, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        dpey = dd(ay, by, cy, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        dpez = dd(az, bz, cz, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                    }

                    // ---------- ionization source sss (cm^-3 s^-1), only for Na+ ----------
                    const double sss = Photo(i, j, k, 0); // cm^-3 s^-1

                    // -------------------- drag frequency vst (1/s) --------------------
                    // Fortran: if(x>=0) vst=sk1 else vst=sk2  （x 为无量纲坐标）
                    // C++ 目前没有直接取物理坐标，这里用 Jac/metric 无法得到 x；
                    // 因此：如果你要严格复现 Fortran 的 x 判定，请从几何场里取坐标场（例如 fid_.fid_coord）再判断。
                    // 这里给一个保守默认：统一用 sk1（你也可以改成 sk2 或按你已有的坐标场实现）
                    // const double vst = sk1;

                    // =====================
                    // species H+  (ls=1)
                    // =====================
                    {
                        const double sns0 = nH_cm;

                        // Fortran-simplified sub:
                        const double subx = (vH - upy) * Bz - (wH - upz) * By;
                        const double suby = (wH - upz) * Bx - (uH - upx) * Bz;
                        const double subz = (uH - upx) * By - (vH - upy) * Bx;

                        const double sjbu = sjbx * uH + sjby * vH + sjbz * wH;
                        const double dpeu = dpex * uH + dpey * vH + dpez * wH;

                        const double subu = subx * uH + suby * vH + subz * wH;

                        RHS_H(i, j, k, 0) += 0.0; // H+ has no mass creation in Fortran here
                        RHS_H(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm);
                        RHS_H(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm);
                        RHS_H(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm);
                        RHS_H(i, j, k, 4) += a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm);

                        // RHS_H(i, j, k, 1) += -sns0 * (dpex / ne_cm);
                        // RHS_H(i, j, k, 2) += -sns0 * (dpey / ne_cm);
                        // RHS_H(i, j, k, 3) += -sns0 * (dpez / ne_cm);

                        // const double roes0 = rhoH_nd;
                        // RHS_H(i, j, k, 1) += -a4 * roes0 * uH * vst;
                        // RHS_H(i, j, k, 2) += -a4 * roes0 * vH * vst;
                        // RHS_H(i, j, k, 3) += -a4 * roes0 * wH * vst;

                        // RHS_H(i, j, k, 4) += a2 * sns0 * subu;
                        // // 若你把动量里加了 -rho*nu*u，这里建议同步加入做功项（否则能量方程会“凭空增温/减温”）
                        // const double roes0 = rhoH_nd;
                        // const double us2 = uH * uH + vH * vH + wH * wH;
                        // RHS_H(i, j, k, 4) += -a4 * roes0 * us2 * vst;
                    }

                    // =====================
                    // species Na+ (ls=2)
                    // =====================
                    {
                        const double sns0 = nNa_cm;

                        const double subx = (vN - upy) * Bz - (wN - upz) * By;
                        const double suby = (wN - upz) * Bx - (uN - upx) * Bz;
                        const double subz = (uN - upx) * By - (vN - upy) * Bx;

                        const double sjbu = sjbx * uN + sjby * vN + sjbz * wN;
                        const double dpeu = dpex * uN + dpey * vN + dpez * wN;

                        const double subu = subx * uN + suby * vN + subz * wN;

                        RHS_Na(i, j, k, 0) += a1_Na * sss; // Na+ mass creation
                        RHS_Na(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm);
                        RHS_Na(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm);
                        RHS_Na(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm);
                        RHS_Na(i, j, k, 4) += a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm);

                        // RHS_Na(i, j, k, 1) += -sns0 * (dpex / ne_cm);
                        // RHS_Na(i, j, k, 2) += -sns0 * (dpey / ne_cm);
                        // RHS_Na(i, j, k, 3) += -sns0 * (dpez / ne_cm);

                        // const double roes0 = rhoNa_nd;
                        // RHS_Na(i, j, k, 1) += -a4 * roes0 * uN * vst;
                        // RHS_Na(i, j, k, 2) += -a4 * roes0 * vN * vst;
                        // RHS_Na(i, j, k, 3) += -a4 * roes0 * wN * vst;

                        // RHS_Na(i, j, k, 4) += a2 * sns0 * subu;

                        // const double roes0 = rhoNa_nd;
                        // const double us2 = uN * uN + vN * vN + wN * wN;
                        // RHS_Na(i, j, k, 4) += -a4 * roes0 * us2 * vst;
                    }
                }
    }
}