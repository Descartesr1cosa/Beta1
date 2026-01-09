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

    const double ne_floor = 1e-30;

    auto derv_a = [&](FieldBlock &Jac, FieldBlock &Axi, FieldBlock &Aet, FieldBlock &Aze,
                      int i, int j, int k, double &ax, double &ay, double &az)
    {
        // Axi(i,j,k) is i+1/2 face -> cell i uses faces (i-1) and (i)
        const double V = std::abs(Jac(i, j, k, 0));
        if (V <= 0.0)
        {
            ax = ay = az = 0.0;
            return;
        }
        ax = 0.5 * (Axi(i - 1, j, k, 0) + Axi(i, j, k, 0)) / V;
        ay = 0.5 * (Axi(i - 1, j, k, 1) + Axi(i, j, k, 1)) / V;
        az = 0.5 * (Axi(i - 1, j, k, 2) + Axi(i, j, k, 2)) / V;
    };

    // curl(U_b) 仅保留 i 方向差分（与 Fortran 保持一致：by/cy 等项在原代码里被注释）
    auto curl_Ub_iOnly = [&](FieldBlock &Jac, FieldBlock &Axi, FieldBlock &Aet, FieldBlock &Aze,
                             FieldBlock &Ub, int i, int j, int k, double J[3])
    {
        double ax, ay, az;
        derv_a(Jac, Axi, Aet, Aze, i, j, k, ax, ay, az);

        const double Bx_p = Ub(i + 1, j, k, 0), Bx_m = Ub(i - 1, j, k, 0);
        const double By_p = Ub(i + 1, j, k, 1), By_m = Ub(i - 1, j, k, 1);
        const double Bz_p = Ub(i + 1, j, k, 2), Bz_m = Ub(i - 1, j, k, 2);

        // Fortran:
        // Bzy = ay*(Bz(i+1)-Bz(i-1))*0.5
        // Byz = az*(By(i+1)-By(i-1))*0.5
        // Bxz = az*(Bx(i+1)-Bx(i-1))*0.5
        // Bzx = ax*(Bz(i+1)-Bz(i-1))*0.5
        // Byx = ax*(By(i+1)-By(i-1))*0.5
        // Bxy = ay*(Bx(i+1)-Bx(i-1))*0.5
        const double Bzy = 0.5 * ay * (Bz_p - Bz_m);
        const double Byz = 0.5 * az * (By_p - By_m);

        const double Bxz = 0.5 * az * (Bx_p - Bx_m);
        const double Bzx = 0.5 * ax * (Bz_p - Bz_m);

        const double Byx = 0.5 * ax * (By_p - By_m);
        const double Bxy = 0.5 * ay * (Bx_p - Bx_m);

        J[0] = Bzy - Byz;
        J[1] = Bxz - Bzx;
        J[2] = Byx - Bxy;
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

                    // ---------- J and JxB ----------
                    double J[3];
                    curl_Ub_iOnly(Jac, Axi, Aet, Aze, Ub, i, j, k, J);

                    const double sjbx = J[1] * Bz - J[2] * By;
                    const double sjby = J[2] * Bx - J[0] * Bz;
                    const double sjbz = J[0] * By - J[1] * Bx;

                    // ---------- grad(pe) i-only like Fortran dpe(:) ----------
                    const double pe = PVH(i, j, k, 3) + PVN(i, j, k, 3);

                    double ax, ay, az;
                    derv_a(Jac, Axi, Aet, Aze, i, j, k, ax, ay, az);

                    const double pe_p = (PVH(i + 1, j, k, 3) + PVN(i + 1, j, k, 3));
                    const double pe_m = (PVH(i - 1, j, k, 3) + PVN(i - 1, j, k, 3));

                    const double dpex = 0.5 * ax * (pe_p - pe_m);
                    const double dpey = 0.5 * ay * (pe_p - pe_m);
                    const double dpez = 0.5 * az * (pe_p - pe_m);

                    // ---------- ionization source sss (cm^-3 s^-1), only for Na+ ----------
                    const double sss = Photo(i, j, k, 0) * NaNeu(i, j, k, 0);

                    // =====================
                    // species H+  (ls=1)
                    // =====================
                    {
                        const double sns0 = nH_cm;

                        // Fortran-simplified sub:
                        const double subx = (vH - upy) * Bz;
                        const double suby = (wH - upz) * Bx;
                        const double subz = (uH - upx) * By;

                        const double sjbu = sjbx * uH + sjby * vH + sjbz * wH;
                        const double dpeu = dpex * uH + dpey * vH + dpez * wH;

                        RHS_H(i, j, k, 0) += 0.0; // H+ has no mass creation in Fortran here
                        RHS_H(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm);
                        RHS_H(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm);
                        RHS_H(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm);
                        RHS_H(i, j, k, 4) += a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm);
                    }

                    // =====================
                    // species Na+ (ls=2)
                    // =====================
                    {
                        const double sns0 = nNa_cm;

                        const double subx = (vN - upy) * Bz;
                        const double suby = (wN - upz) * Bx;
                        const double subz = (uN - upx) * By;

                        const double sjbu = sjbx * uN + sjby * vN + sjbz * wN;
                        const double dpeu = dpex * uN + dpey * vN + dpez * wN;

                        RHS_Na(i, j, k, 0) += a1_Na * sss; // Na+ mass creation
                        RHS_Na(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm);
                        RHS_Na(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm);
                        RHS_Na(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm);
                        RHS_Na(i, j, k, 4) += a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm);
                    }
                }
    }
}