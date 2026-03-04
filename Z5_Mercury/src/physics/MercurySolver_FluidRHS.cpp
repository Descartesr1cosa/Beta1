#include "MercurySolver.h"

void MercurySolver::AssembleRHS_Fluid_()
{
    Scheme_U_();
    AddSourceToRHS_Fluid();
}

void MercurySolver::Scheme_U_()
{
    auto calc_Jac_radius_GCL = [&](double &out,
                                   const double rho, const double u, const double v, const double w, const double p,
                                   const double k1, const double k2, const double k3)
    {
        // k1,k2,k3 are "Jac * contravariant metric components" (GCL form)
        const double uvw = k1 * u + k2 * v + k3 * w;
        const double kk = (k1 * k1 + k2 * k2 + k3 * k3);
        const double c = std::sqrt(std::max(0.0, gamma_ * p / std::max(rho, 1e-30)));
        const double a = std::abs(uvw) + c * std::sqrt(kk);
        out = a;
    };

    auto calc_Jac_Flux_GCL = [&](double *flux,
                                 const double rho, const double u, const double v, const double w, const double p,
                                 const double k1, const double k2, const double k3)
    {
        const double uvw = k1 * u + k2 * v + k3 * w;
        const double rhoe = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

        flux[0] = rho * uvw;
        flux[1] = rho * uvw * u + k1 * p;
        flux[2] = rho * uvw * v + k2 * p;
        flux[3] = rho * uvw * w + k3 * p;
        flux[4] = uvw * (rhoe + p);
    };

    // One-face Rusanov (piecewise constant; you can upgrade to MUSCL later)
    auto Reconstruction_Rusanov = [&](const double *metric, int direction,
                                      FieldBlock &PV, FieldBlock &U,
                                      int i, int j, int k,
                                      double *out_flux)
    {
        const double k1 = metric[0];
        const double k2 = metric[1];
        const double k3 = metric[2];

        // left/right states at the face: (i-1,i) or (j-1,j) or (k-1,k)
        int iL = i, jL = j, kL = k;
        int iR = i, jR = j, kR = k;

        if (direction == 0)
        {
            iL = i - 1;
            iR = i;
        }
        if (direction == 1)
        {
            jL = j - 1;
            jR = j;
        }
        if (direction == 2)
        {
            kL = k - 1;
            kR = k;
        }

        // Left primitive
        const double rhoL = U(iL, jL, kL, 0);
        const double uL = PV(iL, jL, kL, 0);
        const double vL = PV(iL, jL, kL, 1);
        const double wL = PV(iL, jL, kL, 2);
        const double pL = PV(iL, jL, kL, 3);

        // Right primitive
        const double rhoR = U(iR, jR, kR, 0);
        const double uR = PV(iR, jR, kR, 0);
        const double vR = PV(iR, jR, kR, 1);
        const double wR = PV(iR, jR, kR, 2);
        const double pR = PV(iR, jR, kR, 3);

        // Conservative UL/UR
        double UL[5], UR[5];
        UL[0] = rhoL;
        UL[1] = rhoL * uL;
        UL[2] = rhoL * vL;
        UL[3] = rhoL * wL;
        UL[4] = pL / (gamma_ - 1.0) + 0.5 * rhoL * (uL * uL + vL * vL + wL * wL);

        UR[0] = rhoR;
        UR[1] = rhoR * uR;
        UR[2] = rhoR * vR;
        UR[3] = rhoR * wR;
        UR[4] = pR / (gamma_ - 1.0) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);

        // Spectral radius
        double radL = 0.0, radR = 0.0;
        calc_Jac_radius_GCL(radL, rhoL, uL, vL, wL, pL, k1, k2, k3);
        calc_Jac_radius_GCL(radR, rhoR, uR, vR, wR, pR, k1, k2, k3);
        const double rad = std::max(radL, radR);

        // Fluxes
        double FL[5], FR[5];
        calc_Jac_Flux_GCL(FL, rhoL, uL, vL, wL, pL, k1, k2, k3);
        calc_Jac_Flux_GCL(FR, rhoR, uR, vR, wR, pR, k1, k2, k3);

        for (int m = 0; m < 5; ++m)
            out_flux[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * rad * (UR[m] - UL[m]);
    };

    auto do_one_species = [&](int fidU, int fidPV, int fidRHS)
    {
        for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
        {
            FieldBlock &U = fld_->field(fidU, iblk);
            FieldBlock &PV = fld_->field(fidPV, iblk);
            FieldBlock &RHS = fld_->field(fidRHS, iblk);

            FieldBlock &Jac = fld_->field(fid_.fid_Jac, iblk);
            FieldBlock &XI = fld_->field(fid_.fid_metric.xi, iblk);   // FaceXi,3
            FieldBlock &ET = fld_->field(fid_.fid_metric.eta, iblk);  // FaceEt,3
            FieldBlock &ZE = fld_->field(fid_.fid_metric.zeta, iblk); // FaceZe,3

            if (!U.is_allocated() || !PV.is_allocated() || !RHS.is_allocated())
                continue;

            // temp flux fields (shared is OK; we overwrite per species)
            FieldBlock &Fxi = fld_->field("F_xi", iblk);
            FieldBlock &Fet = fld_->field("F_eta", iblk);
            FieldBlock &Fze = fld_->field("F_zeta", iblk);

            double metric[3];
            double Flux[5];

            // ---- xi faces ----
            {
                Int3 sub = Fxi.inner_lo();
                Int3 sup = Fxi.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = XI(i, j, k, 0);
                            metric[1] = XI(i, j, k, 1);
                            metric[2] = XI(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 0, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fxi(i, j, k, m) = Flux[m];
                        }
            }

            // ---- eta faces ----
            {
                Int3 sub = Fet.inner_lo();
                Int3 sup = Fet.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = ET(i, j, k, 0);
                            metric[1] = ET(i, j, k, 1);
                            metric[2] = ET(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 1, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fet(i, j, k, m) = Flux[m];
                        }
            }

            // ---- zeta faces ----
            {
                Int3 sub = Fze.inner_lo();
                Int3 sup = Fze.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = ZE(i, j, k, 0);
                            metric[1] = ZE(i, j, k, 1);
                            metric[2] = ZE(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 2, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fze(i, j, k, m) = Flux[m];
                        }
            }

            // ---- divergence to RHS ----
            {
                Int3 sub = RHS.inner_lo();
                Int3 sup = RHS.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            const double invJ = 1.0 / Jac(i, j, k, 0);
                            for (int m = 0; m < 5; ++m)
                            {
                                RHS(i, j, k, m) -= (Fxi(i + 1, j, k, m) - Fxi(i, j, k, m)) * invJ;
                                RHS(i, j, k, m) -= (Fet(i, j + 1, k, m) - Fet(i, j, k, m)) * invJ;
                                RHS(i, j, k, m) -= (Fze(i, j, k + 1, m) - Fze(i, j, k, m)) * invJ;
                            }
                        }
            }
        }
    };

    // H+ and Na+ (no split, still one function)
    do_one_species(fid_.fid_U_H, fid_.fid_PV_H, fid_.fid_RHS_H);
    do_one_species(fid_.fid_U_Na, fid_.fid_PV_Na, fid_.fid_RHS_Na);
}

// 添加 Fortran source_species 中对流体 dq 的“生效项”
// 依赖字段：U_H/U_Na, PV_H/PV_Na, U_plus, B_cell, U_b(用于curl), Na(neutral), Photo_rate, Jac, metric(Axi/Aet/Aze)
void MercurySolver::AddSourceToRHS_Fluid()
{

    double Rabs_max_l = 0.0;
    double Rrel_max_l = 0.0;
    int Rmax_ib_l = -1, Rmax_i_l = -1, Rmax_j_l = -1, Rmax_k_l = -1;
    double EM_eabs_max_l = 0.0;
    double EM_erel_max_l = 0.0;
    double EM_Find_abs_max_l = 0.0;

    // ---------- constants  ----------
    const double Tn0 = 185.0;  // K
    const double sk1 = 5.0e-5; // 1/s  (day side)
    const double sk2 = 1.0e-5; // 1/s  (night side)

    // coefficients matching Fortran structure (see explanation in my previous message)

    const double a1_Na = (L_ref / U_ref) * (m_Na / rho_ref) * 1e6; // mass source coeff for Na+

    // a4 = L_ref / U_ref(秒)，Fortran : a4 = sl8 / v8
    const double a4 = (L_ref / U_ref);

    const double a5 = (3.0 * L_ref * k_Boltz) / (rho_ref * U_ref * U_ref * U_ref);

    const double a6 = (1.0e6 * L_ref * k_Boltz) / (rho_ref * U_ref * U_ref * U_ref * (gamma_ - 1.0));

    const double ne_floor = 1e-2 * n_ref;

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
        FieldBlock &Bt = fld_->field(fid_.fid_Bcell, ib);

        FieldBlock &Jxi = fld_->field(fid_.fid_J.xi, ib);
        FieldBlock &Jet = fld_->field(fid_.fid_J.eta, ib);
        FieldBlock &Jze = fld_->field(fid_.fid_J.zeta, ib);

        FieldBlock &pinvGTc = fld_->field(fid_.fid_pinvGT_Cell, ib);

        FieldBlock &NaNeu = fld_->field(fld_->field_id("Na"), ib);
        FieldBlock &Photo = fld_->field(fld_->field_id("Photo_rate"), ib);

        FieldBlock &RHS_H = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHS_Na = fld_->field(fid_.fid_RHS_Na, ib);

        double3D &x = fld_->grd->grids(ib).x;

        if (!Jac.is_allocated() || !Axi.is_allocated())
            continue;
        if (!UH.is_allocated() || !UNa.is_allocated())
            continue;
        if (!PVH.is_allocated() || !PVN.is_allocated())
            continue;
        if (!Up.is_allocated() || !Bt.is_allocated())
            continue;
        if (!RHS_H.is_allocated() || !RHS_Na.is_allocated())
            continue;
        if (!NaNeu.is_allocated() || !Photo.is_allocated())
            continue;
        if (!Jxi.is_allocated() || !Jet.is_allocated() || !Jze.is_allocated())
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

                    // ---------- number densities in m^-3 ----------
                    const double rhoH_nd = std::max(UH(i, j, k, 0), 0.0);
                    const double rhoNa_nd = std::max(UNa(i, j, k, 0), 0.0);

                    const double nH_m = (rho_ref * rhoH_nd) / m_H;
                    const double nNa_m = (rho_ref * rhoNa_nd) / m_Na;
                    // const double ne_m = std::max(nH_m + nNa_m, ne_floor);

                    const double ne_true = nH_m + nNa_m;
                    const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor_dimensional * ne_hall_floor_dimensional);
                    const double s = ne_true / (ne_true + ne_hall_cut_dimensional); // 取同一个阈值，或用与 Hall 相同的量纲阈值

                    // -------------------- metric derv(ax..cz) --------------------
                    double ax, ay, az, bx, by, bz, cx, cy, cz;
                    derv_at(Jac, Axi, Aet, Aze, i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

                    // ---------- J and JxB ----------
                    double J[3];

                    {
                        // 1) Edge -> cell (4-pt symmetric average, same spirit as your Hall co-location)
                        const double Jxi_int = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j + 1, k, 0) +
                                                       Jxi(i, j, k + 1, 0) + Jxi(i, j + 1, k + 1, 0));

                        const double Jeta_int = 0.25 * (Jet(i, j, k, 0) + Jet(i + 1, j, k, 0) +
                                                        Jet(i, j, k + 1, 0) + Jet(i + 1, j, k + 1, 0));

                        const double Jzeta_int = 0.25 * (Jze(i, j, k, 0) + Jze(i + 1, j, k, 0) +
                                                         Jze(i, j + 1, k, 0) + Jze(i + 1, j + 1, k, 0));

                        // Reconstruct physical vector J(x,y,z) by: J = pinvGT_cell * w
                        // pinvGT_cell is row-major: [0..2; 3..5; 6..8]
                        J[0] = pinvGTc(i, j, k, 0) * Jxi_int + pinvGTc(i, j, k, 1) * Jeta_int + pinvGTc(i, j, k, 2) * Jzeta_int;
                        J[1] = pinvGTc(i, j, k, 3) * Jxi_int + pinvGTc(i, j, k, 4) * Jeta_int + pinvGTc(i, j, k, 5) * Jzeta_int;
                        J[2] = pinvGTc(i, j, k, 6) * Jxi_int + pinvGTc(i, j, k, 7) * Jeta_int + pinvGTc(i, j, k, 8) * Jzeta_int;
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

                        // dpex = dd(ax, bx, cx, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        // dpey = dd(ay, by, cy, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        // dpez = dd(az, bz, cz, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);

                        dpex = 0.0;
                        dpey = 0.0;
                        dpez = 0.0;
                    }

                    // ---------- ionization source sss (cm^-3 s^-1), only for Na+ ----------
                    const double sss = Photo(i, j, k, 0); // cm^-3 s^-1

                    // -------------------- drag frequency vst (1/s) --------------------
                    // Fortran: if(x>=0) vst=sk1 else vst=sk2  （x 为无量纲坐标）
                    // C++ 目前没有直接取物理坐标，这里用 Jac/metric 无法得到 x；
                    // 因此：如果你要严格复现 Fortran 的 x 判定，请从几何场里取坐标场（例如 fid_.fid_coord）再判断。
                    // 这里给一个保守默认：统一用 sk1（你也可以改成 sk2 或按你已有的坐标场实现）
                    // const double sk1 = 5E-5;
                    // const double sk2 = 1E-5;
                    const double vst = (x(i, j, k) >= 0) ? sk1 : sk2;

                    // b2 = (sm2/(sm1+sm2))*vst ;  b1 = (Tn0 - Ts0)*sm1/(sm1+sm2)
                    const double b2 = (m_Na / (m_H + m_Na)) * vst;

                    // sse = qm1 (Fortran), here Photo is already (cm^-3 s^-1) For electrics
                    const double sse = Photo(i, j, k, 0);

                    // =====================
                    // species H+  (ls=1)
                    // =====================
                    {
                        // Fortran-simplified sub:
                        const double subx = (vH - upy) * Bz - (wH - upz) * By;
                        const double suby = (wH - upz) * Bx - (uH - upx) * Bz;
                        const double subz = (uH - upx) * By - (vH - upy) * Bx;

                        const double sjbu = sjbx * uH + sjby * vH + sjbz * wH;
                        const double dpeu = dpex * uH + dpey * vH + dpez * wH;
                        const double subu = subx * uH + suby * vH + subz * wH;

                        // Ts0 in Kelvin
                        const double Ts0 = PVH(i, j, k, 4) * T_ref;
                        const double us2 = uH * uH + vH * vH + wH * wH;
                        const double b1 = (Tn0 - Ts0) * (m_H / (m_H + m_Na));

                        // RHS_H(i, j, k, 0) += 0.0; // H+ has no mass creation in Fortran here
                        // RHS_H(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm) - sns0 * (dpex / ne_cm) - a4 * rhoH_nd * uH * vst;
                        // RHS_H(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm) - sns0 * (dpey / ne_cm) - a4 * rhoH_nd * vH * vst;
                        // RHS_H(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm) - sns0 * (dpez / ne_cm) - a4 * rhoH_nd * wH * vst;
                        // RHS_H(i, j, k, 4) += a2 * sns0 * subu + a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm) + a4 * rhoH_nd * us2 * b2;

                        // RHS_H(i, j, k, 4) += a5 * sns0 * b1 + a6 * sns0 * sse * Tn0 / ne_cm; //+ a6 * 0.0 * Tn0 as sss = 0 For H+

                        RHS_H(i, j, k, 0) += 0.0;                                                                                    // H+ has no mass creation in Fortran here
                        RHS_H(i, j, k, 1) += momentum_induce_coeff * nH_m * subx + momentum_hall_coeff * s * nH_m * (sjbx / ne_eff); //- a4 * rhoH_nd * uH * vst;
                        RHS_H(i, j, k, 2) += momentum_induce_coeff * nH_m * suby + momentum_hall_coeff * s * nH_m * (sjby / ne_eff); //- a4 * rhoH_nd * vH * vst;
                        RHS_H(i, j, k, 3) += momentum_induce_coeff * nH_m * subz + momentum_hall_coeff * s * nH_m * (sjbz / ne_eff); // - a4 * rhoH_nd * wH * vst;
                        RHS_H(i, j, k, 4) += momentum_induce_coeff * nH_m * subu + momentum_hall_coeff * s * nH_m * (sjbu / ne_eff); // + a4 * rhoH_nd * us2 * b2; // work term for species energy
                    }

                    // =====================
                    // species Na+ (ls=2)
                    // =====================
                    {

                        const double subx = (vN - upy) * Bz - (wN - upz) * By;
                        const double suby = (wN - upz) * Bx - (uN - upx) * Bz;
                        const double subz = (uN - upx) * By - (vN - upy) * Bx;

                        const double sjbu = sjbx * uN + sjby * vN + sjbz * wN;
                        const double dpeu = dpex * uN + dpey * vN + dpez * wN;
                        const double subu = subx * uN + suby * vN + subz * wN;

                        const double Ts0 = PVN(i, j, k, 4) * T_ref;
                        const double us2 = uN * uN + vN * vN + wN * wN;
                        const double b1 = (Tn0 - Ts0) * (m_H / (m_H + m_Na));

                        // RHS_Na(i, j, k, 0) += a1_Na * sss; // Na+ mass creation
                        // RHS_Na(i, j, k, 1) += a2 * sns0 * subx + a3 * sns0 * (sjbx / ne_cm) - sns0 * (dpex / ne_cm) - a4 * rhoNa_nd * uN * vst;
                        // RHS_Na(i, j, k, 2) += a2 * sns0 * suby + a3 * sns0 * (sjby / ne_cm) - sns0 * (dpey / ne_cm) - a4 * rhoNa_nd * vN * vst;
                        // RHS_Na(i, j, k, 3) += a2 * sns0 * subz + a3 * sns0 * (sjbz / ne_cm) - sns0 * (dpez / ne_cm) - a4 * rhoNa_nd * wN * vst;
                        // RHS_Na(i, j, k, 4) += a2 * sns0 * subu + a3 * sns0 * (sjbu / ne_cm) - sns0 * (dpeu / ne_cm) + a4 * rhoNa_nd * us2 * vst;

                        // RHS_Na(i, j, k, 4) += a5 * sns0 * b1 + a6 * sns0 * sse * Tn0 / ne_cm + a6 * sss * Tn0;

                        // RHS_Na(i, j, k, 1) += -a1_Na * sss * uN; // 光致电力产生速度为零，相对流动的Na离子产生动量源项
                        // RHS_Na(i, j, k, 2) += -a1_Na * sss * vN; // 光致电力产生速度为零，相对流动的Na离子产生动量源项
                        // RHS_Na(i, j, k, 3) += -a1_Na * sss * wN; // 光致电力产生速度为零，相对流动的Na离子产生动量源项

                        // RHS_Na(i, j, k, 0) += a1_Na * sss; // Na+ mass creation
                        RHS_Na(i, j, k, 1) += momentum_induce_coeff * nNa_m * subx + momentum_hall_coeff * s * nNa_m * (sjbx / ne_eff); // - a4 * rhoNa_nd * uN * vst;
                        RHS_Na(i, j, k, 2) += momentum_induce_coeff * nNa_m * suby + momentum_hall_coeff * s * nNa_m * (sjby / ne_eff); // - a4 * rhoNa_nd * vN * vst;
                        RHS_Na(i, j, k, 3) += momentum_induce_coeff * nNa_m * subz + momentum_hall_coeff * s * nNa_m * (sjbz / ne_eff); // - a4 * rhoNa_nd * wN * vst;
                        RHS_Na(i, j, k, 4) += momentum_induce_coeff * nNa_m * subu + momentum_hall_coeff * s * nNa_m * (sjbu / ne_eff); // + a4 * rhoNa_nd * us2 * vst;
                    }

                    {
                        // --------- Debug: check charge-weighted consistency of u_plus ---------
                        // R = nH*(uH-u+) + nNa*(uNa-u+)  (should be ~0 if u_plus is charge-weighted)

                        const double duHx = uH - upx;
                        const double duHy = vH - upy;
                        const double duHz = wH - upz;

                        const double duNax = uN - upx;
                        const double duNay = vN - upy;
                        const double duNaz = wN - upz;

                        // Use the same number densities you already computed in m^-3
                        const double Rx = nH_m * duHx + nNa_m * duNax;
                        const double Ry = nH_m * duHy + nNa_m * duNay;
                        const double Rz = nH_m * duHz + nNa_m * duNaz;

                        const double Rabs = std::sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

                        // A relative metric: |R| / (ne * (|uH-u+|+|uNa-u+|) + eps)
                        const double duH_abs = std::sqrt(duHx * duHx + duHy * duHy + duHz * duHz);
                        const double duNa_abs = std::sqrt(duNax * duNax + duNay * duNay + duNaz * duNaz);
                        const double denom = ne_eff * (duH_abs + duNa_abs) + 1e-300;
                        const double Rrel = Rabs / denom;

                        // accumulate maxima (define these outside loops as double locals)
                        Rabs_max_l = std::max(Rabs_max_l, Rabs);
                        Rrel_max_l = std::max(Rrel_max_l, Rrel);

                        // optional: remember where it happens (use ints outside loops)
                        if (Rabs == Rabs_max_l)
                        {
                            Rmax_ib_l = ib;
                            Rmax_i_l = i;
                            Rmax_j_l = j;
                            Rmax_k_l = k;
                        }
                    }

                    {
                        // ===================== Debug: EM force consistency =====================
                        // Build per-species EM force you add via sources (only EM parts):
                        //   F_s = induce_coeff * n_s * ((u_s - u_plus) x B)
                        //       + hall_coeff   * n_s/ne * (J x B)
                        // Then compare sum_s F_s  vs  hall_coeff * (J x B).

                        // (1) already have subx/suby/subz for each species
                        // H+:
                        const double subx = (vH - upy) * Bz - (wH - upz) * By;
                        const double suby = (wH - upz) * Bx - (uH - upx) * Bz;
                        const double subz = (uH - upx) * By - (vH - upy) * Bx;
                        const double FHx = momentum_induce_coeff * nH_m * subx + momentum_hall_coeff * s * nH_m * (sjbx / ne_eff);
                        const double FHy = momentum_induce_coeff * nH_m * suby + momentum_hall_coeff * s * nH_m * (sjby / ne_eff);
                        const double FHz = momentum_induce_coeff * nH_m * subz + momentum_hall_coeff * s * nH_m * (sjbz / ne_eff);

                        // Na+:
                        const double subx_Na = (vN - upy) * Bz - (wN - upz) * By;
                        const double suby_Na = (wN - upz) * Bx - (uN - upx) * Bz;
                        const double subz_Na = (uN - upx) * By - (vN - upy) * Bx;
                        const double FNx = momentum_induce_coeff * nNa_m * subx_Na + momentum_hall_coeff * s * nNa_m * (sjbx / ne_eff);
                        const double FNy = momentum_induce_coeff * nNa_m * suby_Na + momentum_hall_coeff * s * nNa_m * (sjby / ne_eff);
                        const double FNz = momentum_induce_coeff * nNa_m * subz_Na + momentum_hall_coeff * s * nNa_m * (sjbz / ne_eff);

                        // sum
                        const double Fsx = FHx + FNx;
                        const double Fsy = FHy + FNy;
                        const double Fsz = FHz + FNz;

                        // target (total Lorentz force in your nondim form)
                        const double chi = s * (nH_m + nNa_m) / (ne_eff); // 或者 chi=s*(ne_true/ne_eff) 看你如何定义目标
                        const double FJx = momentum_hall_coeff * chi * sjbx;
                        const double FJy = momentum_hall_coeff * chi * sjby;
                        const double FJz = momentum_hall_coeff * chi * sjbz;

                        // absolute and relative errors
                        const double ex = Fsx - FJx;
                        const double ey = Fsy - FJy;
                        const double ez = Fsz - FJz;

                        const double eabs = std::sqrt(ex * ex + ey * ey + ez * ez);
                        const double FJabs = std::sqrt(FJx * FJx + FJy * FJy + FJz * FJz);

                        // relative: |e| / (|FJ| + eps)
                        const double erel = eabs / (FJabs + 1e-300);

                        // accumulate max (locals defined outside loops)
                        EM_eabs_max_l = std::max(EM_eabs_max_l, eabs);
                        EM_erel_max_l = std::max(EM_erel_max_l, erel);

                        // Optional: also track how big the induce-part sum is (should be ~0)
                        const double Findx = momentum_induce_coeff * (nH_m * subx + nNa_m * subx_Na);
                        const double Findy = momentum_induce_coeff * (nH_m * suby + nNa_m * suby_Na);
                        const double Findz = momentum_induce_coeff * (nH_m * subz + nNa_m * subz_Na);
                        const double Find_abs = std::sqrt(Findx * Findx + Findy * Findy + Findz * Findz);
                        EM_Find_abs_max_l = std::max(EM_Find_abs_max_l, Find_abs);
                    }
                }
    }

    double Rabs_max_g = Rabs_max_l;
    double Rrel_max_g = Rrel_max_l;

    // global maxima
    PARALLEL::mpi_max(&Rabs_max_l, &Rabs_max_g, 1);
    PARALLEL::mpi_max(&Rrel_max_l, &Rrel_max_g, 1);

    if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    {
        std::printf("[UplusCheck] step=%d  max|R|=%.3e  maxRel=%.3e\n",
                    run_data_->step, Rabs_max_g, Rrel_max_g);
        std::fflush(stdout);
    }

    double EM_eabs_max_g = EM_eabs_max_l;
    double EM_erel_max_g = EM_erel_max_l;
    double EM_Find_abs_max_g = EM_Find_abs_max_l;

    PARALLEL::mpi_max(&EM_eabs_max_l, &EM_eabs_max_g, 1);
    PARALLEL::mpi_max(&EM_erel_max_l, &EM_erel_max_g, 1);
    PARALLEL::mpi_max(&EM_Find_abs_max_l, &EM_Find_abs_max_g, 1);

    if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    {
        std::printf("[EMForceCheck] step=%d  max|Fsum-FJB|=%.3e  maxRel=%.3e  max|Find_sum|=%.3e\n",
                    run_data_->step, EM_eabs_max_g, EM_erel_max_g, EM_Find_abs_max_g);
        std::fflush(stdout);
    }
}