#include "MercurySolver.h"

void MercurySolver::calc_physical_constant(Param *par)
{
    // -------- constants / reference ----------
    auto cst = par_->GetDou_List("constant");
    auto ref = par_->GetDou_List("REF");

    gamma_ = cst.data["gamma"];
    double R_uni = cst.data["R_uni"];

    const double U_ref = ref.data.at("U");
    const double T_ref = ref.data.at("T");
    M_H = par->GetDou("mole_mass1");
    M_Na = par->GetDou("mole_mass2");

    // 无量纲状态方程系数：p = rho * T * coeff
    // coeff = (R_uni * T_ref) / (M * U_ref^2)
    state_coeff_H = (R_uni * T_ref) / (M_H * U_ref * U_ref);
    state_coeff_Na = (R_uni * T_ref) / (M_Na * U_ref * U_ref);

    CFL = par_->GetDou("CFL");

    hall_coef = ref.data["B_ref"] / (U_ref * cst.data["q_e"] * ref.data["L_ref"] * cst.data["mu_mag"]);

    rho_ref = M_H * ref.data["n"];

    ambi_coef = M_H * U_ref / (cst.data["q_e"] * ref.data["L_ref"] * ref.data["B_ref"]);

    inver_MA2 = ref.data["B_ref"] * ref.data["B_ref"] / (U_ref * U_ref * cst.data["mu_mag"] * rho_ref);
}

void MercurySolver::calc_PV()
{
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    auto fill_one = [&](int fidU, int fidPV, double coeff)
    {
        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            FieldBlock &U = fld_->field(fidU, ib);
            FieldBlock &PV = fld_->field(fidPV, ib);
            if (!U.is_allocated() || !PV.is_allocated())
                continue;

            const Int3 lo = PV.get_lo(); // 含 ghost：[-ng, ...]
            const Int3 hi = PV.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double rho0 = U(i, j, k, 0);
                        const double rho = std::max(rho0, rho_floor);

                        const double u = U(i, j, k, 1) / rho;
                        const double v = U(i, j, k, 2) / rho;
                        const double w = U(i, j, k, 3) / rho;

                        const double E = U(i, j, k, 4);
                        const double ke = 0.5 * rho * (u * u + v * v + w * w);

                        // 压力：p = (gamma-1) * (E - ke)
                        double eint = E - ke;
                        if (eint < 0.0)
                            eint = 0.0;
                        double p = (gamma_ - 1.0) * eint;
                        p = std::max(p, p_floor);

                        // 温度：T = p / (rho * coeff)
                        // 对 Na 用 coeff_Na -> 等价于 Fortran 的 T = p/(n kB)，其中 n= rho/m_Na
                        double T = p / (rho * coeff);

                        PV(i, j, k, 0) = u;
                        PV(i, j, k, 1) = v;
                        PV(i, j, k, 2) = w;
                        PV(i, j, k, 3) = p;
                        PV(i, j, k, 4) = T;
                    }
        }
    };

    fill_one(fid_.fid_U_H, fid_.fid_PV_H, state_coeff_H);
    fill_one(fid_.fid_U_Na, fid_.fid_PV_Na, state_coeff_Na);
}

void MercurySolver::calc_Uplus()
{
    const double rho_floor = 1e-20;
    const double inv23 = M_H / M_Na; // 与 Fortran sm2≈23*sm1 对齐

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib); // 新增

        if (!UH.is_allocated() || !UN.is_allocated() || !Up.is_allocated())
            continue;

        const Int3 lo = Up.get_lo();
        const Int3 hi = Up.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH = std::max(UH(i, j, k, 0), rho_floor);
                    const double rhoNa = std::max(UN(i, j, k, 0), 0.0);

                    const double uH = UH(i, j, k, 1) / rhoH;
                    const double vH = UH(i, j, k, 2) / rhoH;
                    const double wH = UH(i, j, k, 3) / rhoH;

                    double uNa = 0, vNa = 0, wNa = 0;
                    if (rhoNa > rho_floor)
                    {
                        uNa = UN(i, j, k, 1) / rhoNa;
                        vNa = UN(i, j, k, 2) / rhoNa;
                        wNa = UN(i, j, k, 3) / rhoNa;
                    }

                    const double nH = rhoH;
                    const double nNa = rhoNa * inv23;
                    const double nt = nH + nNa;

                    if (nt <= rho_floor)
                    {
                        Up(i, j, k, 0) = 0.0;
                        Up(i, j, k, 1) = 0.0;
                        Up(i, j, k, 2) = 0.0;
                    }
                    else
                    {
                        Up(i, j, k, 0) = (nH * uH + nNa * uNa) / nt;
                        Up(i, j, k, 1) = (nH * vH + nNa * vNa) / nt;
                        Up(i, j, k, 2) = (nH * wH + nNa * wNa) / nt;
                    }
                }
    }
}

void MercurySolver::calc_Bcell()
{
    const int nblock = fld_->num_blocks();

    const double eps = 1e-300;
    const double delta = 1e-15; // same spirit as your reference

    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto norm = [&](const std::array<double, 3> &a)
    {
        return std::sqrt(dot(a, a));
    };

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
        // auto &U = fld_->field(fid_.fid_U, ib);

        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);

        auto &Baddxi = fld_->field(fid_.fid_Badd.xi, ib);
        auto &Baddeta = fld_->field(fid_.fid_Badd.eta, ib);
        auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);

        auto &Jac = fld_->field(fid_.fid_Jac, ib);
        auto &A_xi = fld_->field(fid_.fid_metric.xi, ib);   // JDxi
        auto &A_eta = fld_->field(fid_.fid_metric.eta, ib); // JDet
        auto &A_ze = fld_->field(fid_.fid_metric.zeta, ib); // JDze

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        auto &cx = grd_->grids(ib).dual_x;
        auto &cy = grd_->grids(ib).dual_y;
        auto &cz = grd_->grids(ib).dual_z;

        // face centers (consistent with your Jac construction style)
        auto xfc_xi = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {
                0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
                0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
                0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};
        };
        auto xfc_eta = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {
                0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
                0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
                0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};
        };
        auto xfc_zeta = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {
                0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
                0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
                0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};
        };

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // cell center
                    std::array<double, 3> Xc = {
                        cx(i + 1, j + 1, k + 1),
                        cy(i + 1, j + 1, k + 1),
                        cz(i + 1, j + 1, k + 1)};

                    struct FaceEq
                    {
                        std::array<double, 3> n; // normalized S
                        double phi;              // normalized Phi
                        double w;                // weight
                    };
                    FaceEq eqs[6];
                    int K = 0;

                    auto push = [&](const std::array<double, 3> &S,
                                    double Phi,
                                    const std::array<double, 3> &Xf)
                    {
                        double s_norm = norm(S) + eps;
                        std::array<double, 3> nvec = {S[0] / s_norm, S[1] / s_norm, S[2] / s_norm};
                        double phi_hat = Phi / s_norm;

                        double dx = Xf[0] - Xc[0];
                        double dy = Xf[1] - Xc[1];
                        double dz = Xf[2] - Xc[2];
                        double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                        double w = 1.0 / std::sqrt(dist * dist + delta * delta);

                        eqs[K++] = {nvec, phi_hat, w};
                    };

                    // ---------- xi- face (at i) ----------
                    std::array<double, 3> S_xm = {
                        -A_xi(i, j, k, 0),
                        -A_xi(i, j, k, 1),
                        -A_xi(i, j, k, 2)};
                    double Phi_xm = -Bxi(i, j, k, 0) - Baddxi(i, j, k, 0);
                    push(S_xm, Phi_xm, xfc_xi(i, j, k));

                    // ---------- xi+ face (at i+1) ----------
                    std::array<double, 3> S_xp = {
                        A_xi(i + 1, j, k, 0),
                        A_xi(i + 1, j, k, 1),
                        A_xi(i + 1, j, k, 2)};
                    double Phi_xp = Bxi(i + 1, j, k, 0) + Baddxi(i + 1, j, k, 0);
                    push(S_xp, Phi_xp, xfc_xi(i + 1, j, k));

                    // ---------- eta- face (at j) ----------
                    std::array<double, 3> S_em = {
                        -A_eta(i, j, k, 0),
                        -A_eta(i, j, k, 1),
                        -A_eta(i, j, k, 2)};
                    double Phi_em = -Beta(i, j, k, 0) - Baddeta(i, j, k, 0);
                    push(S_em, Phi_em, xfc_eta(i, j, k));

                    // ---------- eta+ face (at j+1) ----------
                    std::array<double, 3> S_ep = {
                        A_eta(i, j + 1, k, 0),
                        A_eta(i, j + 1, k, 1),
                        A_eta(i, j + 1, k, 2)};
                    double Phi_ep = Beta(i, j + 1, k, 0) + Baddeta(i, j + 1, k, 0);
                    push(S_ep, Phi_ep, xfc_eta(i, j + 1, k));

                    // ---------- zeta- face (at k) ----------
                    std::array<double, 3> S_zm = {
                        -A_ze(i, j, k, 0),
                        -A_ze(i, j, k, 1),
                        -A_ze(i, j, k, 2)};
                    double Phi_zm = -Bzeta(i, j, k, 0) - Baddzeta(i, j, k, 0);
                    push(S_zm, Phi_zm, xfc_zeta(i, j, k));

                    // ---------- zeta+ face (at k+1) ----------
                    std::array<double, 3> S_zp = {
                        A_ze(i, j, k + 1, 0),
                        A_ze(i, j, k + 1, 1),
                        A_ze(i, j, k + 1, 2)};
                    double Phi_zp = Bzeta(i, j, k + 1, 0) + Baddzeta(i, j, k + 1, 0);
                    push(S_zp, Phi_zp, xfc_zeta(i, j, k + 1));

                    // ---------- build normal equations N = A^T W A, r = A^T W phi ----------
                    double N00 = 0, N01 = 0, N02 = 0, N11 = 0, N12 = 0, N22 = 0;
                    double rx = 0, ry = 0, rz = 0;

                    for (int t = 0; t < K; ++t)
                    {
                        double w = eqs[t].w;
                        const auto &n = eqs[t].n;
                        double phi = eqs[t].phi;

                        N00 += w * n[0] * n[0];
                        N01 += w * n[0] * n[1];
                        N02 += w * n[0] * n[2];
                        N11 += w * n[1] * n[1];
                        N12 += w * n[1] * n[2];
                        N22 += w * n[2] * n[2];

                        rx += w * phi * n[0];
                        ry += w * phi * n[1];
                        rz += w * phi * n[2];
                    }

                    auto det3 = [&](double a, double b, double c, double d, double e, double f)
                    {
                        // | a b c |
                        // | b d e |
                        // | c e f |
                        return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
                    };

                    double det = det3(N00, N01, N02, N11, N12, N22);
                    double reg = 1e-14 * (N00 + N11 + N22);

                    if (std::abs(det) < reg)
                    {
                        N00 += reg;
                        N11 += reg;
                        N22 += reg;
                        det = det3(N00, N01, N02, N11, N12, N22);
                    }

                    // cofactors of symmetric matrix
                    double C00 = (N11 * N22 - N12 * N12);
                    double C01 = (N02 * N12 - N01 * N22);
                    double C02 = (N01 * N12 - N02 * N11);
                    double C11 = (N00 * N22 - N02 * N02);
                    double C12 = (N01 * N02 - N00 * N12);
                    double C22 = (N00 * N11 - N01 * N01);

                    double inv = 1.0 / det;

                    double Bx_tot = inv * (C00 * rx + C01 * ry + C02 * rz);
                    double By_tot = inv * (C01 * rx + C11 * ry + C12 * rz);
                    double Bz_tot = inv * (C02 * rx + C12 * ry + C22 * rz);

                    // // energy consistency (keep your existing style)
                    // const double Bx_old = Bcell(i, j, k, 0);
                    // const double By_old = Bcell(i, j, k, 1);
                    // const double Bz_old = Bcell(i, j, k, 2);
                    // const double Delta_Eb =
                    //     0.5 * inver_MA2 * (Bx_tot * Bx_tot + By_tot * By_tot + Bz_tot * Bz_tot) -
                    //     0.5 * inver_MA2 * (Bx_old * Bx_old + By_old * By_old + Bz_old * Bz_old);

                    Bcell(i, j, k, 0) = Bx_tot;
                    Bcell(i, j, k, 1) = By_tot;
                    Bcell(i, j, k, 2) = Bz_tot;
                }
    }
}

void MercurySolver::calc_divB()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &divB = fld_->field(fid_.fid_divB, ib);
        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);
        auto &Jac = fld_->field(fid_.fid_Jac, ib);

        Int3 lo = divB.inner_lo();
        Int3 hi = divB.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    divB(i, j, k, 0) = (Bxi(i + 1, j, k, 0) - Bxi(i, j, k, 0) +
                                        Beta(i, j + 1, k, 0) - Beta(i, j, k, 0) +
                                        Bzeta(i, j, k + 1, 0) - Bzeta(i, j, k, 0)) /
                                       Jac(i, j, k, 0);
                }
    }
}
