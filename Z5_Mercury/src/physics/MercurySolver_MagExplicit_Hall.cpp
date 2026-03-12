#include "MercurySolver.h"
void MercurySolver::AddHallEdgeEMF_()
{
    // 0) 清零 Ehall_edge / Ehall_face
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto zero_scalar_edge = [&](FieldBlock &F)
        {
            if (!F.is_allocated())
                return;
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        F(i, j, k, 0) = 0.0;
        };

        auto zero_vec_face = [&](FieldBlock &F)
        {
            if (!F.is_allocated())
                return;
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 3; ++m)
                            F(i, j, k, m) = 0.0;
        };

        zero_scalar_edge(fld_->field(fid_.fid_Ehall.xi, ib));
        zero_scalar_edge(fld_->field(fid_.fid_Ehall.eta, ib));
        zero_scalar_edge(fld_->field(fid_.fid_Ehall.zeta, ib));

        zero_vec_face(fld_->field(fid_.fid_Eface.xi, ib));
        zero_vec_face(fld_->field(fid_.fid_Eface.eta, ib));
        zero_vec_face(fld_->field(fid_.fid_Eface.zeta, ib));
    }

    // 2) face 上做 Hall-Rusanov
    BuildHallFaceEMF_Rusanov_();

    // 3) sync Hall face field
    mercury_bound_.Sync("Eface");

    // 4) face -> edge
    AssembleEdgeEMF_FromFaceE_Hall_();

    // 5) sync edge Hall EMF
    mercury_bound_.Sync("Ehall");

    // 6) 加到总 E 上
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &E_xi = fld_->field(fid_.fid_E.xi, ib);
        auto &E_et = fld_->field(fid_.fid_E.eta, ib);
        auto &E_ze = fld_->field(fid_.fid_E.zeta, ib);

        auto &EH_xi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &EH_et = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &EH_ze = fld_->field(fid_.fid_Ehall.zeta, ib);

        if (E_xi.is_allocated())
        {
            Int3 lo = E_xi.inner_lo(), hi = E_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_xi(i, j, k, 0) += EH_xi(i, j, k, 0);
        }

        if (E_et.is_allocated())
        {
            Int3 lo = E_et.inner_lo(), hi = E_et.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_et(i, j, k, 0) += EH_et(i, j, k, 0);
        }

        if (E_ze.is_allocated())
        {
            Int3 lo = E_ze.inner_lo(), hi = E_ze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_ze(i, j, k, 0) += EH_ze(i, j, k, 0);
        }
    }
}
// void MercurySolver::AddHallEdgeEMF_()
// {

//     const int nb = fld_->num_blocks();

//     const double hall_coeff = hall_coef;
//     constexpr double eps = 1e-30;

//     // 这里先把人工 Hall 电阻关掉，先只看纯 JxB
//     constexpr double C_eta = 0.0;

//     auto cross_x = [&](double Jx, double Jy, double Jz,
//                        double Bx, double By, double Bz) -> double
//     {
//         return Jy * Bz - Jz * By;
//     };
//     auto cross_y = [&](double Jx, double Jy, double Jz,
//                        double Bx, double By, double Bz) -> double
//     {
//         return Jz * Bx - Jx * Bz;
//     };
//     auto cross_z = [&](double Jx, double Jy, double Jz,
//                        double Bx, double By, double Bz) -> double
//     {
//         return Jx * By - Jy * Bx;
//     };

//     for (int ib = 0; ib < nb; ++ib)
//     {
//         auto &UH = fld_->field(fid_.fid_U_H, ib);
//         auto &UNa = fld_->field(fid_.fid_U_Na, ib);

//         auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
//         auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

//         auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, ib);
//         auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, ib);
//         auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, ib);

//         auto &E_xi = fld_->field(fid_.fid_E.xi, ib);
//         auto &E_eta = fld_->field(fid_.fid_E.eta, ib);
//         auto &E_zeta = fld_->field(fid_.fid_E.zeta, ib);

//         if (!UH.is_allocated() || !UNa.is_allocated() ||
//             !Bcell.is_allocated() || !Jcell.is_allocated() ||
//             !Ehall_xi.is_allocated() || !Ehall_eta.is_allocated() || !Ehall_zeta.is_allocated())
//             continue;

//         auto &x = grd_->grids(ib).x;
//         auto &y = grd_->grids(ib).y;
//         auto &z = grd_->grids(ib).z;
//         auto &cx = grd_->grids(ib).dual_x;
//         auto &cy = grd_->grids(ib).dual_y;
//         auto &cz = grd_->grids(ib).dual_z;

//         // ------------------------------------------------------------
//         // 1) clear Ehall
//         // ------------------------------------------------------------
//         auto zero_inner_scalar = [&](FieldBlock &F)
//         {
//             Int3 lo = F.inner_lo();
//             Int3 hi = F.inner_hi();
//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         F(i, j, k, 0) = 0.0;
//         };

//         zero_inner_scalar(Ehall_xi);
//         zero_inner_scalar(Ehall_eta);
//         zero_inner_scalar(Ehall_zeta);

//         // ------------------------------------------------------------
//         // 2) cell Hall coefficient alpha(i,j,k)
//         // ------------------------------------------------------------
//         auto NumCell = [&](int i, int j, int k) -> double
//         {
//             return UH(i, j, k, 0) / M_H + UNa(i, j, k, 0) / M_Na;
//         };

//         auto radius_cell = [&](int i, int j, int k) -> double
//         {
//             // cell(i,j,k) center -> dual(i+1,j+1,k+1)
//             const double xx = cx(i + 1, j + 1, k + 1);
//             const double yy = cy(i + 1, j + 1, k + 1);
//             const double zz = cz(i + 1, j + 1, k + 1);
//             return std::sqrt(xx * xx + yy * yy + zz * zz);
//         };

//         auto hall_alpha_from_ne = [&](double Num, double r) -> double
//         {
//             const double ne_true = std::max(Num, 0.0);

//             // 平滑 floor
//             const double ne_eff = std::sqrt(ne_true * ne_true +
//                                             ne_hall_floor * ne_hall_floor);

//             // 平滑 cut
//             const double s_ne = ne_true / (ne_true + ne_hall_cut);

//             // 你原来那套径向 taper 先保留
//             if (r <= 1.01)
//                 return 0.0;
//             if (r >= 1.50)
//                 return hall_coeff * s_ne / ne_eff;

//             const double xi = (r - 1.01) / 0.49;
//             const double w = xi * xi * (3.0 - 2.0 * xi); // smoothstep

//             return hall_coeff * s_ne / ne_eff * w;
//         };

//         // ------------------------------------------------------------
//         // 3) evaluate cell Hall electric field vector
//         //    Ecell = alpha * (Jcell x Bcell) + eta_h * Jcell
//         // ------------------------------------------------------------
//         auto eval_Ehall_cell = [&](int i, int j, int k,
//                                    double &Ex, double &Ey, double &Ez)
//         {
//             const double Jx = Jcell(i, j, k, 0);
//             const double Jy = Jcell(i, j, k, 1);
//             const double Jz = Jcell(i, j, k, 2);

//             const double Bx = Bcell(i, j, k, 0);
//             const double By = Bcell(i, j, k, 1);
//             const double Bz = Bcell(i, j, k, 2);

//             const double Num = NumCell(i, j, k);
//             const double rcell = radius_cell(i, j, k);
//             const double alpha = hall_alpha_from_ne(Num, rcell);

//             Ex = alpha * cross_x(Jx, Jy, Jz, Bx, By, Bz) + C_eta * Jx;
//             Ey = alpha * cross_y(Jx, Jy, Jz, Bx, By, Bz) + C_eta * Jy;
//             Ez = alpha * cross_z(Jx, Jy, Jz, Bx, By, Bz) + C_eta * Jz;
//         };

//         // ============================================================
//         // EdgeXi:
//         // surrounding cells = (i,j,k), (i,j-1,k), (i,j,k-1), (i,j-1,k-1)
//         // ============================================================
//         {
//             Int3 lo = Ehall_xi.inner_lo();
//             Int3 hi = Ehall_xi.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         double E0x, E0y, E0z;
//                         double E1x, E1y, E1z;
//                         double E2x, E2y, E2z;
//                         double E3x, E3y, E3z;

//                         eval_Ehall_cell(i, j, k, E0x, E0y, E0z);
//                         eval_Ehall_cell(i, j - 1, k, E1x, E1y, E1z);
//                         eval_Ehall_cell(i, j, k - 1, E2x, E2y, E2z);
//                         eval_Ehall_cell(i, j - 1, k - 1, E3x, E3y, E3z);

//                         const double Ex = 0.25 * (E0x + E1x + E2x + E3x);
//                         const double Ey = 0.25 * (E0y + E1y + E2y + E3y);
//                         const double Ez = 0.25 * (E0z + E1z + E2z + E3z);

//                         const double dx = x(i + 1, j, k) - x(i, j, k);
//                         const double dy = y(i + 1, j, k) - y(i, j, k);
//                         const double dz = z(i + 1, j, k) - z(i, j, k);

//                         Ehall_xi(i, j, k, 0) = Ex * dx + Ey * dy + Ez * dz;
//                     }
//         }

//         // ============================================================
//         // EdgeEta:
//         // surrounding cells = (i,j,k), (i-1,j,k), (i,j,k-1), (i-1,j,k-1)
//         // ============================================================
//         {
//             Int3 lo = Ehall_eta.inner_lo();
//             Int3 hi = Ehall_eta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         double E0x, E0y, E0z;
//                         double E1x, E1y, E1z;
//                         double E2x, E2y, E2z;
//                         double E3x, E3y, E3z;

//                         eval_Ehall_cell(i, j, k, E0x, E0y, E0z);
//                         eval_Ehall_cell(i - 1, j, k, E1x, E1y, E1z);
//                         eval_Ehall_cell(i, j, k - 1, E2x, E2y, E2z);
//                         eval_Ehall_cell(i - 1, j, k - 1, E3x, E3y, E3z);

//                         const double Ex = 0.25 * (E0x + E1x + E2x + E3x);
//                         const double Ey = 0.25 * (E0y + E1y + E2y + E3y);
//                         const double Ez = 0.25 * (E0z + E1z + E2z + E3z);

//                         const double dx = x(i, j + 1, k) - x(i, j, k);
//                         const double dy = y(i, j + 1, k) - y(i, j, k);
//                         const double dz = z(i, j + 1, k) - z(i, j, k);

//                         Ehall_eta(i, j, k, 0) = Ex * dx + Ey * dy + Ez * dz;
//                     }
//         }

//         // ============================================================
//         // EdgeZeta:
//         // surrounding cells = (i,j,k), (i-1,j,k), (i,j-1,k), (i-1,j-1,k)
//         // ============================================================
//         {
//             Int3 lo = Ehall_zeta.inner_lo();
//             Int3 hi = Ehall_zeta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         double E0x, E0y, E0z;
//                         double E1x, E1y, E1z;
//                         double E2x, E2y, E2z;
//                         double E3x, E3y, E3z;

//                         eval_Ehall_cell(i, j, k, E0x, E0y, E0z);
//                         eval_Ehall_cell(i - 1, j, k, E1x, E1y, E1z);
//                         eval_Ehall_cell(i, j - 1, k, E2x, E2y, E2z);
//                         eval_Ehall_cell(i - 1, j - 1, k, E3x, E3y, E3z);

//                         const double Ex = 0.25 * (E0x + E1x + E2x + E3x);
//                         const double Ey = 0.25 * (E0y + E1y + E2y + E3y);
//                         const double Ez = 0.25 * (E0z + E1z + E2z + E3z);

//                         const double dx = x(i, j, k + 1) - x(i, j, k);
//                         const double dy = y(i, j, k + 1) - y(i, j, k);
//                         const double dz = z(i, j, k + 1) - z(i, j, k);

//                         Ehall_zeta(i, j, k, 0) = Ex * dx + Ey * dy + Ez * dz;
//                     }
//         }

//         // 如果你已经给 Ehall 建了边界/halo pattern，可以保留
//         mercury_bound_.Sync("Ehall");

//         // ------------------------------------------------------------
//         // 4) add Ehall into total edge EMF
//         // ------------------------------------------------------------
//         auto add_inner_scalar = [&](FieldBlock &E, FieldBlock &Eh)
//         {
//             Int3 lo = E.inner_lo();
//             Int3 hi = E.inner_hi();
//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         E(i, j, k, 0) += Eh(i, j, k, 0);
//         };

//         add_inner_scalar(E_xi, Ehall_xi);
//         add_inner_scalar(E_eta, Ehall_eta);
//         add_inner_scalar(E_zeta, Ehall_zeta);
//     }
// }

void MercurySolver::BuildHallFaceEMF_Rusanov_()
{
    constexpr double eps = 1e-300;
    constexpr double Cwh = 1.0;   // 先取 1.0，后面可调 0.5~2
    constexpr double C_eta = 0.0; // 先关掉额外人工耗散

    auto cross3 = [](const std::array<double, 3> &a,
                     const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
    };

    auto dot3 = [](const std::array<double, 3> &a,
                   const std::array<double, 3> &b) -> double
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    auto norm3 = [&](const std::array<double, 3> &a) -> double
    {
        return std::sqrt(dot3(a, a));
    };

    auto plus3 = [](const std::array<double, 3> &a,
                    const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };

    auto minus3 = [](const std::array<double, 3> &a,
                     const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };

    auto scale3 = [](const std::array<double, 3> &a, double s) -> std::array<double, 3>
    {
        return {s * a[0], s * a[1], s * a[2]};
    };

    auto proj_tangent = [&](const std::array<double, 3> &E,
                            const std::array<double, 3> &n) -> std::array<double, 3>
    {
        double En = dot3(E, n);
        return {
            E[0] - En * n[0],
            E[1] - En * n[1],
            E[2] - En * n[2]};
    };

    const int nb = fld_->num_blocks();

    auto hall_alpha_from_ne = [&](double ne_true, double r) -> double
    {
        const double ne_pos = std::max(ne_true, 0.0);
        const double ne_eff = std::sqrt(ne_pos * ne_pos + ne_hall_floor * ne_hall_floor);
        const double s_ne = ne_pos / (ne_pos + ne_hall_cut + eps);

        if (r <= 1.01)
            return 0.0;
        if (r >= 1.50)
            return hall_coef * s_ne / ne_eff;

        const double xi = (r - 1.01) / 0.49;
        const double w = xi * xi * (3.0 - 2.0 * xi); // smoothstep
        return hall_coef * s_ne / ne_eff * w;
    };

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);
        auto &Jc = fld_->field("J_cell", ib);               // 3 comps
        auto &Bc = fld_->field(fid_.fid_Bcell, ib);         // 3 comps, total B
        auto &Binduce = fld_->field(fid_.fid_Bindcell, ib); // 3 comps, total B

        auto &JDxi = fld_->field("JDxi", ib);
        auto &JDet = fld_->field("JDet", ib);
        auto &JDze = fld_->field("JDze", ib);

        auto &dlst_xi = fld_->field("dlstar_xi", ib);
        auto &dlst_et = fld_->field("dlstar_eta", ib);
        auto &dlst_ze = fld_->field("dlstar_zeta", ib);

        auto &Efxi = fld_->field(fid_.fid_Eface.xi, ib);
        auto &Efet = fld_->field(fid_.fid_Eface.eta, ib);
        auto &Efze = fld_->field(fid_.fid_Eface.zeta, ib);

        auto &cx = grd_->grids(ib).dual_x;
        auto &cy = grd_->grids(ib).dual_y;
        auto &cz = grd_->grids(ib).dual_z;

        if (!UH.is_allocated() || !UNa.is_allocated() || !Jc.is_allocated() ||
            !Bc.is_allocated() || !Efxi.is_allocated())
            continue;

        auto ne_cell = [&](int i, int j, int k) -> double
        {
            return UH(i, j, k, 0) / M_H + UNa(i, j, k, 0) / M_Na;
        };

        auto r_cell = [&](int i, int j, int k) -> double
        {
            double x0 = cx(i + 1, j + 1, k + 1);
            double y0 = cy(i + 1, j + 1, k + 1);
            double z0 = cz(i + 1, j + 1, k + 1);
            return std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
        };

        auto B_cell = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {Bc(i, j, k, 0), Bc(i, j, k, 1), Bc(i, j, k, 2)};
        };

        auto Bind_cell = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {Binduce(i, j, k, 0), Binduce(i, j, k, 1), Binduce(i, j, k, 2)};
        };

        auto J_cell = [&](int i, int j, int k) -> std::array<double, 3>
        {
            return {Jc(i, j, k, 0), Jc(i, j, k, 1), Jc(i, j, k, 2)};
        };

        auto Ehall_cell = [&](int i, int j, int k) -> std::array<double, 3>
        {
            std::array<double, 3> J = J_cell(i, j, k);
            std::array<double, 3> B = B_cell(i, j, k);
            double ne = ne_cell(i, j, k);
            double rr = r_cell(i, j, k);
            double alpha = hall_alpha_from_ne(ne, rr);

            auto E = cross3(J, B);
            E = scale3(E, alpha);

            if constexpr (true)
            {
                if (C_eta != 0.0)
                {
                    double Bmag = norm3(B);
                    double hloc = 1.0;
                    E = plus3(E, scale3(J, C_eta * Bmag * hloc));
                }
            }
            return E;
        };

        // --------------------------------------------------------
        // xi-face : face(i,j,k) separates cell(i-1,j,k) and cell(i,j,k)
        // --------------------------------------------------------
        {
            Int3 lo = Efxi.inner_lo();
            Int3 hi = Efxi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        int iL = i - 1, iR = i;

                        std::array<double, 3> BL = Bind_cell(iL, j, k);
                        std::array<double, 3> BR = Bind_cell(iR, j, k);
                        std::array<double, 3> BL_all = B_cell(iL, j, k);
                        std::array<double, 3> BR_all = B_cell(iR, j, k);
                        std::array<double, 3> EL = Ehall_cell(iL, j, k);
                        std::array<double, 3> ER = Ehall_cell(iR, j, k);

                        std::array<double, 3> S = {
                            JDxi(i, j, k, 0), JDxi(i, j, k, 1), JDxi(i, j, k, 2)};
                        double Smag = norm3(S) + eps;
                        std::array<double, 3> n = {S[0] / Smag, S[1] / Smag, S[2] / Smag};

                        EL = proj_tangent(EL, n);
                        ER = proj_tangent(ER, n);

                        double neL = ne_cell(iL, j, k);
                        double neR = ne_cell(iR, j, k);
                        double aL = std::abs(hall_alpha_from_ne(neL, r_cell(iL, j, k)));
                        double aR = std::abs(hall_alpha_from_ne(neR, r_cell(iR, j, k)));

                        double h_n = std::max(dlst_xi(i, j, k, 0), eps);
                        double sH = Cwh * std::max(aL * norm3(BL_all), aR * norm3(BR_all)) / h_n;

                        std::array<double, 3> Ecen = scale3(plus3(EL, ER), 0.5);
                        std::array<double, 3> dB = minus3(BR, BL);
                        std::array<double, 3> Ediss = scale3(cross3(n, dB), 0.5 * sH);
                        std::array<double, 3> Ef = plus3(Ecen, Ediss);

                        Efxi(i, j, k, 0) = Ef[0];
                        Efxi(i, j, k, 1) = Ef[1];
                        Efxi(i, j, k, 2) = Ef[2];
                    }
        }

        // --------------------------------------------------------
        // eta-face : face(i,j,k) separates cell(i,j-1,k) and cell(i,j,k)
        // --------------------------------------------------------
        {
            Int3 lo = Efet.inner_lo();
            Int3 hi = Efet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        int jL = j - 1, jR = j;

                        std::array<double, 3> BL = Bind_cell(i, jL, k);
                        std::array<double, 3> BR = Bind_cell(i, jR, k);
                        std::array<double, 3> BL_all = B_cell(i, jL, k);
                        std::array<double, 3> BR_all = B_cell(i, jR, k);
                        std::array<double, 3> EL = Ehall_cell(i, jL, k);
                        std::array<double, 3> ER = Ehall_cell(i, jR, k);

                        std::array<double, 3> S = {
                            JDet(i, j, k, 0), JDet(i, j, k, 1), JDet(i, j, k, 2)};
                        double Smag = norm3(S) + eps;
                        std::array<double, 3> n = {S[0] / Smag, S[1] / Smag, S[2] / Smag};

                        EL = proj_tangent(EL, n);
                        ER = proj_tangent(ER, n);

                        double neL = ne_cell(i, jL, k);
                        double neR = ne_cell(i, jR, k);
                        double aL = std::abs(hall_alpha_from_ne(neL, r_cell(i, jL, k)));
                        double aR = std::abs(hall_alpha_from_ne(neR, r_cell(i, jR, k)));

                        double h_n = std::max(dlst_et(i, j, k, 0), eps);
                        double sH = Cwh * std::max(aL * norm3(BL_all), aR * norm3(BR_all)) / h_n;

                        std::array<double, 3> Ecen = scale3(plus3(EL, ER), 0.5);
                        std::array<double, 3> dB = minus3(BR, BL);
                        std::array<double, 3> Ediss = scale3(cross3(n, dB), 0.5 * sH);
                        std::array<double, 3> Ef = plus3(Ecen, Ediss);

                        Efet(i, j, k, 0) = Ef[0];
                        Efet(i, j, k, 1) = Ef[1];
                        Efet(i, j, k, 2) = Ef[2];
                    }
        }

        // --------------------------------------------------------
        // zeta-face : face(i,j,k) separates cell(i,j,k-1) and cell(i,j,k)
        // --------------------------------------------------------
        {
            Int3 lo = Efze.inner_lo();
            Int3 hi = Efze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        int kL = k - 1, kR = k;

                        std::array<double, 3> BL = Bind_cell(i, j, kL);
                        std::array<double, 3> BR = Bind_cell(i, j, kR);
                        std::array<double, 3> BL_all = B_cell(i, j, kL);
                        std::array<double, 3> BR_all = B_cell(i, j, kR);
                        std::array<double, 3> EL = Ehall_cell(i, j, kL);
                        std::array<double, 3> ER = Ehall_cell(i, j, kR);

                        std::array<double, 3> S = {
                            JDze(i, j, k, 0), JDze(i, j, k, 1), JDze(i, j, k, 2)};
                        double Smag = norm3(S) + eps;
                        std::array<double, 3> n = {S[0] / Smag, S[1] / Smag, S[2] / Smag};

                        EL = proj_tangent(EL, n);
                        ER = proj_tangent(ER, n);

                        double neL = ne_cell(i, j, kL);
                        double neR = ne_cell(i, j, kR);
                        double aL = std::abs(hall_alpha_from_ne(neL, r_cell(i, j, kL)));
                        double aR = std::abs(hall_alpha_from_ne(neR, r_cell(i, j, kR)));

                        double h_n = std::max(dlst_ze(i, j, k, 0), eps);
                        double sH = Cwh * std::max(aL * norm3(BL_all), aR * norm3(BR_all)) / h_n;

                        std::array<double, 3> Ecen = scale3(plus3(EL, ER), 0.5);
                        std::array<double, 3> dB = minus3(BR, BL);
                        std::array<double, 3> Ediss = scale3(cross3(n, dB), 0.5 * sH);
                        std::array<double, 3> Ef = plus3(Ecen, Ediss);

                        Efze(i, j, k, 0) = Ef[0];
                        Efze(i, j, k, 1) = Ef[1];
                        Efze(i, j, k, 2) = Ef[2];
                    }
        }
    }
}

void MercurySolver::AssembleEdgeEMF_FromFaceE_Hall_()
{
    Int3 sub, sup;

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        double3D &x = fld_->grd->grids(iblk).x;
        double3D &y = fld_->grd->grids(iblk).y;
        double3D &z = fld_->grd->grids(iblk).z;

        Vec3 E, dr;

        // -----------------------------------------
        // Exi(edge) from face-eta and face-zeta
        // -----------------------------------------
        {
            auto &Exi = fld_->field(fid_.fid_Ehall.xi, iblk);
            auto &Ef_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_eta(i, j, k, 0) + Ef_eta(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i, j - 1, k, 0));
                        E.vec[1] = 0.25 * (Ef_eta(i, j, k, 1) + Ef_eta(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i, j - 1, k, 1));
                        E.vec[2] = 0.25 * (Ef_eta(i, j, k, 2) + Ef_eta(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i, j - 1, k, 2));

                        dr = {
                            x(i + 1, j, k) - x(i, j, k),
                            y(i + 1, j, k) - y(i, j, k),
                            z(i + 1, j, k) - z(i, j, k)};

                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Eeta(edge) from face-xi and face-zeta
        // -----------------------------------------
        {
            auto &Eeta = fld_->field(fid_.fid_Ehall.eta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i - 1, j, k, 2));

                        dr = {
                            x(i, j + 1, k) - x(i, j, k),
                            y(i, j + 1, k) - y(i, j, k),
                            z(i, j + 1, k) - z(i, j, k)};

                        Eeta(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Ezeta(edge) from face-xi and face-eta
        // -----------------------------------------
        {
            auto &Eze = fld_->field(fid_.fid_Ehall.zeta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_et = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Eze.inner_lo();
            sup = Eze.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j - 1, k, 0) + Ef_et(i, j, k, 0) + Ef_et(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j - 1, k, 1) + Ef_et(i, j, k, 1) + Ef_et(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j - 1, k, 2) + Ef_et(i, j, k, 2) + Ef_et(i - 1, j, k, 2));

                        dr = {
                            x(i, j, k + 1) - x(i, j, k),
                            y(i, j, k + 1) - y(i, j, k),
                            z(i, j, k + 1) - z(i, j, k)};

                        Eze(i, j, k, 0) = E * dr;
                    }
        }
    }
}
// void MercurySolver::AddHallEdgeEMF_()
// {
//     // Add Ehall_xi/eta/zeta（integration along Edge）to E
//     // ------------------------------------------------------------
//     // Hall coefficient:
//     //   E_hall = alpha * (J x B)
//     //   alpha typically = +1/(n e)  (sign can be absorbed here)
//     // ------------------------------------------------------------
//     const double hall_coeff = hall_coef;
//     // const double N_floor = ne_hall_floor; // 1e-300; // 防止除零
//     constexpr double C_eta = 0.0;

//     // Set E_hall to zero
//     for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
//     {
//         auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
//         auto &Ehall_et = fld_->field(fid_.fid_Ehall.eta, iblk);
//         auto &Ehall_ze = fld_->field(fid_.fid_Ehall.zeta, iblk);
//         if (!Ehall_xi.is_allocated() || !Ehall_et.is_allocated() || !Ehall_ze.is_allocated())
//             continue;

//         {
//             Int3 lo = Ehall_xi.inner_lo();
//             Int3 hi = Ehall_xi.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         Ehall_xi(i, j, k, 0) = 0.0;
//         }
//         {
//             Int3 lo = Ehall_et.inner_lo();
//             Int3 hi = Ehall_et.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         Ehall_et(i, j, k, 0) = 0.0;
//         }
//         {
//             Int3 lo = Ehall_ze.inner_lo();
//             Int3 hi = Ehall_ze.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         Ehall_ze(i, j, k, 0) = 0.0;
//         }
//     }

//     // Calculate E_hall
//     for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
//     {
//         auto &UH = fld_->field(fid_.fid_U_H, iblk);
//         auto &UNa = fld_->field(fid_.fid_U_Na, iblk);
//         if (!UH.is_allocated() || !UNa.is_allocated())
//             continue;
//         auto NUM = [&](int i, int j, int k, int m) -> double
//         {
//             return UH(i, j, k, 0) / M_H + UNa(i, j, k, 0) / M_Na;
//         };

//         auto &Bind_xi = fld_->field(fid_.fid_B.xi, iblk);
//         auto &Bind_eta = fld_->field(fid_.fid_B.eta, iblk);
//         auto &Bind_zeta = fld_->field(fid_.fid_B.zeta, iblk);
//         auto &Badd_xi = fld_->field(fid_.fid_Badd.xi, iblk);
//         auto &Badd_eta = fld_->field(fid_.fid_Badd.eta, iblk);
//         auto &Badd_zeta = fld_->field(fid_.fid_Badd.zeta, iblk);
//         auto Bxi = [&](int i, int j, int k, int m) -> double
//         {
//             return Bind_xi(i, j, k, m) + Badd_xi(i, j, k, m);
//         };
//         auto Beta = [&](int i, int j, int k, int m) -> double
//         {
//             return Bind_eta(i, j, k, m) + Badd_eta(i, j, k, m);
//         };
//         auto Bzeta = [&](int i, int j, int k, int m) -> double
//         {
//             return Bind_zeta(i, j, k, m) + Badd_zeta(i, j, k, m);
//         };

//         auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
//         auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
//         auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

//         auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
//         auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
//         auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);

//         // edge cache: 9 comps (row-major 3x3)
//         auto &pinvGT_xi = fld_->field(fid_.fid_pinvGT.xi, iblk);
//         auto &pinvGT_eta = fld_->field(fid_.fid_pinvGT.eta, iblk);
//         auto &pinvGT_zeta = fld_->field(fid_.fid_pinvGT.zeta, iblk);
//         auto &pinvAT_xi = fld_->field(fid_.fid_pinvAT.xi, iblk);
//         auto &pinvAT_eta = fld_->field(fid_.fid_pinvAT.eta, iblk);
//         auto &pinvAT_zeta = fld_->field(fid_.fid_pinvAT.zeta, iblk);

//         auto &dlx = fld_->field("dl_xi", iblk);
//         auto &dle = fld_->field("dl_eta", iblk);
//         auto &dlz = fld_->field("dl_zeta", iblk);

//         auto &x = grd_->grids(iblk).x;
//         auto &y = grd_->grids(iblk).y;
//         auto &z = grd_->grids(iblk).z;

//         // ---------- geometry / radius ----------
//         auto radius_xi = [&](int i, int j, int k) -> double
//         {
//             // 用 xi-edge 中点半径，比直接 radius(i,j,k) 更一致
//             double xx = 0.5 * (x(i + 1, j, k) + x(i, j, k));
//             double yy = 0.5 * (y(i + 1, j, k) + y(i, j, k));
//             double zz = 0.5 * (z(i + 1, j, k) + z(i, j, k));
//             return std::sqrt(xx * xx + yy * yy + zz * zz);
//         };
//         auto radius_eta = [&](int i, int j, int k) -> double
//         {
//             double xx = 0.5 * (x(i, j + 1, k) + x(i, j, k));
//             double yy = 0.5 * (y(i, j + 1, k) + y(i, j, k));
//             double zz = 0.5 * (z(i, j + 1, k) + z(i, j, k));
//             return std::sqrt(xx * xx + yy * yy + zz * zz);
//         };

//         auto radius_zeta = [&](int i, int j, int k) -> double
//         {
//             double xx = 0.5 * (x(i, j, k + 1) + x(i, j, k));
//             double yy = 0.5 * (y(i, j, k + 1) + y(i, j, k));
//             double zz = 0.5 * (z(i, j, k + 1) + z(i, j, k));
//             return std::sqrt(xx * xx + yy * yy + zz * zz);
//         };

//         // ---------- smooth 1/ne ----------
//         // Num : 1/M_H (mol/kg)
//         auto hall_alpha_from_ne = [&](double Num, double r) -> double
//         {
//             const double ne_true = std::max(Num, 0.0);

//             // 平滑 floor，避免硬 max
//             const double ne_eff = std::sqrt(ne_true * ne_true +
//                                             ne_hall_floor * ne_hall_floor);

//             // 平滑 cut，避免低密度区突然爆
//             const double s_ne = ne_true / (ne_true + ne_hall_cut);

//             if (r <= 1.01)
//                 return 0.0;
//             if (r >= 1.50)
//                 return hall_coeff * s_ne / ne_eff;

//             const double xi = (r - 1.01) / 0.49;         // xi in (0,1)
//             const double w = xi * xi * (3.0 - 2.0 * xi); // smoothstep

//             return hall_coeff * s_ne / ne_eff * w;

//             // if (r < 1.01)
//             //     return 0;
//             // else if (1.01 < r && r < 1.1)
//             //     return hall_coeff * s_ne / ne_eff * (r - 1.01) / 0.09;
//             // return hall_coeff * s_ne / ne_eff;
//         };

//         // ---------- radial dissipation ----------
//         auto hall_eta_by_r = [&](double r, double h2) -> double
//         {
//             // 老程序的 cb0(r) 形状
//             const double cb_inner = 10.0;
//             const double cb_outer = 5.0;
//             const double rr00 = 0.5 * (1.1 + 1.2); // 1.15
//             const double delta0 = 0.03;

//             const double cb0 =
//                 (cb_inner - cb_outer) * 0.5 * (1.0 - std::tanh((r - rr00) / delta0)) + cb_outer;

//             return C_eta * cb0 * h2;
//         };

//         // Geometry Tools
//         auto cache_is_valid = [&](FieldBlock &M9, int i, int j, int k) -> bool
//         {
//             double s = 0.0;
//             for (int m = 0; m < 9; ++m)
//                 s += std::abs(M9(i, j, k, m));
//             return s > 0.0;
//         };
//         auto hmin2 = [&](int i, int j, int k) -> double
//         {
//             constexpr double h_eps = 1e-12;
//             double hx = dlx.is_allocated() ? dlx(i, j, k, 0) : 1e100;
//             double he = dle.is_allocated() ? dle(i, j, k, 0) : 1e100;
//             double hz = dlz.is_allocated() ? dlz(i, j, k, 0) : 1e100;
//             double h = std::min(hx, std::min(he, hz));
//             if (h <= h_eps)
//                 return 0.0;
//             return h * h;
//         };

//         auto matvec3 = [&](FieldBlock &M9, int i, int j, int k,
//                            double c0, double c1, double c2) -> Vec3
//         {
//             Vec3 v;
//             v.vec[0] = M9(i, j, k, 0) * c0 + M9(i, j, k, 1) * c1 + M9(i, j, k, 2) * c2;
//             v.vec[1] = M9(i, j, k, 3) * c0 + M9(i, j, k, 4) * c1 + M9(i, j, k, 5) * c2;
//             v.vec[2] = M9(i, j, k, 6) * c0 + M9(i, j, k, 7) * c1 + M9(i, j, k, 8) * c2;
//             return v;
//         };

//         // ============================================================
//         // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
//         // ============================================================
//         {
//             Int3 lo = Ehall_xi.inner_lo();
//             Int3 hi = Ehall_xi.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         if (!cache_is_valid(pinvGT_xi, i, j, k) || !cache_is_valid(pinvAT_xi, i, j, k))
//                         {
//                             Ehall_xi(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double h2 = hmin2(i, j, k);
//                         if (h2 <= 0.0)
//                         {
//                             Ehall_xi(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double r = radius_xi(i, j, k);

//                         // rho at xi-edge: average 4 surrounding cells
//                         double Num = 0.25 * (NUM(i, j, k, 0) +
//                                              NUM(i, j - 1, k, 0) +
//                                              NUM(i, j, k - 1, 0) +
//                                              NUM(i, j - 1, k - 1, 0));

//                         // Hall 系数：hall/(e*n_e)
//                         const double alpha = hall_alpha_from_ne(Num, r);

//                         // Phi (2-form) co-located at xi-edge center
//                         double Phi_xi = 0.0;
//                         for (int di : {0, 1})
//                             for (int dj : {0, -1})
//                                 for (int dk : {0, -1})
//                                     Phi_xi += Bxi(i + di, j + dj, k + dk, 0);
//                         Phi_xi *= 0.125;

//                         double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i, j, k - 1, 0));
//                         double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i, j - 1, k, 0));

//                         // j (1-form) co-located at same xi-edge center
//                         double j_xi = Jxi(i, j, k, 0);
//                         double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i + 1, j, k, 0) + Jeta(i, j - 1, k, 0) + Jeta(i + 1, j - 1, k, 0));
//                         double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i + 1, j, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i + 1, j, k - 1, 0));

//                         // map to physical vectors
//                         Vec3 Jvec = matvec3(pinvGT_xi, i, j, k, j_xi, j_eta, j_zeta);
//                         Vec3 Bvec = matvec3(pinvAT_xi, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

//                         // Hall part
//                         Vec3 Evec = (Jvec ^ Bvec);
//                         Evec *= alpha;

//                         // radial dissipation part
//                         const double eta_h = hall_eta_by_r(r, h2);
//                         Evec.vec[0] += eta_h * Jvec.vec[0];
//                         Evec.vec[1] += eta_h * Jvec.vec[1];
//                         Evec.vec[2] += eta_h * Jvec.vec[2];

//                         Vec3 dr;
//                         dr.vec[0] = x(i + 1, j, k) - x(i, j, k);
//                         dr.vec[1] = y(i + 1, j, k) - y(i, j, k);
//                         dr.vec[2] = z(i + 1, j, k) - z(i, j, k);

//                         Ehall_xi(i, j, k, 0) = Evec * dr;
//                     }
//         }

//         // ============================================================
//         // 2) EdgeEta : Ehall_eta(i,j,k) = (alpha * (J x B) + eta_h J) · dr_eta
//         // ============================================================
//         {
//             Int3 lo = Ehall_eta.inner_lo();
//             Int3 hi = Ehall_eta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         if (!cache_is_valid(pinvGT_eta, i, j, k) || !cache_is_valid(pinvAT_eta, i, j, k))
//                         {
//                             Ehall_eta(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double h2 = hmin2(i, j, k);
//                         if (h2 <= 0.0)
//                         {
//                             Ehall_eta(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double r = radius_eta(i, j, k);

//                         // rho at eta-edge: average 4 surrounding cells (i- and k- directions)
//                         double Num = 0.25 * (NUM(i, j, k, 0) +
//                                              NUM(i - 1, j, k, 0) +
//                                              NUM(i, j, k - 1, 0) +
//                                              NUM(i - 1, j, k - 1, 0));

//                         // Hall 系数：hall/(e*n_e)
//                         const double alpha = hall_alpha_from_ne(Num, r);

//                         // Phi (2-form) co-located at eta-edge center
//                         double Phi_eta = 0.0;
//                         for (int di : {0, -1})
//                             for (int dj : {0, 1})
//                                 for (int dk : {0, -1})
//                                     Phi_eta += Beta(i + di, j + dj, k + dk, 0);
//                         Phi_eta *= 0.125;

//                         double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j, k - 1, 0));
//                         double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i - 1, j, k, 0));

//                         // j (1-form) co-located at same eta-edge center
//                         double j_eta = Jeta(i, j, k, 0);
//                         double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j + 1, k, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j + 1, k, 0));
//                         double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i, j + 1, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i, j + 1, k - 1, 0));

//                         // map to physical vectors
//                         Vec3 Jvec = matvec3(pinvGT_eta, i, j, k, j_xi, j_eta, j_zeta);
//                         Vec3 Bvec = matvec3(pinvAT_eta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

//                         // Hall part
//                         Vec3 Evec = (Jvec ^ Bvec);
//                         Evec *= alpha;

//                         // radial dissipation part
//                         const double eta_h = hall_eta_by_r(r, h2);
//                         Evec.vec[0] += eta_h * Jvec.vec[0];
//                         Evec.vec[1] += eta_h * Jvec.vec[1];
//                         Evec.vec[2] += eta_h * Jvec.vec[2];

//                         Vec3 dr;
//                         dr.vec[0] = x(i, j + 1, k) - x(i, j, k);
//                         dr.vec[1] = y(i, j + 1, k) - y(i, j, k);
//                         dr.vec[2] = z(i, j + 1, k) - z(i, j, k);

//                         Ehall_eta(i, j, k, 0) = Evec * dr;
//                     }
//         }

//         // ============================================================
//         // 3) EdgeZeta : Ehall_zeta(i,j,k) = (alpha * (J x B) + eta_h J) · dr_zeta
//         // ============================================================
//         {
//             Int3 lo = Ehall_zeta.inner_lo();
//             Int3 hi = Ehall_zeta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         if (!cache_is_valid(pinvGT_zeta, i, j, k) || !cache_is_valid(pinvAT_zeta, i, j, k))
//                         {
//                             Ehall_zeta(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double h2 = hmin2(i, j, k);
//                         if (h2 <= 0.0)
//                         {
//                             Ehall_zeta(i, j, k, 0) = 0.0;
//                             continue;
//                         }

//                         const double r = radius_zeta(i, j, k);

//                         // rho at zeta-edge: average 4 surrounding cells (i- and j- directions)
//                         double Num = 0.25 * (NUM(i, j, k, 0) +
//                                              NUM(i - 1, j, k, 0) +
//                                              NUM(i, j - 1, k, 0) +
//                                              NUM(i - 1, j - 1, k, 0));

//                         // Hall 系数：hall/(e*n_e)
//                         const double alpha = hall_alpha_from_ne(Num, r);

//                         // Phi (2-form) co-located at zeta-edge center
//                         double Phi_zeta = 0.0;
//                         for (int di : {0, -1})
//                             for (int dj : {0, -1})
//                                 for (int dk : {0, 1})
//                                     Phi_zeta += Bzeta(i + di, j + dj, k + dk, 0);
//                         Phi_zeta *= 0.125;

//                         double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j - 1, k, 0));
//                         double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i - 1, j, k, 0));

//                         // j (1-form) co-located at same zeta-edge center
//                         double j_zeta = Jzeta(i, j, k, 0);
//                         double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j, k + 1, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j, k + 1, 0));
//                         double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i, j, k + 1, 0) + Jeta(i, j - 1, k, 0) + Jeta(i, j - 1, k + 1, 0));

//                         // map to physical vectors
//                         Vec3 Jvec = matvec3(pinvGT_zeta, i, j, k, j_xi, j_eta, j_zeta);
//                         Vec3 Bvec = matvec3(pinvAT_zeta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

//                         // Hall part
//                         Vec3 Evec = (Jvec ^ Bvec);
//                         Evec *= alpha;

//                         // radial dissipation part
//                         const double eta_h = hall_eta_by_r(r, h2);
//                         Evec.vec[0] += eta_h * Jvec.vec[0];
//                         Evec.vec[1] += eta_h * Jvec.vec[1];
//                         Evec.vec[2] += eta_h * Jvec.vec[2];

//                         Vec3 dr;
//                         dr.vec[0] = x(i, j, k + 1) - x(i, j, k);
//                         dr.vec[1] = y(i, j, k + 1) - y(i, j, k);
//                         dr.vec[2] = z(i, j, k + 1) - z(i, j, k);

//                         Ehall_zeta(i, j, k, 0) = Evec * dr;
//                     }
//         }
//     }

//     // Add boundary conditions to E_hall
//     mercury_bound_.Sync("Ehall");

//     // Add E_hall to E(Edge)
//     for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
//     {
//         auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
//         auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
//         auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);
//         if (!Ehall_xi.is_allocated() || !Ehall_eta.is_allocated() || !Ehall_zeta.is_allocated())
//             continue;

//         auto &E_xi = fld_->field(fid_.fid_E.xi, iblk);
//         auto &E_eta = fld_->field(fid_.fid_E.eta, iblk);
//         auto &E_zeta = fld_->field(fid_.fid_E.zeta, iblk);

//         if (!E_xi.is_allocated() || !E_eta.is_allocated() || !E_zeta.is_allocated())
//             continue;

//         // ============================================================
//         // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
//         // ============================================================
//         {
//             Int3 lo = E_xi.inner_lo();
//             Int3 hi = E_xi.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         E_xi(i, j, k, 0) += Ehall_xi(i, j, k, 0);
//         }

//         // ============================================================
//         // 2) EdgeEt
//         // ============================================================
//         {
//             Int3 lo = E_eta.inner_lo();
//             Int3 hi = E_eta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         E_eta(i, j, k, 0) += Ehall_eta(i, j, k, 0);
//         }

//         // ============================================================
//         // 3) EdgeZe
//         // ============================================================
//         {
//             Int3 lo = E_zeta.inner_lo();
//             Int3 hi = E_zeta.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                         E_zeta(i, j, k, 0) += Ehall_zeta(i, j, k, 0);
//         }
//     }
// }
