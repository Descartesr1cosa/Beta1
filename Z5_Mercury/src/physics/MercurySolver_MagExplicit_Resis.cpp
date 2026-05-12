#include "MercurySolver.h"

void MercurySolver::AddResistiveEdgeEMF_()
{
    if (!resist_control.is_Mercury_resistance)
        return;

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

        auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
        auto &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
        auto &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

        auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

        if (!Jcell.is_allocated())
            continue;

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        auto get_edge_cells = [](int edge_axis, int i, int j, int k)
        {
            std::array<Int3, 4> c;

            if (edge_axis == 0)
            {
                c[0] = {i, j, k};
                c[1] = {i, j - 1, k};
                c[2] = {i, j, k - 1};
                c[3] = {i, j - 1, k - 1};
            }
            else if (edge_axis == 1)
            {
                c[0] = {i, j, k};
                c[1] = {i - 1, j, k};
                c[2] = {i, j, k - 1};
                c[3] = {i - 1, j, k - 1};
            }
            else
            {
                c[0] = {i, j, k};
                c[1] = {i - 1, j, k};
                c[2] = {i, j - 1, k};
                c[3] = {i - 1, j - 1, k};
            }

            return c;
        };

        auto add_one_edge = [&](FieldBlock &E,
                                FieldBlock &dr,
                                int edge_axis,
                                int di,
                                int dj,
                                int dk)
        {
            if (!E.is_allocated() || !dr.is_allocated())
                return;

            Int3 lo = E.inner_lo();
            Int3 hi = E.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i + di, j + dj, k + dk));
                        const double ym = 0.5 * (y(i, j, k) + y(i + di, j + dj, k + dk));
                        const double zm = 0.5 * (z(i, j, k) + z(i + di, j + dj, k + dk));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        const auto cells = get_edge_cells(edge_axis, i, j, k);
                        double Jx = 0.0;
                        double Jy = 0.0;
                        double Jz = 0.0;

                        for (int q = 0; q < 4; ++q)
                        {
                            const int ic = cells[q].i;
                            const int jc = cells[q].j;
                            const int kc = cells[q].k;
                            Jx += Jcell(ic, jc, kc, 0);
                            Jy += Jcell(ic, jc, kc, 1);
                            Jz += Jcell(ic, jc, kc, 2);
                        }

                        Jx *= 0.25;
                        Jy *= 0.25;
                        Jz *= 0.25;

                        const double Jedge =
                            Jx * dr(i, j, k, 0) +
                            Jy * dr(i, j, k, 1) +
                            Jz * dr(i, j, k, 2);

                        E(i, j, k, 0) += (inver_Rem * yita0) * Jedge;
                    }
        };

        // --- EdgeXi: use the edge midpoint of node (i,j,k) -> (i+1,j,k) ---
        {
            add_one_edge(Exi, dr_xi, 0, 1, 0, 0);
        }

        // --- EdgeEt: node (i,j,k) -> (i,j+1,k) ---
        {
            add_one_edge(Eeta, dr_eta, 1, 0, 1, 0);
        }

        // --- EdgeZe: node (i,j,k) -> (i,j,k+1) ---
        {
            add_one_edge(Eze, dr_zeta, 2, 0, 0, 1);
        }
    }
}

void MercurySolver::AddPoleResistiveEdgeEMF_FromJcell_()
{
    const auto &ctrl = resist_control;

    if (!ctrl.enabled)
        return;

    if (ctrl.eta_max <= 0.0)
        return;

    auto smoothstep01 = [](double x) -> double
    {
        x = std::max(0.0, std::min(1.0, x));
        return x * x * (3.0 - 2.0 * x);
    };

    auto smoothstep = [&](double a, double b, double x) -> double
    {
        if (b <= a)
            return x >= a ? 1.0 : 0.0;

        return smoothstep01((x - a) / (b - a));
    };

    auto get_edge_cells = [](int edge_axis, int i, int j, int k)
    {
        std::array<Int3, 4> c;

        if (edge_axis == 0) // Exi
        {
            c[0] = {i, j, k};
            c[1] = {i, j - 1, k};
            c[2] = {i, j, k - 1};
            c[3] = {i, j - 1, k - 1};
        }
        else if (edge_axis == 1) // Eeta
        {
            c[0] = {i, j, k};
            c[1] = {i - 1, j, k};
            c[2] = {i, j, k - 1};
            c[3] = {i - 1, j, k - 1};
        }
        else // Ezeta
        {
            c[0] = {i, j, k};
            c[1] = {i - 1, j, k};
            c[2] = {i, j - 1, k};
            c[3] = {i - 1, j - 1, k};
        }

        return c;
    };

    auto add_one_edge = [&](FieldBlock &E,
                            FieldBlock &dr,
                            FieldBlock &Jcell,
                            int edge_axis)
    {
        if (!E.is_allocated() || !dr.is_allocated() || !Jcell.is_allocated())
            return;

        const Int3 elo = E.inner_lo();
        const Int3 ehi = E.inner_hi();

        for (int i = elo.i; i < ehi.i; ++i)
            for (int j = elo.j; j < ehi.j; ++j)
                for (int k = elo.k; k < ehi.k; ++k)
                {
                    const auto cells = get_edge_cells(edge_axis, i, j, k);

                    double Jx_sum = 0.0;
                    double Jy_sum = 0.0;
                    double Jz_sum = 0.0;

                    double Jmag[4];

                    double Jmax = 0.0;
                    double Jmean = 0.0;

                    for (int q = 0; q < 4; ++q)
                    {
                        const int ic = cells[q].i;
                        const int jc = cells[q].j;
                        const int kc = cells[q].k;

                        const double Jx = Jcell(ic, jc, kc, 0);
                        const double Jy = Jcell(ic, jc, kc, 1);
                        const double Jz = Jcell(ic, jc, kc, 2);

                        Jx_sum += Jx;
                        Jy_sum += Jy;
                        Jz_sum += Jz;

                        Jmag[q] = std::sqrt(Jx * Jx + Jy * Jy + Jz * Jz);

                        Jmax = std::max(Jmax, Jmag[q]);
                        Jmean += Jmag[q];
                    }

                    Jmean *= 0.25;

                    double var = 0.0;
                    for (int q = 0; q < 4; ++q)
                    {
                        const double d = Jmag[q] - Jmean;
                        var += d * d;
                    }
                    var *= 0.25;

                    const double cv =
                        std::sqrt(var) / (Jmean + ctrl.eps);

                    const double Jx_avg = 0.25 * Jx_sum;
                    const double Jy_avg = 0.25 * Jy_sum;
                    const double Jz_avg = 0.25 * Jz_sum;

                    const double Javg_mag =
                        std::sqrt(Jx_avg * Jx_avg +
                                  Jy_avg * Jy_avg +
                                  Jz_avg * Jz_avg);

                    double cancel =
                        1.0 - Javg_mag / (Jmean + ctrl.eps);

                    cancel = std::max(0.0, std::min(1.0, cancel));

                    const double osc = std::max(cv, cancel);

                    const double S_J =
                        smoothstep(ctrl.J_on, ctrl.J_full, Jmax);

                    const double S_osc =
                        smoothstep(ctrl.osc_on, ctrl.osc_full, osc);

                    const double S_bad = S_J * S_osc;

                    const double S_force =
                        (Jmax >= ctrl.J_max_force) ? 1.0 : 0.0;

                    const double S =
                        std::max(S_bad, S_force);

                    if (S <= 0.0)
                        continue;

                    const double eta_eff = ctrl.eta_max * S;

                    const double Jdotdr =
                        Jx_avg * dr(i, j, k, 0) +
                        Jy_avg * dr(i, j, k, 1) +
                        Jz_avg * dr(i, j, k, 2);

                    E(i, j, k, 0) += eta_eff * Jdotdr;
                }
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Exi = fld_->field(fid_.fid_E.xi, ib);
        FieldBlock &Eeta = fld_->field(fid_.fid_E.eta, ib);
        FieldBlock &Eze = fld_->field(fid_.fid_E.zeta, ib);

        FieldBlock &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
        FieldBlock &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
        FieldBlock &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

        FieldBlock &Jcell = fld_->field(fid_.fid_Jcell, ib);

        add_one_edge(Exi, dr_xi, Jcell, 0);
        add_one_edge(Eeta, dr_eta, Jcell, 1);
        add_one_edge(Eze, dr_zeta, Jcell, 2);
    }
}

// void MercurySolver::AddPoleResistiveEdgeEMF_FromJcell_()
// {
//     constexpr int shift_layers = 4;
//     constexpr double radius_limit = 1.2;
//     constexpr double filter_floor = 0.0;

//     auto loc_delta = [](int edge_axis) -> Int3
//     {
//         if (edge_axis == 0)
//             return {1, 0, 0};
//         if (edge_axis == 1)
//             return {0, 1, 0};
//         return {0, 0, 1};
//     };

//     auto get_comp = [](const Int3 &a, int ax) -> int
//     {
//         if (ax == 0)
//             return a.i;
//         if (ax == 1)
//             return a.j;
//         return a.k;
//     };

//     auto set_comp = [](Int3 &a, int ax, int v)
//     {
//         if (ax == 0)
//             a.i = v;
//         else if (ax == 1)
//             a.j = v;
//         else
//             a.k = v;
//     };

//     auto add_one_edge = [&](FieldBlock &E,
//                             FieldBlock &dr,
//                             int edge_axis,
//                             const TOPO::PhysicalPatch &p,
//                             double3D &x,
//                             double3D &y,
//                             double3D &z,
//                             FieldBlock &Jcell)
//     {
//         if (!E.is_allocated() || !dr.is_allocated() || !Jcell.is_allocated())
//             return;

//         const int pole_dir = std::abs(p.direction);
//         if (pole_dir != 1 && pole_dir != 2)
//             return;

//         const int norm_axis = pole_dir - 1;
//         const bool high_side = (p.direction > 0);
//         const Int3 delta = loc_delta(edge_axis);

//         const Int3 elo = E.inner_lo();
//         const Int3 ehi = E.inner_hi();

//         Int3 base_lo = elo;
//         Int3 base_hi = ehi;

//         for (int ax = 0; ax < 3; ++ax)
//         {
//             if (ax == norm_axis)
//                 continue;

//             const int d = get_comp(delta, ax);
//             const int patch_lo = get_comp(p.this_box_node.lo, ax);
//             const int patch_hi = get_comp(p.this_box_node.hi, ax) - d;

//             set_comp(base_lo, ax, std::max(get_comp(elo, ax), patch_lo));
//             set_comp(base_hi, ax, std::min(get_comp(ehi, ax), patch_hi));
//         }

//         const int nlo = get_comp(elo, norm_axis);
//         const int nhi = get_comp(ehi, norm_axis);
//         const int max_layers = std::min(shift_layers, nhi - nlo);

//         if (max_layers <= 0)
//             return;

//         auto edge_mid_radius = [&](int i, int j, int k) -> double
//         {
//             const int ip = i + (edge_axis == 0 ? 1 : 0);
//             const int jp = j + (edge_axis == 1 ? 1 : 0);
//             const int kp = k + (edge_axis == 2 ? 1 : 0);

//             const double xm = 0.5 * (x(i, j, k) + x(ip, jp, kp));
//             const double ym = 0.5 * (y(i, j, k) + y(ip, jp, kp));
//             const double zm = 0.5 * (z(i, j, k) + z(ip, jp, kp));

//             return std::sqrt(xm * xm + ym * ym + zm * zm);
//         };

//         auto jcell_to_edge = [&](int i, int j, int k) -> double
//         {
//             double Jx = 0.0;
//             double Jy = 0.0;
//             double Jz = 0.0;

//             if (edge_axis == 0)
//             {
//                 Jx = 0.25 * (Jcell(i, j - 1, k - 1, 0) + Jcell(i, j - 1, k, 0) +
//                              Jcell(i, j, k - 1, 0) + Jcell(i, j, k, 0));
//                 Jy = 0.25 * (Jcell(i, j - 1, k - 1, 1) + Jcell(i, j - 1, k, 1) +
//                              Jcell(i, j, k - 1, 1) + Jcell(i, j, k, 1));
//                 Jz = 0.25 * (Jcell(i, j - 1, k - 1, 2) + Jcell(i, j - 1, k, 2) +
//                              Jcell(i, j, k - 1, 2) + Jcell(i, j, k, 2));
//             }
//             else if (edge_axis == 1)
//             {
//                 Jx = 0.25 * (Jcell(i - 1, j, k - 1, 0) + Jcell(i - 1, j, k, 0) +
//                              Jcell(i, j, k - 1, 0) + Jcell(i, j, k, 0));
//                 Jy = 0.25 * (Jcell(i - 1, j, k - 1, 1) + Jcell(i - 1, j, k, 1) +
//                              Jcell(i, j, k - 1, 1) + Jcell(i, j, k, 1));
//                 Jz = 0.25 * (Jcell(i - 1, j, k - 1, 2) + Jcell(i - 1, j, k, 2) +
//                              Jcell(i, j, k - 1, 2) + Jcell(i, j, k, 2));
//             }
//             else
//             {
//                 Jx = 0.25 * (Jcell(i - 1, j - 1, k, 0) + Jcell(i - 1, j, k, 0) +
//                              Jcell(i, j - 1, k, 0) + Jcell(i, j, k, 0));
//                 Jy = 0.25 * (Jcell(i - 1, j - 1, k, 1) + Jcell(i - 1, j, k, 1) +
//                              Jcell(i, j - 1, k, 1) + Jcell(i, j, k, 1));
//                 Jz = 0.25 * (Jcell(i - 1, j - 1, k, 2) + Jcell(i - 1, j, k, 2) +
//                              Jcell(i, j - 1, k, 2) + Jcell(i, j, k, 2));
//             }

//             return Jx * dr(i, j, k, 0) + Jy * dr(i, j, k, 1) + Jz * dr(i, j, k, 2);
//         };

//         for (int layer = 0; layer < max_layers; ++layer)
//         {
//             Int3 lo = base_lo;
//             Int3 hi = base_hi;

//             const int n = high_side ? (nhi - 1 - layer) : (nlo + layer);
//             set_comp(lo, norm_axis, n);
//             set_comp(hi, norm_axis, n + 1);

//             if (!(lo.i < hi.i && lo.j < hi.j && lo.k < hi.k))
//                 continue;

//             double weight = 1.0 - static_cast<double>(layer) / static_cast<double>(max_layers);
//             weight = std::max(filter_floor, std::min(1.0, weight));

//             if (weight <= 0.0)
//                 continue;

//             const double eta = 0.05 * weight;

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         if (edge_mid_radius(i, j, k) > radius_limit)
//                             continue;

//                         E(i, j, k, 0) += eta * jcell_to_edge(i, j, k);
//                     }
//         }
//     };

//     const int nb = fld_->num_blocks();
//     for (int ib = 0; ib < nb; ++ib)
//     {
//         FieldBlock &Exi = fld_->field(fid_.fid_E.xi, ib);
//         FieldBlock &Eeta = fld_->field(fid_.fid_E.eta, ib);
//         FieldBlock &Eze = fld_->field(fid_.fid_E.zeta, ib);

//         FieldBlock &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
//         FieldBlock &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
//         FieldBlock &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

//         FieldBlock &Jcell = fld_->field(fid_.fid_Jcell, ib);

//         auto &x = grd_->grids(ib).x;
//         auto &y = grd_->grids(ib).y;
//         auto &z = grd_->grids(ib).z;

//         for (const auto &p : topo_->physical_patches)
//         {
//             if (p.this_block != ib)
//                 continue;
//             if (p.bc_name != "Pole")
//                 continue;

//             add_one_edge(Exi, dr_xi, 0, p, x, y, z, Jcell);
//             add_one_edge(Eeta, dr_eta, 1, p, x, y, z, Jcell);
//             add_one_edge(Eze, dr_zeta, 2, p, x, y, z, Jcell);
//         }
//     }
// }

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
        auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);

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
        auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);

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
