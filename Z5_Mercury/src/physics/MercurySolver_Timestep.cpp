
#include <cmath>
#include <algorithm>
#include <iostream>

#include "MercurySolver.h"

void MercurySolver::Compute_Timestep()
{
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    double dt_local = 1e100;

    auto norm3 = [](double x, double y, double z)
    { return std::sqrt(x * x + y * y + z * z); };

    double dt_mhd_min_local = 1e100;
    double dt_hall_min_local = 1e100;

    auto scan_one = [&](int fidU)
    {
        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            FieldBlock &U = fld_->field(fidU, ib);
            FieldBlock &Bcell = fld_->field(fid_.fid_Bcell, ib);
            FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
            FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
            FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
            FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

            if (!U.is_allocated())
                continue;

            Int3 lo = Jac.inner_lo();
            Int3 hi = Jac.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double V = std::abs(Jac(i, j, k, 0));
                        if (V <= 0.0)
                            continue;

                        double rho = std::max(U(i, j, k, 0), rho_floor);
                        double ux = U(i, j, k, 1) / rho;
                        double uy = U(i, j, k, 2) / rho;
                        double uz = U(i, j, k, 3) / rho;
                        const double E = U(i, j, k, 4);
                        const double ke = 0.5 * rho * (ux * ux + uy * uy + uz * uz);

                        // 压力：p = (gamma-1) * (E - ke)
                        double eint = E - ke;
                        double p = (gamma_ - 1.0) * eint;

                        double cs2 = gamma_ * p / rho;

                        double Bx = Bcell(i, j, k, 0);
                        double By = Bcell(i, j, k, 1);
                        double Bz = Bcell(i, j, k, 2);
                        double B2 = Bx * Bx + By * By + Bz * Bz;

                        auto fast_cf = [&](double vA2, double vAn2) -> double
                        {
                            double term = cs2 + vA2;
                            double disc = term * term - 4.0 * cs2 * vAn2;
                            if (disc < 0.0)
                                disc = 0.0; // 数值保护
                            return std::sqrt(0.5 * (term + std::sqrt(disc)));
                        };

                        auto face_term = [&](double ax, double ay, double az, bool outward_plus)
                        {
                            double A = norm3(ax, ay, az);
                            if (A <= 0.0)
                                return 0.0;
                            double nx = ax / A, ny = ay / A, nz = az / A;
                            if (!outward_plus)
                            {
                                nx = -nx;
                                ny = -ny;
                                nz = -nz;
                            } // outward to minus side
                            double un = ux * nx + uy * ny + uz * nz;
                            double Bn = Bx * nx + By * ny + Bz * nz;
                            double vA2 = inver_MA2 * B2 / rho;
                            double vAn2 = inver_MA2 * (Bn * Bn) / rho;
                            double cf = fast_cf(vA2, vAn2);
                            return (std::abs(un) + cf) * A;
                        };

                        double denom = 0.0;

                        // xi+ at i  : outward = Axi(i)
                        denom += face_term(Axi(i, j, k, 0), Axi(i, j, k, 1), Axi(i, j, k, 2), false);
                        // xi- at i+1: outward = -Axi(i+1)
                        denom += face_term(Axi(i + 1, j, k, 0), Axi(i + 1, j, k, 1), Axi(i + 1, j, k, 2), true);

                        // eta+ at j
                        denom += face_term(Aet(i, j, k, 0), Aet(i, j, k, 1), Aet(i, j, k, 2), false);
                        // eta- at j+1
                        denom += face_term(Aet(i, j + 1, k, 0), Aet(i, j + 1, k, 1), Aet(i, j + 1, k, 2), true);

                        // zeta+ at k
                        denom += face_term(Aze(i, j, k, 0), Aze(i, j, k, 1), Aze(i, j, k, 2), false);
                        // zeta- at k+1
                        denom += face_term(Aze(i, j, k + 1, 0), Aze(i, j, k + 1, 1), Aze(i, j, k + 1, 2), true);

                        if (denom > 0.0)
                        {
                            double dtc = CFL * V / denom;
                            dt_mhd_min_local = std::min(dt_mhd_min_local, dtc);
                            dt_local = std::min(dt_local, dtc);
                        }
                    }
        }
    };

    scan_one(fid_.fid_U_H);
    scan_one(fid_.fid_U_Na);

    // ---------------------------
    // (2) NEW: Hall 显式稳定步长限制（whistler: dt ~ h^2 / (alpha |B|)）
    // ---------------------------
    double dt_hall_min_l = 1e100;
    double h_at_min_l = 0.0, ne_at_min_l = 0.0, B_at_min_l = 0.0, alpha_at_min_l = 0.0;
    int ib_at_min_l = -1, i_at_min_l = -1, j_at_min_l = -1, k_at_min_l = -1, dir_at_min_l = -1;

    // ===========================
    // Hall explicit dt constraint (whistler) with axis-degenerate exclusion for EdgeZeta
    // ===========================
    if (std::abs(hall_coef) > 0.0)
    {
        const double CFL_HALL = 0.05;          // 先用更保守的 0.05；你之前 0.2 太松
        const double ne_floor = ne_hall_floor; // 1e-30;
        const double h_eps_abs = 1e-12;
        const double B_floor = 1e-30;

        // 参考量：rho_ref & particle masses（同你之前）
        const auto &C = par_->GetDou_List("constant").data;
        const auto &R = par_->GetDou_List("REF").data;
        const double NA = C.at("NA");
        const double n_ref = R.at("n");
        const double Mref = R.at("Molecular_mass");
        const double rho_ref = (Mref / NA) * n_ref;

        const double mH = par_->GetDou("mole_mass1");
        const double mNa = par_->GetDou("mole_mass2");

        auto ne_cell = [&](FieldBlock &UH, FieldBlock &UNa, int i, int j, int k) -> double
        {
            const double rhoH_nd = std::max(0.0, UH(i, j, k, 0));
            const double rhoNa_nd = std::max(0.0, UNa(i, j, k, 0));
            // const double nH = (rhoH_nd * rho_ref) / mH;
            // const double nNa = (rhoNa_nd * rho_ref) / mNa;
            const double nH = (rhoH_nd) / mH;
            const double nNa = (rhoNa_nd) / mNa;
            return nH + nNa;
        };

        auto Babs_cell = [&](FieldBlock &Bcell, int i, int j, int k) -> double
        {
            const double bx = Bcell(i, j, k, 0);
            const double by = Bcell(i, j, k, 1);
            const double bz = Bcell(i, j, k, 2);
            return std::sqrt(bx * bx + by * by + bz * bz);
        };

        double hmin_pol_l = 1e100;
        double alpha_max_l = 0.0;
        double Bmax_l = 0.0;

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &UH = fld_->field(fid_.fid_U_H, ib);
            auto &UNa = fld_->field(fid_.fid_U_Na, ib);
            auto &Bcel = fld_->field(fid_.fid_Bcell, ib);

            auto &dlx = fld_->field("dl_xi", ib);
            auto &dle = fld_->field("dl_eta", ib);

            if (!UH.is_allocated() || !UNa.is_allocated() || !Bcel.is_allocated())
                continue;
            if (!dlx.is_allocated() || !dle.is_allocated())
                continue;

            // (a) scan poloidal min edge length
            {
                Int3 lo = dlx.inner_lo(), hi = dlx.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            double h = dlx(i, j, k, 0);
                            if (h > h_eps_abs)
                                hmin_pol_l = std::min(hmin_pol_l, h);
                        }
            }
            {
                Int3 lo = dle.inner_lo(), hi = dle.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            double h = dle(i, j, k, 0);
                            if (h > h_eps_abs)
                                hmin_pol_l = std::min(hmin_pol_l, h);
                        }
            }

            // (b) scan alpha_max and Bmax on cells
            {
                Int3 lo = UH.inner_lo(), hi = UH.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double ne = ne_cell(UH, UNa, i, j, k);
                            const double alpha = std::abs(hall_coef) / (ne + ne_floor);
                            alpha_max_l = std::max(alpha_max_l, alpha);

                            const double Babs = Babs_cell(Bcel, i, j, k);
                            Bmax_l = std::max(Bmax_l, Babs);
                        }
            }
        }

        double hmin_pol_g = hmin_pol_l;
        double alpha_max_g = alpha_max_l;
        double Bmax_g = Bmax_l;

        PARALLEL::mpi_min(&hmin_pol_l, &hmin_pol_g, 1);
        PARALLEL::mpi_max(&alpha_max_l, &alpha_max_g, 1);
        PARALLEL::mpi_max(&Bmax_l, &Bmax_g, 1);

        // if something went wrong
        if (!(hmin_pol_g < 1e99))
            hmin_pol_g = 1.0;

        if (Bmax_g > B_floor && alpha_max_g > 0.0)
        {
            const double dt_hall_global = CFL_HALL * (hmin_pol_g * hmin_pol_g) / (alpha_max_g * Bmax_g + 1e-300);
            dt_hall_min_local = std::min(dt_hall_min_local, dt_hall_global);
            // dt_local = std::min(dt_local, dt_hall_global);
        }

        if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
        {
            std::printf("[HallDtGlobal] step=%d  hmin_pol=%.3e  alpha_max=%.3e  Bmax=%.3e  dt_hall=%.3e\n",
                        io_.Run().step, hmin_pol_g, alpha_max_g, Bmax_g, dt_hall_min_local);
            std::fflush(stdout);
        }
    }

    // if (std::abs(hall_coef) > 0.0) // 只有 Hall 打开时才生效
    // {
    //     // 你可以把 0.2 改成参数；先用保守值稳住
    //     const double CFL_HALL = 0.2;

    //     const double ne_floor = 1e-30; // number density floor (m^-3)
    //     const double B_floor = 1e-30;
    //     const double h_floor = 1e-30;

    //     // 物理常数/参考量（用于从无量纲 rho -> ne）
    //     const auto &C = par_->GetDou_List("constant").data;
    //     const auto &R = par_->GetDou_List("REF").data;
    //     const double NA = C.at("NA");
    //     const double n_ref = R.at("n");             // 1/m^3
    //     const double Mref = R.at("Molecular_mass"); // kg/mol
    //     const double rho_ref = (Mref / NA) * n_ref; // kg/m^3

    //     const double mH = par_->GetDou("mole_mass1") / NA;  // kg/particle
    //     const double mNa = par_->GetDou("mole_mass2") / NA; // kg/particle

    //     auto NUM = [&](FieldBlock &UH, FieldBlock &UNa, int i, int j, int k) -> double
    //     {
    //         // UH/UNa 的 rho_nd * rho_ref => kg/m^3，再除粒子质量 => 1/m^3
    //         const double rhoH_nd = std::max(0.0, UH(i, j, k, 0));
    //         const double rhoNa_nd = std::max(0.0, UNa(i, j, k, 0));
    //         const double nH = (rhoH_nd * rho_ref) / mH;
    //         const double nNa = (rhoNa_nd * rho_ref) / mNa;
    //         return nH + nNa;
    //     };

    //     auto Babs_cell = [&](FieldBlock &Bcell, int i, int j, int k) -> double
    //     {
    //         const double bx = Bcell(i, j, k, 0);
    //         const double by = Bcell(i, j, k, 1);
    //         const double bz = Bcell(i, j, k, 2);
    //         return std::sqrt(bx * bx + by * by + bz * bz);
    //     };

    //     auto update_dt_hall = [&](double h, double ne, double Babs)
    //     {
    //         if (h <= h_floor)
    //             return; // 退化边：line integral ~ 0，这里不让它把 dt 压成 0
    //         if (Babs < B_floor)
    //             return;

    //         const double alpha = std::abs(hall_coef) / (ne + ne_floor);
    //         const double dth = CFL_HALL * (h * h) / (alpha * Babs + 1e-300);
    //         dt_local = std::min(dt_local, dth);
    //     };

    //     auto update_dt_hall_dbg = [&](double h, double ne, double Babs, int ib, int i, int j, int k, int dir)
    //     {
    //         if (h <= 0.0 || Babs <= 0.0)
    //             return;

    //         const double alpha = std::abs(hall_coef) / (ne + ne_floor);
    //         const double dth = CFL_HALL * (h * h) / (alpha * Babs + 1e-300);

    //         if (dth < dt_hall_min_l)
    //         {
    //             dt_hall_min_l = dth;
    //             h_at_min_l = h;
    //             ne_at_min_l = ne;
    //             B_at_min_l = Babs;
    //             alpha_at_min_l = alpha;
    //             ib_at_min_l = ib;
    //             i_at_min_l = i;
    //             j_at_min_l = j;
    //             k_at_min_l = k;
    //             dir_at_min_l = dir;
    //         }

    //         dt_local = std::min(dt_local, dth);
    //     };

    //     const int nb = fld_->num_blocks();
    //     for (int ib = 0; ib < nb; ++ib)
    //     {
    //         auto &UH = fld_->field(fid_.fid_U_H, ib);
    //         auto &UNa = fld_->field(fid_.fid_U_Na, ib);
    //         auto &Bcel = fld_->field(fid_.fid_Bcell, ib);

    //         auto &dl_xi = fld_->field("dl_xi", ib);
    //         auto &dl_eta = fld_->field("dl_eta", ib);
    //         auto &dl_zeta = fld_->field("dl_zeta", ib);

    //         if (!UH.is_allocated() || !UNa.is_allocated() || !Bcel.is_allocated())
    //             continue;
    //         if (!dl_xi.is_allocated() || !dl_eta.is_allocated() || !dl_zeta.is_allocated())
    //             continue;

    //         // EdgeXi: cells (i,j,k),(i,j-1,k),(i,j,k-1),(i,j-1,k-1)
    //         {
    //             Int3 lo = dl_xi.inner_lo();
    //             Int3 hi = dl_xi.inner_hi();
    //             for (int i = lo.i; i < hi.i; ++i)
    //                 for (int j = lo.j; j < hi.j; ++j)
    //                     for (int k = lo.k; k < hi.k; ++k)
    //                     {
    //                         const double h = dl_xi(i, j, k, 0);

    //                         const double ne =
    //                             0.25 * (NUM(UH, UNa, i, j, k) +
    //                                     NUM(UH, UNa, i, j - 1, k) +
    //                                     NUM(UH, UNa, i, j, k - 1) +
    //                                     NUM(UH, UNa, i, j - 1, k - 1));

    //                         const double Babs =
    //                             0.25 * (Babs_cell(Bcel, i, j, k) +
    //                                     Babs_cell(Bcel, i, j - 1, k) +
    //                                     Babs_cell(Bcel, i, j, k - 1) +
    //                                     Babs_cell(Bcel, i, j - 1, k - 1));

    //                         // update_dt_hall(h, ne, Babs);
    //                         update_dt_hall_dbg(h, ne, Babs, ib, i, j, k, 0);
    //                     }
    //         }

    //         // EdgeEta: cells (i,j,k),(i-1,j,k),(i,j,k-1),(i-1,j,k-1)
    //         {
    //             Int3 lo = dl_eta.inner_lo();
    //             Int3 hi = dl_eta.inner_hi();
    //             for (int i = lo.i; i < hi.i; ++i)
    //                 for (int j = lo.j; j < hi.j; ++j)
    //                     for (int k = lo.k; k < hi.k; ++k)
    //                     {
    //                         const double h = dl_eta(i, j, k, 0);

    //                         const double ne =
    //                             0.25 * (NUM(UH, UNa, i, j, k) +
    //                                     NUM(UH, UNa, i - 1, j, k) +
    //                                     NUM(UH, UNa, i, j, k - 1) +
    //                                     NUM(UH, UNa, i - 1, j, k - 1));

    //                         const double Babs =
    //                             0.25 * (Babs_cell(Bcel, i, j, k) +
    //                                     Babs_cell(Bcel, i - 1, j, k) +
    //                                     Babs_cell(Bcel, i, j, k - 1) +
    //                                     Babs_cell(Bcel, i - 1, j, k - 1));

    //                         // update_dt_hall(h, ne, Babs);
    //                         update_dt_hall_dbg(h, ne, Babs, ib, i, j, k, 1);
    //                     }
    //         }

    //         // EdgeZeta: cells (i,j,k),(i-1,j,k),(i,j-1,k),(i-1,j-1,k)
    //         {
    //             Int3 lo = dl_zeta.inner_lo();
    //             Int3 hi = dl_zeta.inner_hi();
    //             for (int i = lo.i; i < hi.i; ++i)
    //                 for (int j = lo.j; j < hi.j; ++j)
    //                     for (int k = lo.k; k < hi.k; ++k)
    //                     {
    //                         const double h = dl_zeta(i, j, k, 0);

    //                         const double ne =
    //                             0.25 * (NUM(UH, UNa, i, j, k) +
    //                                     NUM(UH, UNa, i - 1, j, k) +
    //                                     NUM(UH, UNa, i, j - 1, k) +
    //                                     NUM(UH, UNa, i - 1, j - 1, k));

    //                         const double Babs =
    //                             0.25 * (Babs_cell(Bcel, i, j, k) +
    //                                     Babs_cell(Bcel, i - 1, j, k) +
    //                                     Babs_cell(Bcel, i, j - 1, k) +
    //                                     Babs_cell(Bcel, i - 1, j - 1, k));

    //                         // update_dt_hall(h, ne, Babs);
    //                         update_dt_hall_dbg(h, ne, Babs, ib, i, j, k, 2);
    //                     }
    //         }
    //     }
    // }

    // MPI 全局最小 dt
    double dt_global = dt_local;
    PARALLEL::mpi_min(&dt_local, &dt_global, 1);

    dt = dt_global;

    double dt_mhd_min_global = dt_mhd_min_local;
    double dt_hall_min_global = dt_hall_min_local;
    PARALLEL::mpi_min(&dt_mhd_min_local, &dt_mhd_min_global, 1);
    PARALLEL::mpi_min(&dt_hall_min_local, &dt_hall_min_global, 1);

    if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    {
        std::printf("[dt split] step=%d dt=%.3e  dt_mhd=%.3e  dt_hall=%.3e\n",
                    io_.Run().step, dt, dt_mhd_min_global, dt_hall_min_global);
        std::fflush(stdout);
    }

    dt_hall = dt_hall_min_local;

    // const int myid = par_->GetInt("myid");
    // if (std::abs(dt_hall_min_l - dt_global) <= 1e-12 * (dt_global + 1e-300))
    // {
    //     std::printf("[HallDtMin rank %d] dt=%.3e dir=%d ib=%d (%d,%d,%d)  h=%.3e  ne=%.3e  |B|=%.3e  alpha=%.3e\n",
    //                 myid, dt_hall_min_l, dir_at_min_l, ib_at_min_l, i_at_min_l, j_at_min_l, k_at_min_l,
    //                 h_at_min_l, ne_at_min_l, B_at_min_l, alpha_at_min_l);
    //     std::fflush(stdout);
    // }
}