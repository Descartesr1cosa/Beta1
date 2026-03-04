
#include <cmath>
#include <algorithm>
#include <iostream>

#include "MercurySolver.h"

void MercurySolver::Compute_Timestep()
{
    const double rho_floor = 1e-2;
    const double p_floor = 1e-8;

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

    // double dt_mhd_temp_H;
    // PARALLEL::mpi_min(&dt_mhd_min_local, &dt_mhd_temp_H, 1);
    // if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    // {
    //     std::printf("[MHDdtGlobal] step=%d  dt_min_H=%.3e \n",
    //                 io_.Run().step, dt_mhd_temp_H);
    //     std::fflush(stdout);
    // }
    // scan_one(fid_.fid_U_Na);

    // double dt_mhd_temp_H_NA;
    // PARALLEL::mpi_min(&dt_mhd_min_local, &dt_mhd_temp_H_NA, 1);
    // if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    // {
    //     std::printf("[MHDdtGlobal] step=%d  dt_min_H=%.3e \n",
    //                 io_.Run().step, dt_mhd_temp_H_NA);
    //     std::fflush(stdout);
    // }

    if (std::abs(hall_coef) > 0.0)
    {
        const double CFL_HALL = 1.0;
        // const double ne_floor = ne_hall_floor; // same "prime" unit as your rho/M
        const double h_eps = 1e-12;
        const double B_floor = 1e-30;

        double dt_hall_min_l = 1e100;

        auto NUM_cell = [&](FieldBlock &UH, FieldBlock &UNa, int i, int j, int k) -> double
        {
            const double rhoH = std::max(0.0, UH(i, j, k, 0));
            const double rhoN = std::max(0.0, UNa(i, j, k, 0));
            return rhoH / M_H + rhoN / M_Na; // this is your n_e' (mol/kg), consistent with hall_coef
        };

        auto Bcell_vec = [&](FieldBlock &Bcel, int i, int j, int k, double &bx, double &by, double &bz)
        {
            bx = Bcel(i, j, k, 0);
            by = Bcel(i, j, k, 1);
            bz = Bcel(i, j, k, 2); // you said Bcell already includes Badd
        };

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &UH = fld_->field(fid_.fid_U_H, ib);
            auto &UNa = fld_->field(fid_.fid_U_Na, ib);
            auto &Bcel = fld_->field(fid_.fid_Bcell, ib);

            auto &dlx = fld_->field("dl_xi", ib);
            auto &dle = fld_->field("dl_eta", ib);
            auto &dlz = fld_->field("dl_zeta", ib); // may be absent

            if (!UH.is_allocated() || !UNa.is_allocated() || !Bcel.is_allocated())
                continue;
            if (!dlx.is_allocated() || !dle.is_allocated())
                continue;

            auto hmin2 = [&](int i, int j, int k) -> double
            {
                double hx = dlx.is_allocated() ? dlx(i, j, k, 0) : 1e100;
                double he = dle.is_allocated() ? dle(i, j, k, 0) : 1e100;
                double hz = dlz.is_allocated() ? dlz(i, j, k, 0) : 1e100;
                double h = std::min(hx, std::min(he, hz));
                if (h <= h_eps)
                    return 0.0;
                return h * h;
            };

            // ---------- EdgeXi loop (use dlx index space, but h uses local min) ----------
            {
                Int3 lo = dlx.inner_lo(), hi = dlx.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double h2 = hmin2(i, j, k);
                            if (h2 <= 0.0)
                                continue;

                            // ne' at xi-edge: avg 4 surrounding cells in (j,k)
                            const double ne = 0.25 * (NUM_cell(UH, UNa, i, j, k) + NUM_cell(UH, UNa, i, j - 1, k) + NUM_cell(UH, UNa, i, j, k - 1) + NUM_cell(UH, UNa, i, j - 1, k - 1));
                            // const double ne_eff = ne + ne_floor;
                            const double ne_true = ne;
                            // 1) 平滑 floor（避免 max 的硬拐点）
                            const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                            // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                            const double s = ne_true / (ne_true + ne_hall_cut);

                            // B at same edge location: avg same 4 cells (vector)
                            double bx0, by0, bz0, bx1, by1, bz1, bx2, by2, bz2, bx3, by3, bz3;
                            Bcell_vec(Bcel, i, j, k, bx0, by0, bz0);
                            Bcell_vec(Bcel, i, j - 1, k, bx1, by1, bz1);
                            Bcell_vec(Bcel, i, j, k - 1, bx2, by2, bz2);
                            Bcell_vec(Bcel, i, j - 1, k - 1, bx3, by3, bz3);
                            const double bx = 0.25 * (bx0 + bx1 + bx2 + bx3);
                            const double by = 0.25 * (by0 + by1 + by2 + by3);
                            const double bz = 0.25 * (bz0 + bz1 + bz2 + bz3);
                            const double Babs = std::sqrt(bx * bx + by * by + bz * bz);
                            if (Babs <= B_floor)
                                continue;

                            // dt_hall' = CFL * h^2 * (ne'+floor) / (|hall_coef| * |B|)
                            const double dt_try = CFL_HALL * h2 / s * ne_eff / (std::abs(hall_coef) * Babs + 1e-300);
                            dt_hall_min_l = std::min(dt_hall_min_l, dt_try);
                        }
            }

            // ---------- EdgeEta loop ----------
            {
                Int3 lo = dle.inner_lo(), hi = dle.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double h2 = hmin2(i, j, k);
                            if (h2 <= 0.0)
                                continue;

                            const double ne = 0.25 * (NUM_cell(UH, UNa, i, j, k) + NUM_cell(UH, UNa, i - 1, j, k) + NUM_cell(UH, UNa, i, j, k - 1) + NUM_cell(UH, UNa, i - 1, j, k - 1));
                            // const double ne_eff = ne + ne_floor;
                            const double ne_true = ne;
                            // 1) 平滑 floor（避免 max 的硬拐点）
                            const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                            // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                            const double s = ne_true / (ne_true + ne_hall_cut);

                            double bx0, by0, bz0, bx1, by1, bz1, bx2, by2, bz2, bx3, by3, bz3;
                            Bcell_vec(Bcel, i, j, k, bx0, by0, bz0);
                            Bcell_vec(Bcel, i - 1, j, k, bx1, by1, bz1);
                            Bcell_vec(Bcel, i, j, k - 1, bx2, by2, bz2);
                            Bcell_vec(Bcel, i - 1, j, k - 1, bx3, by3, bz3);
                            const double bx = 0.25 * (bx0 + bx1 + bx2 + bx3);
                            const double by = 0.25 * (by0 + by1 + by2 + by3);
                            const double bz = 0.25 * (bz0 + bz1 + bz2 + bz3);
                            const double Babs = std::sqrt(bx * bx + by * by + bz * bz);
                            if (Babs <= B_floor)
                                continue;

                            const double dt_try = CFL_HALL * h2 / s * ne_eff / (std::abs(hall_coef) * Babs + 1e-300);
                            dt_hall_min_l = std::min(dt_hall_min_l, dt_try);
                        }
            }

            // ---------- EdgeZeta loop (only if dlz exists) ----------
            if (dlz.is_allocated())
            {
                Int3 lo = dlz.inner_lo(), hi = dlz.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double h2 = hmin2(i, j, k);
                            if (h2 <= 0.0)
                                continue;

                            const double ne = 0.25 * (NUM_cell(UH, UNa, i, j, k) + NUM_cell(UH, UNa, i - 1, j, k) + NUM_cell(UH, UNa, i, j - 1, k) + NUM_cell(UH, UNa, i - 1, j - 1, k));
                            // const double ne_eff = ne + ne_floor;
                            const double ne_true = ne;
                            // 1) 平滑 floor（避免 max 的硬拐点）
                            const double ne_eff = std::sqrt(ne_true * ne_true + ne_hall_floor * ne_hall_floor);
                            // 2) 平滑 taper（替代 hard cut；ne_cut_hall 控制过渡宽度）
                            const double s = ne_true / (ne_true + ne_hall_cut);

                            double bx0, by0, bz0, bx1, by1, bz1, bx2, by2, bz2, bx3, by3, bz3;
                            Bcell_vec(Bcel, i, j, k, bx0, by0, bz0);
                            Bcell_vec(Bcel, i - 1, j, k, bx1, by1, bz1);
                            Bcell_vec(Bcel, i, j - 1, k, bx2, by2, bz2);
                            Bcell_vec(Bcel, i - 1, j - 1, k, bx3, by3, bz3);
                            const double bx = 0.25 * (bx0 + bx1 + bx2 + bx3);
                            const double by = 0.25 * (by0 + by1 + by2 + by3);
                            const double bz = 0.25 * (bz0 + bz1 + bz2 + bz3);
                            const double Babs = std::sqrt(bx * bx + by * by + bz * bz);
                            if (Babs <= B_floor)
                                continue;

                            const double dt_try = CFL_HALL * h2 / s * ne_eff / (std::abs(hall_coef) * Babs + 1e-300);
                            dt_hall_min_l = std::min(dt_hall_min_l, dt_try);
                        }
            }
        }

        double dt_hall_min_g = dt_hall_min_l;
        PARALLEL::mpi_min(&dt_hall_min_l, &dt_hall_min_g, 1);

        dt_hall_min_local = std::min(dt_hall_min_local, dt_hall_min_g);
    }

    // MPI 全局最小 dt
    double dt_global = dt_local;
    PARALLEL::mpi_min(&dt_local, &dt_global, 1);

    dt = dt_global;

    // double dt_mhd_min_global = dt_mhd_min_local;
    double dt_hall_min_global = dt_hall_min_local;
    // PARALLEL::mpi_min(&dt_mhd_min_local, &dt_mhd_min_global, 1);
    PARALLEL::mpi_min(&dt_hall_min_local, &dt_hall_min_global, 1);

    // if (par_->GetInt("myid") == 0 && (run_data_->step % par_->GetInt("output_residual") == 0))
    // {
    //     std::printf("[dt split] step=%d dt=%.3e  dt_mhd=%.3e  dt_hall=%.3e\n",
    //                 io_.Run().step, dt, dt_mhd_min_global, dt_hall_min_global);
    //     std::fflush(stdout);
    // }

    dt_hall = dt_hall_min_global;

    // const int myid = par_->GetInt("myid");
    // if (std::abs(dt_hall_min_l - dt_global) <= 1e-12 * (dt_global + 1e-300))
    // {
    //     std::printf("[HallDtMin rank %d] dt=%.3e dir=%d ib=%d (%d,%d,%d)  h=%.3e  ne=%.3e  |B|=%.3e  alpha=%.3e\n",
    //                 myid, dt_hall_min_l, dir_at_min_l, ib_at_min_l, i_at_min_l, j_at_min_l, k_at_min_l,
    //                 h_at_min_l, ne_at_min_l, B_at_min_l, alpha_at_min_l);
    //     std::fflush(stdout);
    // }
}