
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

    auto scan_one = [&](int fidU)
    {
        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            FieldBlock &U = fld_->field(fidU, ib);
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

                        double c = std::sqrt(gamma_ * p / rho);

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
                            return (std::abs(un) + c) * A;
                        };

                        double denom = 0.0;

                        // xi+ at i  : outward = Axi(i)
                        denom += face_term(Axi(i, j, k, 0), Axi(i, j, k, 1), Axi(i, j, k, 2), true);
                        // xi- at i+1: outward = -Axi(i+1)
                        denom += face_term(Axi(i - 1, j, k, 0), Axi(i - 1, j, k, 1), Axi(i - 1, j, k, 2), false);

                        // eta+ at j
                        denom += face_term(Aet(i, j, k, 0), Aet(i, j, k, 1), Aet(i, j, k, 2), true);
                        // eta- at j+1
                        denom += face_term(Aet(i, j - 1, k, 0), Aet(i, j - 1, k, 1), Aet(i, j - 1, k, 2), false);

                        // zeta+ at k
                        denom += face_term(Aze(i, j, k, 0), Aze(i, j, k, 1), Aze(i, j, k, 2), true);
                        // zeta- at k+1
                        denom += face_term(Aze(i, j, k - 1, 0), Aze(i, j, k - 1, 1), Aze(i, j, k - 1, 2), false);

                        if (denom > 0.0)
                        {
                            double dtc = CFL * V / denom;
                            dt_local = std::min(dt_local, dtc);
                        }
                    }
        }
    };

    scan_one(fid_.fid_U_H);
    scan_one(fid_.fid_U_Na);

    // MPI 全局最小 dt
    double dt_global = dt_local;
    PARALLEL::mpi_min(&dt_local, &dt_global, 1);

    dt = dt_global;
}