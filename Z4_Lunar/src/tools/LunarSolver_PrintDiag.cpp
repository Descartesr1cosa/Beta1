#include <limits>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "LunarSolver.h"

void LunarSolver::PrintMinMaxDiagnostics_()
{
    double rho_min = std::numeric_limits<double>::infinity();
    double rho_max = -std::numeric_limits<double>::infinity();
    double p_min = std::numeric_limits<double>::infinity();
    double p_max = -std::numeric_limits<double>::infinity();
    double bx_min = std::numeric_limits<double>::infinity();
    double bx_max = -std::numeric_limits<double>::infinity();
    double by_min = std::numeric_limits<double>::infinity();
    double by_max = -std::numeric_limits<double>::infinity();
    double bz_min = std::numeric_limits<double>::infinity();
    double bz_max = -std::numeric_limits<double>::infinity();
    double b2_min = std::numeric_limits<double>::infinity();
    double b2_max = -std::numeric_limits<double>::infinity();
    double divb_max = 0.0;

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        FieldBlock &U = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &PV = fld_->field(fid_.fid_PV_H, ib);
        FieldBlock &B = fld_->field(fid_.fid_Bcell, ib);
        FieldBlock &divB = fld_->field(fid_.fid_divB, ib);
        if (!U.is_allocated() || !PV.is_allocated() ||
            !B.is_allocated() || !divB.is_allocated())
            continue;

        const Int3 lo = U.inner_lo();
        const Int3 hi = U.inner_hi();
        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rho = U(i, j, k, 0);
                    const double p = PV(i, j, k, 3);
                    const double bx = B(i, j, k, 0);
                    const double by = B(i, j, k, 1);
                    const double bz = B(i, j, k, 2);
                    const double b2 = 0.5 * (bx * bx + by * by + bz * bz);

                    rho_min = std::min(rho_min, rho);
                    rho_max = std::max(rho_max, rho);
                    p_min = std::min(p_min, p);
                    p_max = std::max(p_max, p);
                    bx_min = std::min(bx_min, bx);
                    bx_max = std::max(bx_max, bx);
                    by_min = std::min(by_min, by);
                    by_max = std::max(by_max, by);
                    bz_min = std::min(bz_min, bz);
                    bz_max = std::max(bz_max, bz);
                    b2_min = std::min(b2_min, b2);
                    b2_max = std::max(b2_max, b2);
                    divb_max = std::max(divb_max, std::abs(divB(i, j, k, 0)));
                }
    }

    double mins_local[6] = {rho_min, p_min, bx_min, by_min, bz_min, b2_min};
    double maxs_local[7] = {rho_max, p_max, bx_max, by_max, bz_max, b2_max, divb_max};
    double mins[6], maxs[7];
    PARALLEL::mpi_min(mins_local, mins, 6);
    PARALLEL::mpi_max(maxs_local, maxs, 7);

    if (par_->GetInt("myid") == 0)
    {
        const double b_abs_max = std::sqrt(std::max(0.0, 2.0 * maxs[5]));
        std::printf("\n           rho=[%.3e, %.3e]  p=[%.3e, %.3e]\n",
                    mins[0], maxs[0], mins[1], maxs[1]);
        std::printf("           Bx =[%.3e, %.3e]  By=[%.3e, %.3e]  Bz=[%.3e, %.3e]\n",
                    mins[2], maxs[2], mins[3], maxs[3], mins[4], maxs[4]);
        std::printf("           Pmag=[%.3e, %.3e]  |B|max=%.3e  divB_max=%.3e\n\n",
                    mins[5], maxs[5], b_abs_max, maxs[6]);
        std::fflush(stdout);
    }
}
