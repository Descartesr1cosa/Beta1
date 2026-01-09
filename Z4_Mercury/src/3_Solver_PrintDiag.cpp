#include <limits>
#include <algorithm>
#include <cstdio>
#include "MercurySolver.h"

// 统计 rho/p (H, Na) 与 Bx/By/Bz 的范围，并打印（rank0）
void MercurySolver::PrintMinMaxDiagnostics_()
{
    // --- local init ---
    double rhoH_min_l = std::numeric_limits<double>::infinity();
    double rhoH_max_l = -std::numeric_limits<double>::infinity();
    double pH_min_l = std::numeric_limits<double>::infinity();
    double pH_max_l = -std::numeric_limits<double>::infinity();

    double rhoNa_min_l = std::numeric_limits<double>::infinity();
    double rhoNa_max_l = -std::numeric_limits<double>::infinity();
    double pNa_min_l = std::numeric_limits<double>::infinity();
    double pNa_max_l = -std::numeric_limits<double>::infinity();

    double bx_min_l = std::numeric_limits<double>::infinity();
    double bx_max_l = -std::numeric_limits<double>::infinity();
    double by_min_l = std::numeric_limits<double>::infinity();
    double by_max_l = -std::numeric_limits<double>::infinity();
    double bz_min_l = std::numeric_limits<double>::infinity();
    double bz_max_l = -std::numeric_limits<double>::infinity();
    double b2_min_l = std::numeric_limits<double>::infinity();
    double b2_max_l = -std::numeric_limits<double>::infinity();

    const int nblock = fld_->num_blocks();
    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);
        auto &PVH = fld_->field(fid_.fid_PV_H, ib);
        auto &PVNa = fld_->field(fid_.fid_PV_Na, ib);
        auto &Bcel = fld_->field(fid_.fid_Bcell, ib);

        if (!UH.is_allocated() || !UNa.is_allocated() || !PVH.is_allocated() || !PVNa.is_allocated() || !Bcel.is_allocated())
            continue;

        const Int3 lo = UH.inner_lo();
        const Int3 hi = UH.inner_hi();

        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rhoH = UH(i, j, k, 0);
                    const double pH = PVH(i, j, k, 3);

                    const double rhoNa = UNa(i, j, k, 0);
                    const double pNa = PVNa(i, j, k, 3);

                    const double bx = Bcel(i, j, k, 0);
                    const double by = Bcel(i, j, k, 1);
                    const double bz = Bcel(i, j, k, 2);

                    rhoH_min_l = std::min(rhoH_min_l, rhoH);
                    rhoH_max_l = std::max(rhoH_max_l, rhoH);
                    pH_min_l = std::min(pH_min_l, pH);
                    pH_max_l = std::max(pH_max_l, pH);

                    rhoNa_min_l = std::min(rhoNa_min_l, rhoNa);
                    rhoNa_max_l = std::max(rhoNa_max_l, rhoNa);
                    pNa_min_l = std::min(pNa_min_l, pNa);
                    pNa_max_l = std::max(pNa_max_l, pNa);

                    bx_min_l = std::min(bx_min_l, bx);
                    bx_max_l = std::max(bx_max_l, bx);
                    by_min_l = std::min(by_min_l, by);
                    by_max_l = std::max(by_max_l, by);
                    bz_min_l = std::min(bz_min_l, bz);
                    bz_max_l = std::max(bz_max_l, bz);

                    b2_min_l = std::min(b2_min_l, 0.5 * (bx * bx + by * by + bz * bz));
                    b2_max_l = std::max(b2_max_l, 0.5 * (bx * bx + by * by + bz * bz));
                }
    }

    // --- MPI reduction ---
    double mins_l[8] = {rhoH_min_l, pH_min_l, rhoNa_min_l, pNa_min_l, bx_min_l, by_min_l, bz_min_l, b2_min_l};
    double mins_g[8];

    double maxs_l[8] = {rhoH_max_l, pH_max_l, rhoNa_max_l, pNa_max_l, bx_max_l, by_max_l, bz_max_l, b2_max_l};
    double maxs_g[8];

    PARALLEL::mpi_min(mins_l, mins_g, 8);
    PARALLEL::mpi_max(maxs_l, maxs_g, 8);

    const int myid = par_->GetInt("myid");
    if (myid == 0)
    {
        // 与输出对齐："[  Diag  ] " 长度是 11，所以这里用 11 个空格缩进
        // std::printf("           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]  rhoNa=[%.3e, %.3e]  pNa=[%.3e, %.3e]\n",
        //             mins_g[0], maxs_g[0], mins_g[1], maxs_g[1], mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);
        std::printf("\n           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]\n",
                    mins_g[0], maxs_g[0], mins_g[1], maxs_g[1]);
        std::printf("           rhoN=[%.3e, %.3e]  pN=[%.3e, %.3e]\n",
                    mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);

        std::printf("           Bx  =[%.3e, %.3e]  By=[%.3e, %.3e]  Bz   =[%.3e, %.3e]  \n",
                    mins_g[4], maxs_g[4], mins_g[5], maxs_g[5], mins_g[6], maxs_g[6]);
        std::printf("           Bmag=[%.3e, %.3e]\n\n\n",
                    mins_g[7], maxs_g[7]);

        std::fflush(stdout);
    }
}
