#include "MercurySolver.h"

void MercurySolver::AddResistiveEdgeEMF_()
{
    const int nb = fld_->num_blocks();

    const double r_cut_in = 0.8;
    const double r_cut_out = 1.0;
    const double r0 = 0.8;
    const double r1 = 1.0;
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

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        // --- EdgeXi: use the edge midpoint of node (i,j,k) -> (i+1,j,k) ---
        {
            Int3 lo = Exi.inner_lo();
            Int3 hi = Exi.inner_hi();

            // 如果你后续在扩散项里不需要 i±1/j±1/k±1，就不必“严格 inner 收缩”；
            // 这里我们只用本点的 Jxi 和本边的坐标，所以不用收缩。
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i + 1, j, k));
                        const double ym = 0.5 * (y(i, j, k) + y(i + 1, j, k));
                        const double zm = 0.5 * (z(i, j, k) + z(i + 1, j, k));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Exi(i, j, k, 0) += (inver_Rem * yita0) * Jxi(i, j, k, 0);
                    }
        }

        // --- EdgeEt: node (i,j,k) -> (i,j+1,k) ---
        {
            Int3 lo = Eeta.inner_lo();
            Int3 hi = Eeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i, j + 1, k));
                        const double ym = 0.5 * (y(i, j, k) + y(i, j + 1, k));
                        const double zm = 0.5 * (z(i, j, k) + z(i, j + 1, k));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Eeta(i, j, k, 0) += (inver_Rem * yita0) * Jeta(i, j, k, 0);
                    }
        }

        // --- EdgeZe: node (i,j,k) -> (i,j,k+1) ---
        {
            Int3 lo = Eze.inner_lo();
            Int3 hi = Eze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i, j, k + 1));
                        const double ym = 0.5 * (y(i, j, k) + y(i, j, k + 1));
                        const double zm = 0.5 * (z(i, j, k) + z(i, j, k + 1));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        Eze(i, j, k, 0) += (inver_Rem * yita0) * Jzeta(i, j, k, 0);
                    }
        }
    }
}