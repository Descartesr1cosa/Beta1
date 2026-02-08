#include "1_Boundary.h"

void MercuryBoundary::BC_UH_Farfield_H_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = bc_state_.qinf[0];
                U(i, j, k, 1) = bc_state_.qinf[1];
                U(i, j, k, 2) = bc_state_.qinf[2];
                U(i, j, k, 3) = bc_state_.qinf[3];
                U(i, j, k, 4) = bc_state_.qinf[4];
            }
}

void MercuryBoundary::BC_UH_Farfield_Na_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{

    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = bc_state_.qinfs[0];
                U(i, j, k, 1) = bc_state_.qinfs[1];
                U(i, j, k, 2) = bc_state_.qinfs[2];
                U(i, j, k, 3) = bc_state_.qinfs[3];
                U(i, j, k, 4) = bc_state_.qinfs[4];
            }
}

void MercuryBoundary::BC_Solid_Surface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write
    const Box3 &inner = r.inner_slab;                                                     // 1-layer inner slab reference

    const int ax = std::abs(r.direction);        // 1/2/3
    const int sgn = (r.direction > 0) ? +1 : -1; // outward normal sign

    auto pressure_from_cons = [&](double rho, double mx, double my, double mz, double E)
    {
        const double rho_floor = 1e-30;
        rho = std::max(rho, rho_floor);
        const double ke = 0.5 * (mx * mx + my * my + mz * mz) / rho;
        double p = (bc_state_.gamma - 1.0) * (E - ke);
        const double p_floor = 1e-30;
        return std::max(p, p_floor);
    };

    auto ref_index = [&](int i, int j, int k, int &ii, int &jj, int &kk, int &d)
    {
        // pick the inner reference cell on the slab (thickness = 1)
        ii = i;
        jj = j;
        kk = k;

        if (ax == 1)
        {
            const int i_ref = (sgn < 0) ? inner.lo.i : (inner.hi.i - 1);
            d = (sgn < 0) ? (i_ref - i) : (i - i_ref);
            ii = i_ref;
        }
        else if (ax == 2)
        {
            const int j_ref = (sgn < 0) ? inner.lo.j : (inner.hi.j - 1);
            d = (sgn < 0) ? (j_ref - j) : (j - j_ref);
            jj = j_ref;
        }
        else
        {
            const int k_ref = (sgn < 0) ? inner.lo.k : (inner.hi.k - 1);
            d = (sgn < 0) ? (k_ref - k) : (k - k_ref);
            kk = k_ref;
        }

        if (d < 1)
            d = 1; // ghost layer depth starts at 1
    };

    const double rho_min_wall = 0.05; // Fortran: if <=0.05 then set to 0.05

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                int ii, jj, kk, d;
                ref_index(i, j, k, ii, jj, kk, d);

                const double rho_ref = U(ii, jj, kk, 0);
                const double mx_ref = U(ii, jj, kk, 1);
                const double my_ref = U(ii, jj, kk, 2);
                const double mz_ref = U(ii, jj, kk, 3);
                const double E_ref = U(ii, jj, kk, 4);

                const double p_ref = pressure_from_cons(rho_ref, mx_ref, my_ref, mz_ref, E_ref);

                // rho = 0.5^d * rho_ref  (generalization of Fortran two layers)
                double rho = rho_ref * std::pow(0.5, double(d));
                if (rho < rho_min_wall)
                    rho = rho_min_wall;

                // no-slip u=v=w=0  => momentum = 0
                U(i, j, k, 0) = rho;
                U(i, j, k, 1) = 0.0;
                U(i, j, k, 2) = 0.0;
                U(i, j, k, 3) = 0.0;

                // energy = p/(gamma-1) + 0.5*rho*0
                U(i, j, k, 4) = p_ref / (bc_state_.gamma - 1.0);
            }
}

void MercuryBoundary::BC_Solid_Surface_Eface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{

    // const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write
    const Box3 &g = r.inner_slab;

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = 0.0;
                U(i, j, k, 1) = 0.0;
                U(i, j, k, 2) = 0.0;
            }
}
