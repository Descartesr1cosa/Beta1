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

void MercuryBoundary::BC_Solid_Surface_(FieldBlock &U, Field *fld,
                                        const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;            // 1-layer inner slab (cell-centered)
    const int ax = std::abs(r.direction);        // 1/2/3
    const int sgn = (r.direction > 0) ? +1 : -1; // +: high face, -: low face (outward)

    // outward step (inner -> ghost)
    const Int3 cyc = (ax == 1)   ? Int3{sgn, 0, 0}
                     : (ax == 2) ? Int3{0, sgn, 0}
                                 : Int3{0, 0, sgn};

    // face shift: high face uses +1, low face uses 0 (same as your reference)
    const int shift = (sgn > 0) ? 1 : 0;

    // ---- identify species from descriptor.name  ----
    const std::string name = U.descriptor().name;
    const bool is_Na = (name.find("U_Na") != std::string::npos) || (name.find("Na") != std::string::npos);

    // ---- floors: keep consistent with your calc_PV floors to avoid U/PV mismatch ----
    const double rho_floor = 1e-12; // same order as calc_PV
    const double p_floor = 1e-12;

    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto scal = [&](const std::array<double, 3> &a, double s)
    {
        return std::array<double, 3>{a[0] * s, a[1] * s, a[2] * s};
    };
    auto sub = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };

    auto pressure_from_cons = [&](double rho, double mx, double my, double mz, double E)
    {
        rho = std::max(rho, rho_floor);
        const double ke = 0.5 * (mx * mx + my * my + mz * mz) / rho;
        double eint = E - ke;
        if (eint < 0.0)
            eint = 0.0;
        double p = (bc_state_.gamma - 1.0) * eint;
        return std::max(p, p_floor);
    };

    // clamp indices into allocated range (including ghost)
    const Int3 ulo = U.get_lo();
    const Int3 uhi = U.get_hi();
    auto clamp_idx = [&](int &i, int &j, int &k)
    {
        i = std::min(std::max(i, ulo.i), uhi.i - 1);
        j = std::min(std::max(j, ulo.j), uhi.j - 1);
        k = std::min(std::max(k, ulo.k), uhi.k - 1);
    };

    // ---- metric face area vector field ----
    const int ib = r.this_block;

    FieldBlock *JD = nullptr;
    if (ax == 1)
        JD = &fld->field("JDxi", ib);
    if (ax == 2)
        JD = &fld->field("JDet", ib);
    if (ax == 3)
        JD = &fld->field("JDze", ib);

    // make normal always point to FLUID side (same logic as your apply_cell_wall)
    const double sign_to_fluid = (sgn > 0) ? -1.0 : +1.0;

    // loop over inner boundary-adjacent cells
    for (int ii = inner.lo.i; ii < inner.hi.i; ++ii)
        for (int jj = inner.lo.j; jj < inner.hi.j; ++jj)
            for (int kk = inner.lo.k; kk < inner.hi.k; ++kk)
            {
                // wall unit normal from metric face area vector at the boundary face
                const int fi = ii + cyc.i * shift;
                const int fj = jj + cyc.j * shift;
                const int fk = kk + cyc.k * shift;

                std::array<double, 3> A = {(*JD)(fi, fj, fk, 0),
                                           (*JD)(fi, fj, fk, 1),
                                           (*JD)(fi, fj, fk, 2)};

                const double An2 = std::max(dot(A, A), 1e-30);
                const double An = std::sqrt(An2);
                std::array<double, 3> n_hat = scal(A, sign_to_fluid / An); // points into FLUID

                // fill each ghost layer: ghost index = inner + ng*cyc
                // reference index = inner - (ng-1)*cyc (keeps extension smooth for multi-ghost)
                for (int ng = 1; ng <= ngh; ++ng)
                {
                    int ir = ii - (ng - 1) * cyc.i;
                    int jr = jj - (ng - 1) * cyc.j;
                    int kr = kk - (ng - 1) * cyc.k;
                    clamp_idx(ir, jr, kr);

                    double rho_r = std::max(U(ir, jr, kr, 0), rho_floor);
                    double mx_r = U(ir, jr, kr, 1);
                    double my_r = U(ir, jr, kr, 2);
                    double mz_r = U(ir, jr, kr, 3);
                    double E_r = U(ir, jr, kr, 4);

                    double p_r = pressure_from_cons(rho_r, mx_r, my_r, mz_r, E_r);

                    std::array<double, 3> v = {mx_r / rho_r, my_r / rho_r, mz_r / rho_r};
                    const double vn = dot(v, n_hat);

                    // diode: block outflow from surface to fluid (vn>0)
                    std::array<double, 3> v_g = v;
                    v_g = sub(v, scal(n_hat, vn));
                    // if (vn > 0.0)
                    // {
                    //     // remove only the outward-normal component; tangential unchanged
                    //     v_g = sub(v, scal(n_hat, vn)); // vn_g = 0
                    //     // If you want a *stronger sink* like some papers do, uncomment:
                    //     // rho_r = rho_floor;
                    //     // p_r   = p_floor;
                    // }
                    // v_g = sub(v, scal(n_hat, 2.0 * vn));
                    // v_g = {0.0, 0.0, 0.0};

                    int ig = ii + ng * cyc.i;
                    int jg = jj + ng * cyc.j;
                    int kg = kk + ng * cyc.k;
                    clamp_idx(ig, jg, kg);

                    // write ghost
                    U(ig, jg, kg, 0) = rho_r;
                    U(ig, jg, kg, 1) = rho_r * v_g[0];
                    U(ig, jg, kg, 2) = rho_r * v_g[1];
                    U(ig, jg, kg, 3) = rho_r * v_g[2];

                    // energy consistent (Mercury U_H/U_Na excludes magnetic energy)
                    const double ke_g = 0.5 * rho_r * dot(v_g, v_g);
                    U(ig, jg, kg, 4) = p_r / (bc_state_.gamma - 1.0) + ke_g;
                }
            }
}

// void MercuryBoundary::BC_Solid_Surface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
// {
//     const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write
//     const Box3 &inner = r.inner_slab;                                                     // 1-layer inner slab reference

//     const int ax = std::abs(r.direction);        // 1/2/3
//     const int sgn = (r.direction > 0) ? +1 : -1; // outward normal sign

//     std::string fname = U.descriptor().name;
//     const bool is_Na = (fname == "U_Na");

//     const double rho_floor_global = 1e-12;
//     const double p_floor_global = 1e-12;

//     auto pressure_from_cons = [&](double rho, double mx, double my, double mz, double E)
//     {
//         rho = std::max(rho, rho_floor_global);
//         const double ke = 0.5 * (mx * mx + my * my + mz * mz) / rho;
//         double p = (bc_state_.gamma - 1.0) * (E - ke);
//         return std::max(p, p_floor_global);
//     };

//     auto ref_index = [&](int i, int j, int k, int &ii, int &jj, int &kk, int &d)
//     {
//         // pick the inner reference cell on the slab (thickness = 1)
//         ii = i;
//         jj = j;
//         kk = k;

//         if (ax == 1)
//         {
//             const int i_ref = (sgn < 0) ? inner.lo.i : (inner.hi.i - 1);
//             d = (sgn < 0) ? (i_ref - i) : (i - i_ref);
//             ii = i_ref;
//         }
//         else if (ax == 2)
//         {
//             const int j_ref = (sgn < 0) ? inner.lo.j : (inner.hi.j - 1);
//             d = (sgn < 0) ? (j_ref - j) : (j - j_ref);
//             jj = j_ref;
//         }
//         else
//         {
//             const int k_ref = (sgn < 0) ? inner.lo.k : (inner.hi.k - 1);
//             d = (sgn < 0) ? (k_ref - k) : (k - k_ref);
//             kk = k_ref;
//         }

//         if (d < 1)
//             d = 1; // ghost layer depth starts at 1
//     };

//     auto mirror_ref_index = [&](int i, int j, int k, int &ii, int &jj, int &kk, int &d)
//     {
//         ii = i;
//         jj = j;
//         kk = k;
//         d = 1;

//         if (ax == 1)
//         {
//             const int i0 = (sgn < 0) ? inner.lo.i : (inner.hi.i - 1); // 第一层内点（通常是0或hi-1）
//             d = (sgn < 0) ? (i0 - i) : (i - i0);                      // ghost 深度：-1 =>1, -2=>2 ...
//             if (d < 1)
//                 d = 1;
//             ii = i0 + (-sgn) * (d - 1); // 往域内走 d-1 格做镜像参考
//         }
//         else if (ax == 2)
//         {
//             const int j0 = (sgn < 0) ? inner.lo.j : (inner.hi.j - 1);
//             d = (sgn < 0) ? (j0 - j) : (j - j0);
//             if (d < 1)
//                 d = 1;
//             jj = j0 + (-sgn) * (d - 1);
//         }
//         else // ax==3
//         {
//             const int k0 = (sgn < 0) ? inner.lo.k : (inner.hi.k - 1);
//             d = (sgn < 0) ? (k0 - k) : (k - k0);
//             if (d < 1)
//                 d = 1;
//             kk = k0 + (-sgn) * (d - 1);
//         }

//         // 安全夹紧（防止 ngh 过大或极小 block）
//         const Int3 ulo = U.get_lo();
//         const Int3 uhi = U.get_hi();
//         ii = std::min(std::max(ii, ulo.i), uhi.i - 1);
//         jj = std::min(std::max(jj, ulo.j), uhi.j - 1);
//         kk = std::min(std::max(kk, ulo.k), uhi.k - 1);
//     };

//     for (int i = g.lo.i; i < g.hi.i; ++i)
//         for (int j = g.lo.j; j < g.hi.j; ++j)
//             for (int k = g.lo.k; k < g.hi.k; ++k)
//             {
//                 int ii, jj, kk, d;
//                 // ref_index(i, j, k, ii, jj, kk, d);
//                 mirror_ref_index(i, j, k, ii, jj, kk, d);

//                 double rho_ref = U(ii, jj, kk, 0);
//                 const double mx_ref = U(ii, jj, kk, 1);
//                 const double my_ref = U(ii, jj, kk, 2);
//                 const double mz_ref = U(ii, jj, kk, 3);
//                 const double E_ref = U(ii, jj, kk, 4);

//                 rho_ref = std::max(rho_ref, rho_floor_global);
//                 const double p_ref = pressure_from_cons(rho_ref, mx_ref, my_ref, mz_ref, E_ref);

//                 // Slip wall: tangential same, normal antisymmetric (odd reflection)
//                 double mx_g = mx_ref, my_g = my_ref, mz_g = mz_ref;
//                 if (ax == 1)
//                     mx_g = -mx_r;
//                 if (ax == 2)
//                     my_g = -my_r;
//                 if (ax == 3)
//                     mz_g = -mz_r;

//                 // Energy consistent with (rho_g, m_g, p_g)
//                 const double rho_g = rho_ref;
//                 const double ke_g = 0.5 * (mx_g * mx_g + my_g * my_g + mz_g * mz_g) / rho_g;
//                 const double E_g = p_ref / (bc_state_.gamma - 1.0) + ke_g;

//                 // write ghost
//                 U(i, j, k, 0) = rho_g;
//                 U(i, j, k, 1) = mx_g;
//                 U(i, j, k, 2) = my_g;
//                 U(i, j, k, 3) = mz_g;
//                 U(i, j, k, 4) = E_g;
//             }
// }

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

void MercuryBoundary::BC_Solid_Surface_Eedge_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{

    // const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write
    const Box3 &g = r.inner_slab;

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = 0.0;
            }
}

void MercuryBoundary::BC_Pole_Eedge_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            double temp_E = 0.0;
            double num_d = 0.0;

            for (int k = inner.lo.k; k < inner.hi.k - 1; ++k) // avoid recounting: k=0 <=> k=inner.hi.k-1
            {
                temp_E += U(i, j, k, 0);
                num_d += 1.0;
            }
            temp_E /= num_d;

            for (int k = inner.lo.k; k < inner.hi.k; ++k)
                U(i, j, k, 0) = temp_E;
        }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Pole_Eedge_Zero(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
                U(i, j, k, 0) = 0.0;
        }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Pole_Cell_(FieldBlock &U, Field *fld,
                                    const BOUND::PhysicalRegion &r, int ngh)
{
    if (!U.is_allocated())
        return;
    const Box3 &inner = r.inner_slab;

    // 这里把 ncomp() 换成你 FieldBlock 实际的分量数接口
    const int ncomp = U.descriptor().ncomp;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
            for (int n = 0; n < ncomp; ++n)
            {
                double temp_U = 0.0;
                double num_d = 0.0;

                // 若 k=0 与 k=inner.hi.k-1 物理等价，则避免重复计数
                for (int k = inner.lo.k; k < inner.hi.k - 1; ++k)
                {
                    temp_U += U(i, j, k, n);
                    num_d += 1.0;
                }

                if (num_d > 0.0)
                    temp_U /= num_d;
                else
                    temp_U = 0.0;

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    U(i, j, k, n) = temp_U;
            }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Pole_Eedge_RegulateKAndCopyGhost_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;

    // 当前 handler 只应被用于：
    //   E_xi  with ax == 1
    //   E_eta with ax == 2
    // 若以后扩展 E_zeta，也可 ax == 3。
    Int3 nrm{0, 0, 0};
    if (ax == 1)
        nrm = Int3{sgn, 0, 0};
    else if (ax == 2)
        nrm = Int3{0, sgn, 0};
    else
        nrm = Int3{0, 0, sgn};

    const Int3 ulo = U.get_lo();
    const Int3 uhi = U.get_hi();

    auto in_range = [&](int i, int j, int k) -> bool
    {
        return (i >= ulo.i && i < uhi.i &&
                j >= ulo.j && j < uhi.j &&
                k >= ulo.k && k < uhi.k);
    };

    // ------------------------------------------------------------
    // 取当前 edge family 对应的物理边向量 dr。
    //
    // E_xi   -> dr_xi
    // E_eta  -> dr_eta
    // E_zeta -> dr_zeta
    // ------------------------------------------------------------

    const int ib = r.this_block;

    FieldBlock *dr_ptr = nullptr;

    if (ax == 1)
        dr_ptr = &((*fld).field("dr_xi")[ib]);
    else if (ax == 2)
        dr_ptr = &((*fld).field("dr_eta")[ib]);
    else
        dr_ptr = &((*fld).field("dr_zeta")[ib]);

    FieldBlock &dr = *dr_ptr;

    auto get_dr = [&](int i, int j, int k) -> std::array<double, 3>
    {
        return {dr(i, j, k, 0), dr(i, j, k, 1), dr(i, j, k, 2)};
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

    auto solve3x3 = [](double A[3][3], double b[3], double x[3]) -> bool
    {
        double M[3][4] = {
            {A[0][0], A[0][1], A[0][2], b[0]},
            {A[1][0], A[1][1], A[1][2], b[1]},
            {A[2][0], A[2][1], A[2][2], b[2]}};

        for (int c = 0; c < 3; ++c)
        {
            int piv = c;
            double amax = std::abs(M[c][c]);

            for (int r0 = c + 1; r0 < 3; ++r0)
            {
                const double av = std::abs(M[r0][c]);
                if (av > amax)
                {
                    amax = av;
                    piv = r0;
                }
            }

            if (amax < 1.0e-300)
                return false;

            if (piv != c)
            {
                for (int q = c; q < 4; ++q)
                    std::swap(M[c][q], M[piv][q]);
            }

            const double inv_piv = 1.0 / M[c][c];

            for (int q = c; q < 4; ++q)
                M[c][q] *= inv_piv;

            for (int r0 = 0; r0 < 3; ++r0)
            {
                if (r0 == c)
                    continue;

                const double fac = M[r0][c];

                for (int q = c; q < 4; ++q)
                    M[r0][q] -= fac * M[c][q];
            }
        }

        x[0] = M[0][3];
        x[1] = M[1][3];
        x[2] = M[2][3];

        return true;
    };

    constexpr double eps_len = 1.0e-14;
    constexpr double eps_reg = 1.0e-12;

    // ============================================================
    // 1. 对第一层物理边界做：
    //
    //      edge 1-form values -> unique Cartesian E vector -> edge 1-form values
    //
    // 对每个固定的 (i,j)，收集所有 k 上的当前 edge：
    //
    //      e_k = U(i,j,k,0)
    //      dr_k = dr_axis(i,j,k,:)
    //
    // 解：
    //
    //      e_k / |dr_k| ≈ E · (dr_k / |dr_k|)
    //
    // 然后：
    //
    //      U(i,j,k,0) = E · dr_k
    //
    // 这样保证同一个 Pole 点的一圈 edge 由唯一 E_vec 生成。
    // ============================================================

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
    {
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            double A[3][3] = {
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0}};

            double b[3] = {0.0, 0.0, 0.0};

            int cnt = 0;

            for (int k = inner.lo.k; k < inner.hi.k - 1; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const auto dl = get_dr(i, j, k);
                const double L = norm3(dl);

                // 对 collapsed edge，不参与矢量重构。
                // 后面会投影成 0。
                if (L < eps_len)
                    continue;

                const double invL = 1.0 / L;

                const double tau[3] = {
                    dl[0] * invL,
                    dl[1] * invL,
                    dl[2] * invL};

                const double y = U(i, j, k, 0) * invL;

                for (int a = 0; a < 3; ++a)
                {
                    b[a] += y * tau[a];

                    for (int c = 0; c < 3; ++c)
                        A[a][c] += tau[a] * tau[c];
                }

                ++cnt;
            }

            if (cnt <= 0)
            {
                // 所有 edge 都 collapsed，则直接置零这一圈。
                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (in_range(i, j, k))
                        U(i, j, k, 0) = 0.0;
                }
                continue;
            }

            const double trA = A[0][0] + A[1][1] + A[2][2];
            const double reg = eps_reg * std::max(1.0, trA);

            A[0][0] += reg;
            A[1][1] += reg;
            A[2][2] += reg;

            double Evec[3] = {0.0, 0.0, 0.0};

            const bool ok = solve3x3(A, b, Evec);

            if (!ok)
            {
                // 极端退化时，退回到代数平均，但仍然按 1-form 处理。
                // 这里不建议直接保留原值，因为原值可能包含周向噪声。
                double avg = 0.0;
                int navg = 0;

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (!in_range(i, j, k))
                        continue;

                    avg += U(i, j, k, 0);
                    ++navg;
                }

                if (navg > 0)
                    avg /= static_cast<double>(navg);

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (in_range(i, j, k))
                        U(i, j, k, 0) = avg;
                }

                continue;
            }

            // 用唯一 Evec 重新投影回每个 edge 的 1-form。
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const auto dl = get_dr(i, j, k);
                const double L = norm3(dl);

                if (L < eps_len)
                {
                    U(i, j, k, 0) = 0.0;
                }
                else
                {
                    U(i, j, k, 0) =
                        Evec[0] * dl[0] +
                        Evec[1] * dl[1] +
                        Evec[2] * dl[2];
                }
            }
        }
    }

    // ============================================================
    // 2. ghost 从 regulation 后的第一层 copy。
    //
    // 注意：
    // 这里 copy 的是已经由 Evec 投影回来的 1-form。
    // ============================================================

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
    {
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const double val = U(i, j, k, 0);

                for (int g = 1; g <= ngh; ++g)
                {
                    const int ig = i + g * nrm.i;
                    const int jg = j + g * nrm.j;
                    const int kg = k + g * nrm.k;

                    if (!in_range(ig, jg, kg))
                        continue;

                    U(ig, jg, kg, 0) = val;
                }
            }
        }
    }
}
