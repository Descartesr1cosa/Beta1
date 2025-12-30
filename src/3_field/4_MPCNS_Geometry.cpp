#include "3_field/2_MPCNS_Field.h"
#include "array"

void Field::build_geometry()
{
    register_field(FieldDescriptor{"Jac", StaggerLocation::Cell, 1, 0});
    register_field(FieldDescriptor{"JDxi", StaggerLocation::FaceXi, 3, 1});
    register_field(FieldDescriptor{"JDet", StaggerLocation::FaceEt, 3, 1});
    register_field(FieldDescriptor{"JDze", StaggerLocation::FaceZe, 3, 1});

    // --- NEW: covariant basis vectors at cell centers ---
    register_field(FieldDescriptor{"a_xi", StaggerLocation::Cell, 3, 0});
    register_field(FieldDescriptor{"a_eta", StaggerLocation::Cell, 3, 0});
    register_field(FieldDescriptor{"a_zeta", StaggerLocation::Cell, 3, 0});

    register_field(FieldDescriptor{"pinvGT_xi", StaggerLocation::EdgeXi, 9, 0});
    register_field(FieldDescriptor{"pinvAT_xi", StaggerLocation::EdgeXi, 9, 0});
    register_field(FieldDescriptor{"pinvGT_eta", StaggerLocation::EdgeEt, 9, 0});
    register_field(FieldDescriptor{"pinvAT_eta", StaggerLocation::EdgeEt, 9, 0});
    register_field(FieldDescriptor{"pinvGT_zeta", StaggerLocation::EdgeZe, 9, 0});
    register_field(FieldDescriptor{"pinvAT_zeta", StaggerLocation::EdgeZe, 9, 0});

    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    auto cross = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]};
    };
    auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto plus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };

    {
        // Calc_ JDxi  = Jac\nabla\xi = Area
        auto &data_ = field("JDxi");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                        ori_yp = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                        ori_xyp = {x(i, j + 1, k + 1), y(i, j + 1, k + 1), z(i, j + 1, k + 1)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ JDet  = Jac\nabla\eta = Area
        auto &data_ = field("JDet");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                        ori_yp = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                        ori_xyp = {x(i + 1, j, k + 1), y(i + 1, j, k + 1), z(i + 1, j, k + 1)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ JDze  = Jac\nabla\zeta = Area
        auto &data_ = field("JDze");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            std::array<double, 3> ori, ori_xp, ori_yp, ori_xyp, rx, ry, rx_, ry_, Area;
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        ori = {x(i, j, k), y(i, j, k), z(i, j, k)};
                        ori_xp = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                        ori_yp = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                        ori_xyp = {x(i + 1, j + 1, k), y(i + 1, j + 1, k), z(i + 1, j + 1, k)};
                        rx = minus(ori_xp, ori);
                        ry = minus(ori_yp, ori);
                        rx_ = minus(ori_yp, ori_xyp);
                        ry_ = minus(ori_xp, ori_xyp);
                        Area = plus(cross(rx, ry), cross(rx_, ry_));
                        data(i, j, k, 0) = 0.5 * Area[0];
                        data(i, j, k, 1) = 0.5 * Area[1];
                        data(i, j, k, 2) = 0.5 * Area[2];
                    }
        }
    }

    {
        // Calc_ Jac = 1 / 3* \SUM Area \cdot dr
        auto &Jac_ = field("Jac");
        auto &Axi_ = field("JDxi");
        auto &Aeta_ = field("JDet");
        auto &Azeta_ = field("JDze");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &Jac = Jac_[ib];
            auto &A_xi = Axi_[ib];
            auto &A_et = Aeta_[ib];
            auto &A_ze = Azeta_[ib];

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto &cx = grd->grids(ib).dual_x; // cell center
            auto &cy = grd->grids(ib).dual_y;
            auto &cz = grd->grids(ib).dual_z;

            Int3 lo = Jac.inner_lo();
            Int3 hi = Jac.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        std::array<double, 3> xc = {cx(i + 1, j + 1, k + 1),
                                                    cy(i + 1, j + 1, k + 1),
                                                    cz(i + 1, j + 1, k + 1)};

                        double V = 0.0;

                        auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
                        {
                            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                        };
                        auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
                        {
                            return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
                        };

                        // 下面要小心：对 minus 面，把 A 反号变成“对 cell 的外法向”
                        // 并且 face center 要用对应那 4 个顶点的平均

                        // --- xi- face at i ---
                        std::array<double, 3> Af_xm = {-A_xi(i, j, k, 0),
                                                       -A_xi(i, j, k, 1),
                                                       -A_xi(i, j, k, 2)};
                        std::array<double, 3> xf_xm = {
                            0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
                            0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
                            0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};
                        V += dot(Af_xm, minus(xf_xm, xc));

                        // --- xi+ face at i+1 ---
                        std::array<double, 3> Af_xp = {A_xi(i + 1, j, k, 0),
                                                       A_xi(i + 1, j, k, 1),
                                                       A_xi(i + 1, j, k, 2)};
                        std::array<double, 3> xf_xp = {
                            0.25 * (x(i + 1, j, k) + x(i + 1, j + 1, k) + x(i + 1, j, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i + 1, j, k) + y(i + 1, j + 1, k) + y(i + 1, j, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i + 1, j, k) + z(i + 1, j + 1, k) + z(i + 1, j, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_xp, minus(xf_xp, xc));

                        // --- eta- face at j ---
                        std::array<double, 3> Af_em = {-A_et(i, j, k, 0),
                                                       -A_et(i, j, k, 1),
                                                       -A_et(i, j, k, 2)};
                        std::array<double, 3> xf_em = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};
                        V += dot(Af_em, minus(xf_em, xc));

                        // --- eta+ face at j+1 ---
                        std::array<double, 3> Af_ep = {A_et(i, j + 1, k, 0),
                                                       A_et(i, j + 1, k, 1),
                                                       A_et(i, j + 1, k, 2)};
                        std::array<double, 3> xf_ep = {
                            0.25 * (x(i, j + 1, k) + x(i + 1, j + 1, k) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j + 1, k) + y(i + 1, j + 1, k) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j + 1, k) + z(i + 1, j + 1, k) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_ep, minus(xf_ep, xc));

                        // --- zeta- face at k ---
                        std::array<double, 3> Af_zm = {-A_ze(i, j, k, 0),
                                                       -A_ze(i, j, k, 1),
                                                       -A_ze(i, j, k, 2)};
                        std::array<double, 3> xf_zm = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};
                        V += dot(Af_zm, minus(xf_zm, xc));

                        // --- zeta+ face at k+1 ---
                        std::array<double, 3> Af_zp = {A_ze(i, j, k + 1, 0),
                                                       A_ze(i, j, k + 1, 1),
                                                       A_ze(i, j, k + 1, 2)};
                        std::array<double, 3> xf_zp = {
                            0.25 * (x(i, j, k + 1) + x(i + 1, j, k + 1) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j, k + 1) + y(i + 1, j, k + 1) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};
                        V += dot(Af_zp, minus(xf_zp, xc));

                        Jac(i, j, k, 0) = V / 3.0;
                    }
        }
    }

    {
        // Calc_ covariant basis vectors a_xi, a_eta, a_zeta at cell centers
        auto &axi_ = field("a_xi");
        auto &aeta_ = field("a_eta");
        auto &azeta_ = field("a_zeta");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &axi = axi_[ib];
            auto &aeta = aeta_[ib];
            auto &azeta = azeta_[ib];

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 lo = axi.inner_lo();
            Int3 hi = axi.inner_hi();

            auto node = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {x(i, j, k), y(i, j, k), z(i, j, k)};
            };

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // 8 corner nodes of cell (i,j,k)
                        auto r000 = node(i, j, k);
                        auto r100 = node(i + 1, j, k);
                        auto r010 = node(i, j + 1, k);
                        auto r110 = node(i + 1, j + 1, k);

                        auto r001 = node(i, j, k + 1);
                        auto r101 = node(i + 1, j, k + 1);
                        auto r011 = node(i, j + 1, k + 1);
                        auto r111 = node(i + 1, j + 1, k + 1);

                        // ---- a_xi: average of 4 xi-directed edges ----
                        auto ex0 = minus(r100, r000);
                        auto ex1 = minus(r110, r010);
                        auto ex2 = minus(r101, r001);
                        auto ex3 = minus(r111, r011);
                        auto exs = plus(plus(ex0, ex1), plus(ex2, ex3));

                        axi(i, j, k, 0) = 0.25 * exs[0];
                        axi(i, j, k, 1) = 0.25 * exs[1];
                        axi(i, j, k, 2) = 0.25 * exs[2];

                        // ---- a_eta: average of 4 eta-directed edges ----
                        auto ey0 = minus(r010, r000);
                        auto ey1 = minus(r110, r100);
                        auto ey2 = minus(r011, r001);
                        auto ey3 = minus(r111, r101);
                        auto eys = plus(plus(ey0, ey1), plus(ey2, ey3));

                        aeta(i, j, k, 0) = 0.25 * eys[0];
                        aeta(i, j, k, 1) = 0.25 * eys[1];
                        aeta(i, j, k, 2) = 0.25 * eys[2];

                        // ---- a_zeta: average of 4 zeta-directed edges ----
                        auto ez0 = minus(r001, r000);
                        auto ez1 = minus(r101, r100);
                        auto ez2 = minus(r011, r010);
                        auto ez3 = minus(r111, r110);
                        auto ezs = plus(plus(ez0, ez1), plus(ez2, ez3));

                        azeta(i, j, k, 0) = 0.25 * ezs[0];
                        azeta(i, j, k, 1) = 0.25 * ezs[1];
                        azeta(i, j, k, 2) = 0.25 * ezs[2];
                    }
        }
    }

    {
        // ===== Build edge cache: pinvGT_{xi,eta,ze} and pinvAT_{xi,eta,ze} =====
        auto &pinvGT_xi_ = field("pinvGT_xi");
        auto &pinvAT_xi_ = field("pinvAT_xi");
        auto &pinvGT_et_ = field("pinvGT_eta");
        auto &pinvAT_et_ = field("pinvAT_eta");
        auto &pinvGT_ze_ = field("pinvGT_zeta");
        auto &pinvAT_ze_ = field("pinvAT_zeta");

        auto &Axi_ = field("JDxi");
        auto &Aeta_ = field("JDet");
        auto &Azeta_ = field("JDze");

        auto scale = [&](const std::array<double, 3> &a, double s) -> std::array<double, 3>
        {
            return {s * a[0], s * a[1], s * a[2]};
        };

        auto getvec3 = [&](FieldBlock &fb, int i, int j, int k) -> std::array<double, 3>
        {
            return {fb(i, j, k, 0), fb(i, j, k, 1), fb(i, j, k, 2)};
        };

        auto in_range = [&](const Int3 &lo, const Int3 &hi, int i, int j, int k) -> bool
        {
            return (i >= lo.i && i < hi.i &&
                    j >= lo.j && j < hi.j &&
                    k >= lo.k && k < hi.k);
        };

        // invert 3x3 matrix (row-major). return false if singular.
        auto inv3x3 = [&](const double M[9], double invM[9]) -> bool
        {
            const double a00 = M[0], a01 = M[1], a02 = M[2];
            const double a10 = M[3], a11 = M[4], a12 = M[5];
            const double a20 = M[6], a21 = M[7], a22 = M[8];

            const double c00 = a11 * a22 - a12 * a21;
            const double c01 = a02 * a21 - a01 * a22;
            const double c02 = a01 * a12 - a02 * a11;

            const double c10 = a12 * a20 - a10 * a22;
            const double c11 = a00 * a22 - a02 * a20;
            const double c12 = a02 * a10 - a00 * a12;

            const double c20 = a10 * a21 - a11 * a20;
            const double c21 = a01 * a20 - a00 * a21;
            const double c22 = a00 * a11 - a01 * a10;

            const double det = a00 * c00 + a01 * c10 + a02 * c20;
            if (std::fabs(det) < 1e-300)
                return false;

            const double invdet = 1.0 / det;

            invM[0] = c00 * invdet;
            invM[1] = c01 * invdet;
            invM[2] = c02 * invdet;
            invM[3] = c10 * invdet;
            invM[4] = c11 * invdet;
            invM[5] = c12 * invdet;
            invM[6] = c20 * invdet;
            invM[7] = c21 * invdet;
            invM[8] = c22 * invdet;
            return true;
        };

        // pinv(X^T) with Tikhonov: (X X^T + reg I)^-1 X
        auto build_pinv_transpose = [&](const double X[9], double pinv[9])
        {
            // S = X X^T (3x3)
            double S[9] = {0.0};
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += X[r * 3 + kk] * X[c * 3 + kk];
                    S[r * 3 + c] = s;
                }

            const double tr = S[0] + S[4] + S[8];
            double reg = 1e-14 * (tr + 1e-300); // scale-aware regularization

            double Sinv[9];
            bool ok = false;

            for (int it = 0; it < 8; ++it)
            {
                double Stmp[9] = {
                    S[0] + reg, S[1], S[2],
                    S[3], S[4] + reg, S[5],
                    S[6], S[7], S[8] + reg};

                if (inv3x3(Stmp, Sinv))
                {
                    ok = true;
                    break;
                }
                reg *= 10.0;
            }

            if (!ok)
            {
                for (int i = 0; i < 9; ++i)
                    pinv[i] = 0.0;
                return;
            }

            // pinv = Sinv * X
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += Sinv[r * 3 + kk] * X[kk * 3 + c];
                    pinv[r * 3 + c] = s;
                }
        };

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto node = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {x(i, j, k), y(i, j, k), z(i, j, k)};
            };

            auto &A_xi = Axi_[ib];
            auto &A_et = Aeta_[ib];
            auto &A_ze = Azeta_[ib];

            const Int3 loAx = A_xi.inner_lo(), hiAx = A_xi.inner_hi();
            const Int3 loAe = A_et.inner_lo(), hiAe = A_et.inner_hi();
            const Int3 loAz = A_ze.inner_lo(), hiAz = A_ze.inner_hi();

            auto d_node_xi = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii + 1, jj, kk), node(ii - 1, jj, kk)), 0.5);
            };

            auto d_node_eta = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii, jj + 1, kk), node(ii, jj - 1, kk)), 0.5);
            };

            auto d_node_ze = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii, jj, kk + 1), node(ii, jj, kk - 1)), 0.5);
            };

            // -------------------------
            // EdgeXi: pinvGT_xi / pinvAT_xi
            // -------------------------
            {
                auto &pinvGT = pinvGT_xi_[ib];
                auto &pinvAT = pinvAT_xi_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need j-1, k-1 for symmetric stencils
                            // bool okA =
                            //     in_range(loAe, hiAe, i, j, k) && in_range(loAe, hiAe, i, j, k - 1) &&
                            //     in_range(loAz, hiAz, i, j, k) && in_range(loAz, hiAz, i, j - 1, k);

                            // // 8-point xi-face average around edge
                            // for (int di : {0, 1})
                            //     for (int dj : {0, -1})
                            //         for (int dk : {0, -1})
                            //             okA = okA && in_range(loAx, hiAx, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            // G columns: g_xi, g_eta, g_zeta (all at the xi-edge center)
                            auto r000 = node(i, j, k);
                            auto r100 = node(i + 1, j, k);

                            auto g_xi = minus(r100, r000);
                            auto g_eta = scale(plus(d_node_eta(i, j, k), d_node_eta(i + 1, j, k)), 0.5); // avg endpoints
                            auto g_ze = scale(plus(d_node_ze(i, j, k), d_node_ze(i + 1, j, k)), 0.5);    // avg endpoints

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns: A_xi, A_eta, A_zeta (all co-located at the xi-edge center)
                            auto A_eta = scale(plus(getvec3(A_et, i, j, k), getvec3(A_et, i, j, k - 1)), 0.5);
                            auto A_zeta = scale(plus(getvec3(A_ze, i, j, k), getvec3(A_ze, i, j - 1, k)), 0.5);

                            std::array<double, 3> A_xi_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, 1})
                                for (int dj : {0, -1})
                                    for (int dk : {0, -1})
                                        A_xi_edge = plus(A_xi_edge, getvec3(A_xi, i + di, j + dj, k + dk));
                            A_xi_edge = scale(A_xi_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta[0], A_zeta[0],
                                A_xi_edge[1], A_eta[1], A_zeta[1],
                                A_xi_edge[2], A_eta[2], A_zeta[2]};

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }

            // -------------------------
            // EdgeEt: pinvGT_eta / pinvAT_eta
            // -------------------------
            {
                auto &pinvGT = pinvGT_et_[ib];
                auto &pinvAT = pinvAT_et_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need i-1, k-1, and j+1 (edge direction)
                            // bool okA =
                            //     in_range(loAx, hiAx, i, j, k) && in_range(loAx, hiAx, i, j, k - 1) &&
                            //     in_range(loAz, hiAz, i, j, k) && in_range(loAz, hiAz, i - 1, j, k);

                            // // 8-point eta-face average around edge
                            // for (int di : {0, -1})
                            //     for (int dj : {0, 1})
                            //         for (int dk : {0, -1})
                            //             okA = okA && in_range(loAe, hiAe, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            auto r000 = node(i, j, k);
                            auto r010 = node(i, j + 1, k);

                            auto g_eta = minus(r010, r000);
                            auto g_xi = scale(plus(d_node_xi(i, j, k), d_node_xi(i, j + 1, k)), 0.5);
                            auto g_ze = scale(plus(d_node_ze(i, j, k), d_node_ze(i, j + 1, k)), 0.5); // avg endpoints

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns co-located at eta-edge center
                            auto A_xi_edge = scale(plus(getvec3(A_xi, i, j, k), getvec3(A_xi, i, j, k - 1)), 0.5);
                            auto A_zeta = scale(plus(getvec3(A_ze, i, j, k), getvec3(A_ze, i - 1, j, k)), 0.5);

                            std::array<double, 3> A_eta_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, -1})
                                for (int dj : {0, 1})
                                    for (int dk : {0, -1})
                                        A_eta_edge = plus(A_eta_edge, getvec3(A_et, i + di, j + dj, k + dk));
                            A_eta_edge = scale(A_eta_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta_edge[0], A_zeta[0],
                                A_xi_edge[1], A_eta_edge[1], A_zeta[1],
                                A_xi_edge[2], A_eta_edge[2], A_zeta[2]};

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }
            // -------------------------
            // EdgeZe: pinvGT_ze / pinvAT_ze
            // -------------------------
            {
                auto &pinvGT = pinvGT_ze_[ib];
                auto &pinvAT = pinvAT_ze_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need i-1, j-1, and k+1 (edge direction)
                            // bool okA =
                            //     in_range(loAx, hiAx, i, j, k) && in_range(loAx, hiAx, i, j - 1, k) &&
                            //     in_range(loAe, hiAe, i, j, k) && in_range(loAe, hiAe, i - 1, j, k);

                            // // 8-point zeta-face average around edge
                            // for (int di : {0, -1})
                            //     for (int dj : {0, -1})
                            //         for (int dk : {0, 1})
                            //             okA = okA && in_range(loAz, hiAz, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            auto r000 = node(i, j, k);
                            auto r001 = node(i, j, k + 1);

                            auto g_ze = minus(r001, r000);
                            auto g_xi = scale(plus(d_node_xi(i, j, k), d_node_xi(i, j, k + 1)), 0.5);
                            auto g_eta = scale(plus(d_node_eta(i, j, k), d_node_eta(i, j, k + 1)), 0.5);

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns co-located at zeta-edge center
                            auto A_xi_edge = scale(plus(getvec3(A_xi, i, j, k), getvec3(A_xi, i, j - 1, k)), 0.5);
                            auto A_eta = scale(plus(getvec3(A_et, i, j, k), getvec3(A_et, i - 1, j, k)), 0.5);

                            std::array<double, 3> A_ze_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, -1})
                                for (int dj : {0, -1})
                                    for (int dk : {0, 1})
                                        A_ze_edge = plus(A_ze_edge, getvec3(A_ze, i + di, j + dj, k + dk));
                            A_ze_edge = scale(A_ze_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta[0], A_ze_edge[0],
                                A_xi_edge[1], A_eta[1], A_ze_edge[1],
                                A_xi_edge[2], A_eta[2], A_ze_edge[2]};

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }
        }
    }
    // GCL Test
    // {
    //     auto &Jac_ = field("Jac");
    //     auto &Axi_ = field("JDxi");
    //     auto &Aeta_ = field("JDet");
    //     auto &Azeta_ = field("JDze");
    //     for (int ib = 0; ib < grd->nblock; ++ib)
    //     {
    //         auto &Jac = Jac_[ib];
    //         auto &A_xi = Axi_[ib];
    //         auto &A_et = Aeta_[ib];
    //         auto &A_ze = Azeta_[ib];
    //         Int3 lo = Jac.inner_lo();
    //         Int3 hi = Jac.inner_hi();
    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     double error0 = A_xi(i + 1, j, k, 0) - A_xi(i, j, k, 0);
    //                     double error1 = A_xi(i + 1, j, k, 1) - A_xi(i, j, k, 1);
    //                     double error2 = A_xi(i + 1, j, k, 2) - A_xi(i, j, k, 2);
    //                     error0 += A_et(i, j + 1, k, 0) - A_et(i, j, k, 0);
    //                     error1 += A_et(i, j + 1, k, 1) - A_et(i, j, k, 1);
    //                     error2 += A_et(i, j + 1, k, 2) - A_et(i, j, k, 2);
    //                     error0 += A_ze(i, j, k + 1, 0) - A_ze(i, j, k, 0);
    //                     error1 += A_ze(i, j, k + 1, 1) - A_ze(i, j, k, 1);
    //                     error2 += A_ze(i, j, k + 1, 2) - A_ze(i, j, k, 2);
    //                     std::cout << sqrt(error0 * error0 + error1 * error1 + error2 * error2) << "\n";
    //                 }
    //     }
    // }
}