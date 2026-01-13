// DipoleBCell.h
#pragma once

#include <vector>
#include <string>
#include <cmath>

#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/2_MPCNS_Field.h"
#include "9_Vector_tool.h"

// 直接在 cell center 计算偶极子磁场 B，并写入（或累加到）BCell(3分量)
// 参数约定：
//   - L_ref [m], B_ref [T], constant.mu_mag = mu0 [H/m]
//   - num_of_dipole = N
//   - Dipole{i} list: x0,y0,z0 [m], mx,my,mz [A*m^2], a [m]
//     可选：r_cut [m] （若不提供则默认 0）
//   - B_add_x, B_add_y, B_add_z：均匀背景场（无量纲，已经除以 B_ref）
//     （如果你更喜欢用 nT 或 T，这里也可以改成读取后自行归一化）
class Dipoles
{
public:
    struct Dipole
    {
        Vec3 r0_nd;     // 无量纲位置：r0/L_ref
        Vec3 m_nd;      // 无量纲偶极矩：mu0/(4pi) * m_si / (B_ref * L_ref^3)
        double a_nd;    // 无量纲softening：a/L_ref
        double rcut_nd; // 无量纲截断半径：r_cut/L_ref（可选）
    };

    Dipoles() = default;

    void load_from_param(Param *par)
    {
        dips_.clear();

        const int n = par->GetInt("num_of_dipole");
        dips_.reserve(n);

        const double L_ref = par->GetDou_List("REF").data["L_ref"];
        const double B_ref = par->GetDou_List("REF").data["B_ref"];
        const double mu0 = par->GetDou_List("constant").data["mu_mag"];

        // 无量纲化系数：m_nd = mu0/(4pi) * m_si / (B_ref * L_ref^3)
        const double factor = (mu0 / (4.0 * M_PI)) / (B_ref * L_ref * L_ref * L_ref);

        for (int i = 1; i <= n; ++i)
        {
            const std::string key = "Dipole" + std::to_string(i);
            auto dl = par->GetDou_List(key);

            const double x0 = dl.data["x0"];
            const double y0 = dl.data["y0"];
            const double z0 = dl.data["z0"];
            const double mx = dl.data["mx"];
            const double my = dl.data["my"];
            const double mz = dl.data["mz"];
            const double a = dl.data["a"];

            // r_cut 可选；若你的 GetDou_List 对缺key会抛异常，
            // 请在参数文件里显式给 r_cut=0
            double r_cut = 0.0;
            auto it = dl.data.find("r_cut");
            if (it != dl.data.end())
                r_cut = it->second;

            Dipole d;
            d.r0_nd = Vec3(x0, y0, z0) / L_ref;
            d.a_nd = a / L_ref;
            d.rcut_nd = r_cut / L_ref;

            // m_si (A*m^2) -> m_nd
            d.m_nd = Vec3(mx, my, mz) * factor;

            dips_.push_back(d);
        }
    }

    void Build_Badd(Grid *grd_, Field *fld_, int Badd_id) const
    {
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Badd = fld_->field(Badd_id, ib);

            auto &x = grd_->grids(ib).x; // 节点坐标（无量纲或物理？取决于你的网格scale）
            auto &y = grd_->grids(ib).y;
            auto &z = grd_->grids(ib).z;

            const Int3 &lo = Badd.get_lo();
            const Int3 &hi = Badd.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double bx = 0.0;
                        double by = 0.0;
                        double bz = 0.0;

                        // 1) cell center 坐标：用8个corner节点平均（对曲线网格更稳）
                        // 注意：这要求 x/y/z 节点数组至少能访问到 i+1,j+1,k+1
                        const double xc = 0.125 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k) +
                                                   x(i, j, k + 1) + x(i + 1, j, k + 1) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1));
                        const double yc = 0.125 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k) +
                                                   y(i, j, k + 1) + y(i + 1, j, k + 1) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1));
                        const double zc = 0.125 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k) +
                                                   z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1));

                        const Vec3 loc(xc, yc, zc);

                        // 2) 叠加所有偶极子贡献（无量纲 B/B_ref）
                        for (const auto &d : dips_)
                        {
                            Vec3 Bdip;
                            B_of_softened_dipole(loc, d, Bdip);
                            bx += Bdip[0];
                            by += Bdip[1];
                            bz += Bdip[2];
                        }

                        // 写入 Badd（累加或覆盖取决于 clear_first）
                        Badd(i, j, k, 0) += bx;
                        Badd(i, j, k, 1) += by;
                        Badd(i, j, k, 2) += bz;
                    }
        }
    }

private:
    // 与你原 A=(m×r)/R^3 软化形式严格一致的 B=curl(A)（已把 mu0/4pi 吸收进 m_nd）
    // R^2 = r^2 + a^2
    // B = [ 3 (m·r) r + (2 a^2 - r^2) m ] / R^5
    static inline void B_of_softened_dipole(const Vec3 &loc, const Dipole &d, Vec3 &Bout)
    {
        Bout[0] = Bout[1] = Bout[2] = 0.0;

        const Vec3 r = loc - d.r0_nd;
        const double r2 = r.norm2();

        if (d.rcut_nd > 0.0 && r2 < d.rcut_nd * d.rcut_nd)
            return;

        const double a2 = d.a_nd * d.a_nd;
        const double R2 = r2 + a2;
        const double R = std::sqrt(R2);
        const double invR5 = 1.0 / (R2 * R2 * R); // 1/R^5

        const double mdotr = d.m_nd * r; // dot

        // term = 3 (m·r) r + (2 a^2 - r^2) m
        const Vec3 term = (3.0 * mdotr) * r + (2.0 * a2 - r2) * d.m_nd;

        Bout = term * invR5;
    }

private:
    std::vector<Dipole> dips_;
    Vec3 Buni_nd_{0.0, 0.0, 0.0};
};