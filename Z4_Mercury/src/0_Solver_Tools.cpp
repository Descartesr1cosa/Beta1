#include <cmath>
#include <algorithm>
#include <iostream>

#include "MercurySolver.h"

void MercurySolver::calc_physical_constant(Param *par)
{
    // -------- constants / reference ----------
    auto cst = par_->GetDou_List("constant");
    auto ref = par_->GetDou_List("REF");

    gamma_ = cst.data["gamma"];
    double R_uni = cst.data["R_uni"];

    const double U_ref = ref.data.at("U");
    const double T_ref = ref.data.at("T");
    M_H = par->GetDou("mole_mass1");
    M_Na = par->GetDou("mole_mass2");

    // 无量纲状态方程系数：p = rho * T * coeff
    // coeff = (R_uni * T_ref) / (M * U_ref^2)
    state_coeff_H = (R_uni * T_ref) / (M_H * U_ref * U_ref);
    state_coeff_Na = (R_uni * T_ref) / (M_Na * U_ref * U_ref);

    CFL = par_->GetDou("CFL");
}

void MercurySolver::calc_Bcell()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Bcell = fld_->field(fid_.fid_Bcell, ib);
        if (!Bcell.is_allocated())
            continue;

        FieldBlock &Badd = fld_->field(fid_.fid_Badd, ib);
        FieldBlock &Ub = fld_->field(fid_.fid_U_b, ib);

        const Int3 lo = Bcell.get_lo();
        const Int3 hi = Bcell.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    for (int m = 0; m < 3; ++m)
                    {
                        double v = Badd(i, j, k, m);
                        if (Ub.is_allocated())
                            v += Ub(i, j, k, m);
                        Bcell(i, j, k, m) = v;
                    }
                }
    }
}

void MercurySolver::calc_PV()
{
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    auto fill_one = [&](int fidU, int fidPV, double coeff)
    {
        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            FieldBlock &U = fld_->field(fidU, ib);
            FieldBlock &PV = fld_->field(fidPV, ib);
            if (!U.is_allocated() || !PV.is_allocated())
                continue;

            const Int3 lo = PV.get_lo(); // 含 ghost：[-ng, ...]
            const Int3 hi = PV.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double rho0 = U(i, j, k, 0);
                        const double rho = std::max(rho0, rho_floor);

                        const double u = U(i, j, k, 1) / rho;
                        const double v = U(i, j, k, 2) / rho;
                        const double w = U(i, j, k, 3) / rho;

                        const double E = U(i, j, k, 4);
                        const double ke = 0.5 * rho * (u * u + v * v + w * w);

                        // 压力：p = (gamma-1) * (E - ke)
                        double eint = E - ke;
                        if (eint < 0.0)
                            eint = 0.0;
                        double p = (gamma_ - 1.0) * eint;
                        p = std::max(p, p_floor);

                        // 温度：T = p / (rho * coeff)
                        // 对 Na 用 coeff_Na -> 等价于 Fortran 的 T = p/(n kB)，其中 n= rho/m_Na
                        double T = p / (rho * coeff);

                        PV(i, j, k, 0) = u;
                        PV(i, j, k, 1) = v;
                        PV(i, j, k, 2) = w;
                        PV(i, j, k, 3) = p;
                        PV(i, j, k, 4) = T;
                    }
        }
    };

    fill_one(fid_.fid_U_H, fid_.fid_PV_H, state_coeff_H);
    fill_one(fid_.fid_U_Na, fid_.fid_PV_Na, state_coeff_Na);
}

void MercurySolver::calc_Uplus()
{
    const double rho_floor = 1e-20;
    const double inv23 = M_H / M_Na; // 与 Fortran sm2≈23*sm1 对齐

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib); // 新增

        if (!UH.is_allocated() || !UN.is_allocated() || !Up.is_allocated())
            continue;

        const Int3 lo = Up.get_lo();
        const Int3 hi = Up.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH = std::max(UH(i, j, k, 0), rho_floor);
                    const double rhoNa = std::max(UN(i, j, k, 0), 0.0);

                    const double uH = UH(i, j, k, 1) / rhoH;
                    const double vH = UH(i, j, k, 2) / rhoH;
                    const double wH = UH(i, j, k, 3) / rhoH;

                    double uNa = 0, vNa = 0, wNa = 0;
                    if (rhoNa > rho_floor)
                    {
                        uNa = UN(i, j, k, 1) / rhoNa;
                        vNa = UN(i, j, k, 2) / rhoNa;
                        wNa = UN(i, j, k, 3) / rhoNa;
                    }

                    const double nH = rhoH;
                    const double nNa = rhoNa * inv23;
                    const double nt = nH + nNa;

                    if (nt <= rho_floor)
                    {
                        Up(i, j, k, 0) = 0.0;
                        Up(i, j, k, 1) = 0.0;
                        Up(i, j, k, 2) = 0.0;
                    }
                    else
                    {
                        Up(i, j, k, 0) = (nH * uH + nNa * uNa) / nt;
                        Up(i, j, k, 1) = (nH * vH + nNa * vNa) / nt;
                        Up(i, j, k, 2) = (nH * wH + nNa * wNa) / nt;
                    }
                }
    }
}
