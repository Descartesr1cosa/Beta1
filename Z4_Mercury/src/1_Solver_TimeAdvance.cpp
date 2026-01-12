
#include "MercurySolver.h"

void MercurySolver::Time_Advance()
{
    // 1) dq/db 清零（对应 Fortran dq=0, db=0）
    ZeroRHS_();

    // 2) tvdrhs_ns：流体 RHS（先做 H、Na 各一次）
    Scheme_U_();

    // 3) RHS_B：磁场 RHS（使用 B_cell + U_plus，写进 dB）
    Scheme_B_();

    // 4) source
    AddSourceToRHS_B();     // 对应 source_species 中对 db 的三项补充
    AddSourceToRHS_Fluid(); // 对 RHS_H/RHS_Na 加 source_species 的 dq 部分

    // 5) 人工黏性 cf/cb（Fortran 在更新前加到 dq/db）
    // add_artificial_viscosity_();

    // 6) 显式 Euler 更新（对应 q+=dt*dq, qb+=dt*db）
    ApplyUpdate_Euler_();

    // 7) 低密度/负压修复（尽量按 Fortran：邻域平均 + 重建 E）
    // RepairNonPhysical_();
}

void MercurySolver::ZeroRHS_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB = fld_->field(fid_.fid_RHS_b, ib);

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        if (RHSH.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSH(i, j, k, m) = 0.0;

        if (RHSN.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSN(i, j, k, m) = 0.0;

        if (RHSB.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 3; ++m)
                            RHSB(i, j, k, m) = 0.0;
    }
}

void MercurySolver::ApplyUpdate_Euler_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Ub = fld_->field(fid_.fid_U_b, ib);

        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB = fld_->field(fid_.fid_RHS_b, ib);

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        if (RHSH.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UH(i, j, k, m) += dt * RHSH(i, j, k, m);

        if (RHSN.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UN(i, j, k, m) += dt * RHSN(i, j, k, m);

        if (RHSB.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 3; ++m)
                            Ub(i, j, k, m) += dt * RHSB(i, j, k, m);
    }
}