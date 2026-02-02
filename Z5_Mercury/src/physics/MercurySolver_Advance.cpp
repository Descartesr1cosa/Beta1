
#include "MercurySolver.h"

void MercurySolver::Time_Advance()
{
    // 1) dq/db set to ZERO
    ZeroRHS_();

    // 2) RHS for Flow
    AssembleRHS_Fluid_(); // = Scheme_U_ + AddSourceToRHS_Fluid()

    // 3) RHS for Magnetic fields：CT (Edge EMF -> Curl -> Face RHS)
    AssembleRHS_Induction_CT_();

    // 4) Euler 1st order Advance
    ApplyUpdate_Euler_(); // U += dt*RHS, B_face += dt*RHS_Bface

    // 5) 低密度/负压修复（尽量按 Fortran：邻域平均 + 重建 E）
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
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_et = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_ze = fld_->field(fid_.fid_RHS_b.zeta, ib);

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

        if (RHSB_xi.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_xi(i, j, k, 0) = 0.0;

        if (RHSB_et.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_et(i, j, k, 0) = 0.0;

        if (RHSB_ze.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_ze(i, j, k, 0) = 0.0;
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
        FieldBlock &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        FieldBlock &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        FieldBlock &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

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

        if (RHSB_xi.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt * RHSB_xi(i, j, k, 0);

        if (RHSB_eta.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt * RHSB_eta(i, j, k, 0);

        if (RHSB_zeta.is_allocated())
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt * RHSB_zeta(i, j, k, 0);
    }
}