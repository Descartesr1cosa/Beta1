#include "MercurySolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
    void clear_triplet(Field *fld, const IdTriplet &fid)
    {
        if (!fld)
            return;

        const int ids[3] = {fid.xi, fid.eta, fid.zeta};
        for (int ib = 0; ib < fld->num_blocks(); ++ib)
        {
            for (int q = 0; q < 3; ++q)
            {
                auto &F = fld->field(ids[q], ib);
                if (!F.is_allocated())
                    continue;

                const Int3 lo = F.get_lo();
                const Int3 hi = F.get_hi();
                const int ncomp = F.descriptor().ncomp;
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < ncomp; ++m)
                                F(i, j, k, m) = 0.0;
            }
        }
    }

    int edge_fid_from_dir(const IdTriplet &fid, int dir)
    {
        if (dir == 1)
            return fid.xi;
        if (dir == 2)
            return fid.eta;
        if (dir == 3)
            return fid.zeta;
        throw std::runtime_error("edge_fid_from_dir: invalid edge direction.");
    }

    void setup_like_face_snapshot(Scalar &buf, FieldBlock &F)
    {
        if (!F.is_allocated())
            return;

        const Int3 lo = F.get_lo();
        const Int3 hi = F.get_hi();
        const int ghost = -lo.i;
        if (lo.i != -ghost || lo.j != -ghost || lo.k != -ghost)
            throw std::runtime_error("implicit resistive Bstar scratch: non-standard ghost indexing.");

        buf.SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, ghost);
    }
} // namespace

void MercurySolver::SetupImplicitResistiveDiffusion_()
{
    if (implicit_resistive_ready_)
        return;

    if (!topo_equiv_ || !edge_owner_pat_)
        throw std::runtime_error("SetupImplicitResistiveDiffusion_: topology equivalence is required.");

    PetscBool is_init = PETSC_FALSE;
    PetscInitialized(&is_init);
    if (!is_init)
        throw std::runtime_error("PETSc must be initialized before implicit resistive diffusion setup.");

    BuildImplicitResistiveEdgeDofMap_();

    const int nb = fld_->num_blocks();
    resist_Bstar_xi_.resize(nb);
    resist_Bstar_eta_.resize(nb);
    resist_Bstar_ze_.resize(nb);
    for (int ib = 0; ib < nb; ++ib)
    {
        setup_like_face_snapshot(resist_Bstar_xi_[ib], fld_->field(fid_.fid_B.xi, ib));
        setup_like_face_snapshot(resist_Bstar_eta_[ib], fld_->field(fid_.fid_B.eta, ib));
        setup_like_face_snapshot(resist_Bstar_ze_[ib], fld_->field(fid_.fid_B.zeta, ib));
    }

    const PetscInt nloc = static_cast<PetscInt>(implicit_resistive_dofs_.size()); // 本进程需要隐式求解的DOF
    PetscInt nglb = 0;
    MPI_Allreduce(&nloc, &nglb, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); // 获取全局需要隐式求解的DOF

    if (nglb > 0)
    {
        PetscCallAbort(PETSC_COMM_WORLD, VecCreateMPI(PETSC_COMM_WORLD, nloc, nglb, &implicit_resistive_x_));
        PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(implicit_resistive_x_, &implicit_resistive_b_));
        // build unknown implicit_resistive_x_ implicit_resistive_b_

        PetscCallAbort(PETSC_COMM_WORLD,
                       MatCreateShell(PETSC_COMM_WORLD, nloc, nloc, nglb, nglb,
                                      this, &implicit_resistive_A_));
        PetscCallAbort(PETSC_COMM_WORLD,
                       MatShellSetOperation(implicit_resistive_A_, MATOP_MULT,
                                            (void (*)(void))&MercurySolver::MatMultImplicitResistive_));
        // 矩阵向量乘法 A \cdot x, 采用无矩阵形式

        PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &implicit_resistive_ksp_)); // 创建 PETSc Krylov 求解器对象
        PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(implicit_resistive_ksp_,
                                                         implicit_resistive_A_,
                                                         implicit_resistive_A_));        // 告诉 KSP 系数矩阵，预处理矩阵
        PetscCallAbort(PETSC_COMM_WORLD, KSPSetType(implicit_resistive_ksp_, KSPGMRES)); // 选择 GMRES 求解器

        PC pc = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(implicit_resistive_ksp_, &pc));
        PetscCallAbort(PETSC_COMM_WORLD, PCSetType(pc, PCNONE)); // 不使用Precondition

        PetscCallAbort(PETSC_COMM_WORLD,
                       KSPSetTolerances(implicit_resistive_ksp_,
                                        resist_control.implicit_ksp_rtol,
                                        resist_control.implicit_ksp_atol,
                                        PETSC_DEFAULT,
                                        resist_control.implicit_ksp_max_it));
        // rtol--相对残差容差 atol--绝对残差容差 dtol--发散容差,这里用PETSC_DEFAULT max_it--最大迭代次数
        PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(implicit_resistive_ksp_)); // 允许运行程序时用 PETSc 命令行参数覆盖设置
    }

    implicit_resistive_ready_ = true;
}

void MercurySolver::DestroyImplicitResistiveDiffusion_()
{
    if (implicit_resistive_b_)
    {
        VecDestroy(&implicit_resistive_b_);
        implicit_resistive_b_ = nullptr;
    }
    if (implicit_resistive_x_)
    {
        VecDestroy(&implicit_resistive_x_);
        implicit_resistive_x_ = nullptr;
    }
    if (implicit_resistive_A_)
    {
        MatDestroy(&implicit_resistive_A_);
        implicit_resistive_A_ = nullptr;
    }
    if (implicit_resistive_ksp_)
    {
        KSPDestroy(&implicit_resistive_ksp_);
        implicit_resistive_ksp_ = nullptr;
    }
    implicit_resistive_ready_ = false;
}

void MercurySolver::BuildImplicitResistiveEdgeDofMap_()
{
    resist_owner_edges_sorted_.clear(); // 本进程的Edge自由度
    implicit_resistive_dofs_.clear();   // Edge自由度+Edge的扩散系数eta，且要求eta>0.0

    HALO_OWNER::gather_local_owner_edges_sorted(*topo_equiv_, resist_owner_edges_sorted_);

    for (const auto &e : resist_owner_edges_sorted_)
    {
        const double eta = ImplicitResistiveEtaAtEdge_(e);
        if (eta <= 0.0)
            continue;

        implicit_resistive_dofs_.push_back({e, eta}); // 构建Edge自由度+Edge的扩散系数eta
    }

    implicit_resistive_local_.assign(implicit_resistive_dofs_.size(), 0.0);
}

double MercurySolver::ImplicitResistiveEtaAtEdge_(const TOPO::EdgeLocalID &e) const
{
    const int fid_edge = edge_fid_from_dir(fid_.fid_Eres, e.dir);
    auto &E = fld_->field(fid_edge, e.gblock);
    if (!E.is_allocated())
        return 0.0;

    const int di = (e.dir == 1) ? 1 : 0;
    const int dj = (e.dir == 2) ? 1 : 0;
    const int dk = (e.dir == 3) ? 1 : 0;

    auto &x = grd_->grids(e.gblock).x;
    auto &y = grd_->grids(e.gblock).y;
    auto &z = grd_->grids(e.gblock).z;

    const double xm = 0.5 * (x(e.i, e.j, e.k) + x(e.i + di, e.j + dj, e.k + dk));
    const double ym = 0.5 * (y(e.i, e.j, e.k) + y(e.i + di, e.j + dj, e.k + dk));
    const double zm = 0.5 * (z(e.i, e.j, e.k) + z(e.i + di, e.j + dj, e.k + dk));
    const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

    constexpr double r_cut_in = 0.84;
    constexpr double r_cut_out = 1.04;
    constexpr double r0 = 0.84;
    constexpr double r1 = 1.04;
    constexpr double w = 0.02;

    if (r <= r_cut_in || r >= r_cut_out)
        return 0.0;

    const double yita0 = 0.5 * (std::tanh((r - r0) / w) - std::tanh((r - r1) / w));
    return inver_Rem * yita0;
}

void MercurySolver::SnapshotImplicitResistiveBstar_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto copy_one = [](FieldBlock &F, Scalar &buf)
        {
            if (!F.is_allocated())
                return;

            const Int3 lo = F.get_lo();
            const Int3 hi = F.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        buf(i, j, k) = F(i, j, k, 0);
        };

        copy_one(fld_->field(fid_.fid_B.xi, ib), resist_Bstar_xi_[ib]);
        copy_one(fld_->field(fid_.fid_B.eta, ib), resist_Bstar_eta_[ib]);
        copy_one(fld_->field(fid_.fid_B.zeta, ib), resist_Bstar_ze_[ib]);
    }
}

void MercurySolver::RestoreImplicitResistiveBstar_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto restore_one = [](FieldBlock &F, Scalar &buf)
        {
            if (!F.is_allocated())
                return;

            const Int3 lo = F.get_lo();
            const Int3 hi = F.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        F(i, j, k, 0) = buf(i, j, k);
        };

        restore_one(fld_->field(fid_.fid_B.xi, ib), resist_Bstar_xi_[ib]);
        restore_one(fld_->field(fid_.fid_B.eta, ib), resist_Bstar_eta_[ib]);
        restore_one(fld_->field(fid_.fid_B.zeta, ib), resist_Bstar_ze_[ib]);
    }
}

void MercurySolver::UnpackVecToImplicitEres_(Vec X)
{
    clear_triplet(fld_, fid_.fid_Eres);

    const PetscScalar *xarr = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, VecGetArrayRead(X, &xarr));

    for (PetscInt lid = 0; lid < static_cast<PetscInt>(implicit_resistive_dofs_.size()); ++lid)
    {
        const auto &dof = implicit_resistive_dofs_[static_cast<size_t>(lid)];
        const int fid_edge = edge_fid_from_dir(fid_.fid_Eres, dof.edge.dir);
        auto &E = fld_->field(fid_edge, dof.edge.gblock);
        E(dof.edge.i, dof.edge.j, dof.edge.k, 0) = static_cast<double>(xarr[lid]);
    }

    PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArrayRead(X, &xarr));

    HALO_OWNER::sync_edge_1form(*fld_, fid_.fid_Eres, *edge_owner_pat_);
    mercury_bound_.Sync("Eres1form");
}

void MercurySolver::CalcImplicitDeltaJEdgeFromDeltaB_()
{
    clear_triplet(fld_, fid_.fid_dJ);

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_dB.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_dB.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_dB.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_dJ.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_dJ.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_dJ.zeta, iblk);

        if (!Bxi.is_allocated())
            continue;

        auto &beta_xi = fld_->field(fid_.Face_beta.xi, iblk);
        auto &beta_eta = fld_->field(fid_.Face_beta.eta, iblk);
        auto &beta_zeta = fld_->field(fid_.Face_beta.zeta, iblk);

        auto &alpha_xi = fld_->field(fid_.Edge_alpha.xi, iblk);
        auto &alpha_eta = fld_->field(fid_.Edge_alpha.eta, iblk);
        auto &alpha_zeta = fld_->field(fid_.Edge_alpha.zeta, iblk);

        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       beta_xi, beta_eta, beta_zeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);

        auto scale_one = [](FieldBlock &J, FieldBlock &alpha)
        {
            const Int3 lo = J.inner_lo();
            const Int3 hi = J.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        J(i, j, k, 0) *= alpha(i, j, k, 0);
        };

        scale_one(Jxi, alpha_xi);
        scale_one(Jeta, alpha_eta);
        scale_one(Jzeta, alpha_zeta);
    }

    mercury_bound_.Sync("dJ");
}

void MercurySolver::PackImplicitJedgeToVec_(const IdTriplet &fid_Jedge,
                                            Vec Y,
                                            bool multiply_eta,
                                            double x_shift,
                                            Vec X)
{
    const PetscScalar *xarr = nullptr;
    PetscScalar *yarr = nullptr;

    if (X)
        PetscCallAbort(PETSC_COMM_WORLD, VecGetArrayRead(X, &xarr));
    PetscCallAbort(PETSC_COMM_WORLD, VecGetArray(Y, &yarr));

    for (PetscInt lid = 0; lid < static_cast<PetscInt>(implicit_resistive_dofs_.size()); ++lid)
    {
        const auto &dof = implicit_resistive_dofs_[static_cast<size_t>(lid)];
        const int fid_edge = edge_fid_from_dir(fid_Jedge, dof.edge.dir);
        auto &J = fld_->field(fid_edge, dof.edge.gblock);
        const double eta = multiply_eta ? dof.eta : 1.0;
        const double base = X ? x_shift * static_cast<double>(xarr[lid]) : 0.0;
        yarr[lid] = base + eta * J(dof.edge.i, dof.edge.j, dof.edge.k, 0);
    }

    PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArray(Y, &yarr));
    if (X)
        PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArrayRead(X, &xarr));
}

PetscErrorCode MercurySolver::MatMultImplicitResistive_(Mat A, Vec X, Vec Y)
{
    void *ctx = nullptr;
    PetscCall(MatShellGetContext(A, &ctx));
    auto *S = static_cast<MercurySolver *>(ctx);
    PetscCheck(S, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL,
               "Null MercurySolver context in implicit resistive MatShell.");

    try
    {
        S->UnpackVecToImplicitEres_(X);

        clear_triplet(S->fld_, S->fid_.fid_dB);
        for (int ib = 0; ib < S->fld_->num_blocks(); ++ib)
        {
            auto &Exi = S->fld_->field(S->fid_.fid_Eres.xi, ib);
            auto &Eet = S->fld_->field(S->fid_.fid_Eres.eta, ib);
            auto &Eze = S->fld_->field(S->fid_.fid_Eres.zeta, ib);

            auto &Bxi = S->fld_->field(S->fid_.fid_dB.xi, ib);
            auto &Bet = S->fld_->field(S->fid_.fid_dB.eta, ib);
            auto &Bze = S->fld_->field(S->fid_.fid_dB.zeta, ib);

            if (!Exi.is_allocated())
                continue;

            CTOperators::CurlEdgeToFace(ib, Exi, Eet, Eze, Bxi, Bet, Bze,
                                        -S->implicit_resistive_dt_);
        }

        S->mercury_bound_.Sync("dB");
        S->CalcImplicitDeltaJEdgeFromDeltaB_();

        S->PackImplicitJedgeToVec_(S->fid_.fid_dJ, Y,
                                   /*multiply_eta=*/true,
                                   /*x_shift=*/1.0,
                                   X);
        PetscCall(VecScale(Y, -1.0));
        PetscCall(VecAXPY(Y, 2.0, X)); // Y = X - eta*Jdelta after previous Y = -(X+eta*Jdelta)+2X
    }
    catch (const std::exception &e)
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Implicit resistive MatMult failed: %s", e.what());
    }

    return 0;
}

void MercurySolver::SolveImplicitResistiveDiffusion_(double dt_step)
{
    if (!resist_control.is_Mercury_resistance)
        return;

    if (!implicit_resistive_ready_)
        SetupImplicitResistiveDiffusion_();

    mercury_bound_.Sync("Bface");
    SnapshotImplicitResistiveBstar_();

    if (implicit_resistive_dofs_.empty())
    {
        clear_triplet(fld_, fid_.fid_Eres);
        return;
    }

    implicit_resistive_dt_ = dt_step;

    Calc_J_Edge();
    PackImplicitJedgeToVec_(fid_.fid_J, implicit_resistive_b_,
                            /*multiply_eta=*/true,
                            /*x_shift=*/0.0,
                            nullptr);

    PetscCallAbort(PETSC_COMM_WORLD, VecSet(implicit_resistive_x_, 0.0));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSetInitialGuessNonzero(implicit_resistive_ksp_, PETSC_FALSE));
    PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(implicit_resistive_ksp_,
                                              implicit_resistive_b_,
                                              implicit_resistive_x_));

    if (control_.if_outres && par_->GetInt("myid") == 0)
    {
        KSPConvergedReason reason;
        PetscInt its = 0;
        PetscReal rnorm = 0.0;
        KSPGetConvergedReason(implicit_resistive_ksp_, &reason);
        KSPGetIterationNumber(implicit_resistive_ksp_, &its);
        KSPGetResidualNorm(implicit_resistive_ksp_, &rnorm);
        std::cout << "[ImplicitMercuryResistive] reason=" << reason
                  << " its=" << its
                  << " rnorm=" << rnorm << std::endl;
    }

    UnpackVecToImplicitEres_(implicit_resistive_x_);
    ApplyImplicitResistiveUpdate_(dt_step);
}

void MercurySolver::ApplyImplicitResistiveUpdate_(double dt_step)
{
    mercury_bound_.Sync("Eres1form");

    clear_triplet(fld_, fid_.fid_RHS_b_res);
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_Eres.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Eres.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Eres.zeta, ib);

        auto &Rxi = fld_->field(fid_.fid_RHS_b_res.xi, ib);
        auto &Ret = fld_->field(fid_.fid_RHS_b_res.eta, ib);
        auto &Rze = fld_->field(fid_.fid_RHS_b_res.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        CTOperators::CurlEdgeToFace(ib, Exi, Eet, Eze, Rxi, Ret, Rze,
                                    /*multiper=*/-1.0);
    }

    RestoreImplicitResistiveBstar_();
    ApplyUpdate_Euler_BfaceOnly_(dt_step, fid_.fid_RHS_b_res);
    mercury_bound_.Sync("Bface");
}
