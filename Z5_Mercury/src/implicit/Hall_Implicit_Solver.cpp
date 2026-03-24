#include "4_Hall_Implicit.h"

#ifdef HALL_IMPLICIT

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iomanip>

#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"
#include "1_Boundary.h"

ImplicitHallSolver::~ImplicitHallSolver()
{
    FinalizePetsc();
}

void ImplicitHallSolver::Setup(Grid *grd,
                               TOPO::Topology *topo,
                               Field *fld,
                               Halo *halo,
                               Param *par,
                               MercuryBoundary *bound,
                               const SolverFields &fid,
                               const TOPO::TopologyEquiv &equiv,
                               const HALO_OWNER::EdgeOwnerSyncPattern &owner_pat)
{
    grd_ = grd;
    topo_ = topo;
    fld_ = fld;
    halo_ = halo;
    par_ = par;
    bound_ = bound;

    fid_ = fid;
    equiv_ = equiv;
    owner_pat_ = owner_pat;

    if (!fld_ || !bound_)
        throw std::runtime_error("ImplicitHallSolver::Setup: null core pointer.");

    x_local_.resize(static_cast<size_t>(equiv_.n_local_edge_owner), 0.0);
    eh_pred_local_.resize(static_cast<size_t>(equiv_.n_local_edge_owner), 0.0);

    HALO_OWNER::gather_local_owner_edges_sorted(equiv_, owner_edges_sorted_);
    if (static_cast<int>(owner_edges_sorted_.size()) != equiv_.n_local_edge_owner)
    {
        throw std::runtime_error(
            "ImplicitHallSolver::Setup: owner_edges_sorted_ size mismatch.");
    }
}

void ImplicitHallSolver::CheckReady_() const
{
    if (!fld_ || !bound_)
        throw std::runtime_error("ImplicitHallSolver not setup.");
    if (!cb_.sync_Bface || !cb_.sync_Ehalledge || !cb_.calc_PV ||
        !cb_.calc_Uplus || !cb_.build_Ehall_from_current_B)
        throw std::runtime_error("ImplicitHallSolver callbacks are not fully bound.");
}

void ImplicitHallSolver::InitializePetsc()
{
    PetscBool is_init = PETSC_FALSE;
    PetscInitialized(&is_init);
    if (!is_init)
    {
        throw std::runtime_error("PETSc not initialized before ImplicitHallSolver::InitializePetsc()");
    }

    if (petsc_ready_)
        return;
    CheckReady_();
    CreatePetscObjects_();
    petsc_ready_ = true;
}

void ImplicitHallSolver::FinalizePetsc()
{
    DestroyPetscObjects_();
    petsc_ready_ = false;
}

void ImplicitHallSolver::CreatePetscObjects_()
{
    const PetscInt nloc = static_cast<PetscInt>(equiv_.n_local_edge_owner);
    const PetscInt nglb = static_cast<PetscInt>(equiv_.n_global_edge_owner);

    VecCreateMPI(PETSC_COMM_WORLD, nloc, nglb, &X_);
    VecDuplicate(X_, &F_);

    SNESCreate(PETSC_COMM_WORLD, &snes_);
    SNESSetFunction(snes_, F_, &ImplicitHallSolver::FormFunction_, this);

    // matrix-free Jacobian
    MatCreateSNESMF(snes_, &Jmf_);
    SNESSetJacobian(snes_, Jmf_, Jmf_, MatMFFDComputeJacobian, nullptr);

    SNESSetType(snes_, SNESNEWTONLS);

    KSP ksp = nullptr;
    SNESGetKSP(snes_, &ksp);
    KSPSetType(ksp, KSPGMRES);

    PC pc = nullptr;
    KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCJACOBI);
    PCSetType(pc, PCNONE);

    SNESSetFromOptions(snes_);
}

void ImplicitHallSolver::DestroyPetscObjects_()
{
    if (Jmf_)
    {
        MatDestroy(&Jmf_);
        Jmf_ = nullptr;
    }
    if (F_)
    {
        VecDestroy(&F_);
        F_ = nullptr;
    }
    if (X_)
    {
        VecDestroy(&X_);
        X_ = nullptr;
    }
    if (snes_)
    {
        SNESDestroy(&snes_);
        snes_ = nullptr;
    }
}

void ImplicitHallSolver::SolveOneStep(double dt, bool if_outres)
{
    CheckReady_();
    if (!petsc_ready_)
        InitializePetsc();

    dt_ = dt;

    // 1) snapshot B*
    SnapshotCurrentBface_();

    // 2) 初值：拿当前显式/上一步的 Ehall 当 guess
    HALO_OWNER::pack_owner_edge_1form_local(
        *fld_, fid_.fid_Ehall, equiv_, owner_edges_sorted_, x_local_);

    PetscScalar *xarr = nullptr;
    VecGetArray(X_, &xarr);
    for (PetscInt i = 0; i < static_cast<PetscInt>(x_local_.size()); ++i)
        xarr[i] = x_local_[static_cast<size_t>(i)];
    VecRestoreArray(X_, &xarr);

    // 3) solve
    SNESSolve(snes_, nullptr, X_);

    // 4) 把最终解写回 Ehall
    UnpackVecToEhallField_(X_);

    // 5) 用最终 Ehall 做整步 CT 更新：
    //    B^{n+1} = B* + dt * RHS(Ehall), RHS=-curl(Ehall)
    RestoreCurrentBfaceFromSnapshot_();

    cb_.sync_Ehalledge();

    ClearFaceTriplet_(fid_.fid_RHS_b);

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);
        auto &Rxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &Ret = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &Rze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        if (!Exi.is_allocated())
            continue;
        CTOperators::CurlEdgeToFace(ib, Exi, Eet, Eze, Rxi, Ret, Rze, -1.0);
    }

    if (if_outres)
    {
        SNESConvergedReason reason;
        PetscInt its;
        PetscReal fnorm = 0.0;

        SNESGetConvergedReason(snes_, &reason);
        SNESGetIterationNumber(snes_, &its);
        SNESGetFunctionNorm(snes_, &fnorm);

        auto max_abs_face_delta = [&](int fid,
                                      const std::vector<std::vector<double>> &snap) -> double
        {
            double local_max = 0.0;

            for (int ib = 0; ib < fld_->num_blocks(); ++ib)
            {
                auto &F = fld_->field(fid, ib);
                if (!F.is_allocated())
                    continue;

                Int3 lo = F.inner_lo(), hi = F.inner_hi();
                const auto &buf = snap[static_cast<size_t>(ib)];

                size_t t = 0;
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k, ++t)
                            local_max = std::max(local_max, std::abs(F(i, j, k, 0) - buf[t]));
            }

            double global_max = 0.0;
            MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return global_max;
        };

        auto max_dB_hall = [&]() -> double
        {
            const double dxi = max_abs_face_delta(fid_.fid_B.xi, Bstar_xi_);
            const double det = max_abs_face_delta(fid_.fid_B.eta, Bstar_eta_);
            const double dze = max_abs_face_delta(fid_.fid_B.zeta, Bstar_ze_);
            return std::max(dxi, std::max(det, dze));
        };

        const double maxEhall = MaxAbsTriplet_(fid_.fid_Ehall);
        const double maxB = MaxAbsTriplet_(fid_.fid_B);
        const double maxRHSB = MaxAbsTriplet_(fid_.fid_RHS_b);
        const double maxdBhall_est = dt_ * maxRHSB;
        int rank = par_->GetInt("myid");
        if (rank == 0)
        {
            std::cout << std::scientific << std::setprecision(6)
                      << "[HallImplicit] reason=" << reason
                      << "  its=" << its
                      << "  |F|_2=" << fnorm
                      << "  max|Ehall|=" << maxEhall
                      << "  dt*max|RHS_B|=" << maxdBhall_est
                      << std::endl;
        }
    }

    AddFaceInnerFromRHS_(fid_.fid_B.xi, fid_.fid_RHS_b.xi, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.eta, fid_.fid_RHS_b.eta, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.zeta, fid_.fid_RHS_b.zeta, dt_);

    cb_.sync_Bface();
}

double ImplicitHallSolver::MaxAbsTriplet_(const IdTriplet &fid_triplet)
{
    double local_max = 0.0;

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        for (int fid : {fid_triplet.xi, fid_triplet.eta, fid_triplet.zeta})
        {
            auto &F = fld_->field(fid, ib);
            if (!F.is_allocated())
                continue;

            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        local_max = std::max(local_max, std::abs(F(i, j, k, 0)));
        }
    }

    double global_max = 0.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_max;
}

#endif