#pragma once

#include "00_Mercury_Const.h"

#ifdef HALL_IMPLICIT

#include <functional>
#include <vector>
#include <stdexcept>
#include <petscsnes.h>

#include "0_SolverFields.h"
#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "4_halo/1_MPCNS_Halo_EdgeOwner.h"
#include "operators/CTOperators.h"

// forward decl
class Grid;
namespace TOPO
{
    class Topology;
}
class Field;
class Halo;
class Param;
class MercuryBoundary;
class FieldBlock;

class ImplicitHallSolver
{
public:
    struct Callbacks
    {
        // 把当前 Bface 补齐到可用于 stencil / curl(B)
        std::function<void()> sync_Bface;

        // 把当前 Eedge 补齐到可用于 curl(E)
        std::function<void()> sync_Eedge;

        // 从当前 Bface 计算派生量
        std::function<void()> calc_PV;
        std::function<void()> calc_Uplus;

        // 当前 Bface -> Jedge/Jcell -> Ehall
        // 约定：执行后，fid_.fid_Ehall 中存放预测的 Hall edge EMF
        std::function<void()> build_Ehall_from_current_B;
    };

public:
    ImplicitHallSolver() = default;
    ~ImplicitHallSolver();

    void Setup(Grid *grd,
               TOPO::Topology *topo,
               Field *fld,
               Halo *halo,
               Param *par,
               MercuryBoundary *bound,
               const SolverFields &fid,
               const TOPO::TopologyEquiv &equiv,
               const HALO_OWNER::EdgeOwnerSyncPattern &owner_pat);

    void SetCallbacks(const Callbacks &cb) { cb_ = cb; }

    void SetTheta(double theta) { theta_ = theta; } // 1.0: BE, 0.5: midpoint
    void SetVerbose(bool x) { verbose_ = x; }

    void InitializePetsc();
    void FinalizePetsc();

    // 解 Hall 子步；输入/输出都作用在当前 fld_->B_xi/eta/zeta 上
    void SolveOneStep(double dt);

private:
    Grid *grd_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Field *fld_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};
    MercuryBoundary *bound_{nullptr};

    SolverFields fid_;
    TOPO::TopologyEquiv equiv_;
    HALO_OWNER::EdgeOwnerSyncPattern owner_pat_;
    Callbacks cb_;

    bool petsc_ready_{false};
    bool verbose_{true};
    double dt_{0.0};
    double theta_{0.5};

    // PETSc
    SNES snes_{nullptr};
    Vec X_{nullptr}, F_{nullptr};
    Mat Jmf_{nullptr};

    // local owner-edge buffers
    std::vector<double> x_local_;
    std::vector<double> eh_pred_local_;

    // B* snapshot + RHS scratch
    std::vector<std::vector<double>> Bstar_xi_;
    std::vector<std::vector<double>> Bstar_eta_;
    std::vector<std::vector<double>> Bstar_ze_;

    std::vector<TOPO::EdgeLocalID> owner_edges_sorted_;

private:
    static PetscErrorCode FormFunction_(SNES snes, Vec X, Vec F, void *ctx);

    void CheckReady_() const;
    void CreatePetscObjects_();
    void DestroyPetscObjects_();

    void SnapshotCurrentBface_();
    void RestoreCurrentBfaceFromSnapshot_();
    void BuildTrialBfaceFromUnknownE_();

    void UnpackVecToEhallField_(Vec X);
    void PackPredictedEhallToLocal_();

    void CopyEhallToE_();
    void ClearEdgeTriplet_(const IdTriplet &fid_triplet);
    void ClearFaceTriplet_(const IdTriplet &fid_triplet);

    void EvaluatePredictedEhallFromTrialB_();
    void WriteResidual_(Vec X, Vec F);

    // helpers
    void PackFaceInner_(int fid, std::vector<std::vector<double>> &buf);
    void RestoreFaceInner_(int fid, const std::vector<std::vector<double>> &buf);
    void AddFaceInnerFromRHS_(int fid_B, int fid_RHS, double factor);

public:
    double MaxAbsTriplet_(const IdTriplet &fid_triplet);
};

#endif