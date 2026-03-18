#pragma once

#include "5_io/IOModule.h"

#include "1_Boundary.h"
#include "0_SolverFields.h"
#include "2_Initial.h"
#include "3_Control.h"

// ---- forward declarations (avoid heavy includes in header) ----
class Grid;
namespace TOPO
{
    class Topology;
}
class Field;
class Halo;
class Param;
class FieldBlock;

class MercurySolver
{
public:
    MercurySolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par);

    void Advance();

private:
    // --- core pointers ---
    Grid *grd_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Field *fld_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};

    // --- components ---
    Control control_;
    MercuryBoundary mercury_bound_;
    IOModule io_;
    Mercury_Initial initial_;
    RunData *run_data_;
    RuntimeMonitor *runtime_data_;

    // --- cached field ids  ---
    SolverFields fid_;

    // --- constants ---
    double gamma_{0.0};
    double NA{0.0};
    double R_uni{0.0};
    double k_Boltz{0.0};
    double q_e{0.0};
    double mu0{0.0};
    double dt{0.0};
    double dt_hall{0.0};
    double dt_sub{0.0};
    double ne_hall_floor{0.0};
    double ne_hall_floor_dimensional{0.0};

    double U_ref{0.0};
    double L_ref{0.0};
    double B_ref{0.0};
    double T_ref{0.0};
    double n_ref{0.0};
    double rho_ref{0.0};
    double M_ref{0.0};
    double M_H{0.0};
    double M_Na{0.0};
    double m_H{0.0};
    double m_Na{0.0};
    double state_coeff_H{0.0};
    double state_coeff_Na{0.0};
    double CFL{0.0};
    double hall_coef{0.0};
    double ambi_coef{0.0};

    double momentum_induce_coeff{0.0};
    double momentum_hall_coeff{0.0};

    double inver_MA2{0.0};
    double inver_Rem{0.0};

private:
    //=========================================================================
    // TOOLS
    void calc_Bcell();
    void calc_Jcell();
    void calc_divB();
    void calc_PV();
    void calc_Uplus();
    void calc_physical_constant(Param *par);
    void PrintMinMaxDiagnostics_();
    double HallAlpha_Coeffient(double ne_true, double r);
    void Hall_Num_Limiter(double rhoH, double rhoNa, double *num);
    //=========================================================================

    //=========================================================================
    bool StepOnce();
    //---------------------------------------------------------------
    void Compute_Timestep();
    bool UpdateControlAndOutput();
    //=========================================================================

    //=========================================================================
    void Time_Advance();
    //---------------------------------------------------------------
    void ZeroRHS_();
    void AssembleRHS_Fluid_();
    void AssembleRHS_Induction_CT_();
    void ApplyUpdate_Euler_();
    //---------------------------------------------------------------
    // For Fluid
    void Scheme_U_();
    void AddSourceToRHS_Fluid();
    //---------------------------------------------------------------
    // For Magnetic
    void Build_E_explicit_edge_();
    void AddResistiveEdgeEMF_();
    void AddIdealEdgeEMF_();
    void AddHallEdgeEMF_();
    void AddAmbipolarEdgeEMF_();
    void Calc_J_Edge();

    // 只组装 Hall 的 RHS_b（不动 U 的 RHS）
    void AssembleRHS_Induction_CT_HallOnly_();
    // 只更新 Bface: Bface += dt_sub * RHS_b
    void ApplyUpdate_Euler_BfaceOnly_(double dt_sub);

    void BuildHallFaceEMF_Rusanov_();
    void AssembleEdgeEMF_FromFaceE_Hall_();
    //--------------------------------
    //  For Ideal
    void AssembleOneDirectionEMF_(int iblk, int dir, FieldBlock &E_face, FieldBlock &B_face, FieldBlock &B_face_add, FieldBlock &Bcell, FieldBlock &metricField, FieldBlock &Uplus, FieldBlock &UH, FieldBlock &UN);
    void AssembleEdgeEMF_FromFaceE_Ideal_();
    void ReconstructionEMF_(double *metric, int32_t direction,
                            FieldBlock &Uplus, FieldBlock &UH, FieldBlock &UN, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k,
                            double *out_flux);
    //---------------------------------------------------------------

    //=========================================================================

private:
    // //=========================================================================
    // void BC_UH_Farfield_Na(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    // // void BC_UH_Farfield_b(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh) {};
    // void BC_UH_Farfield_H(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    // void BC_Solid_Surface(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);

    // // void Time_Advance();
    // //---------------------------------------------------------------
    // void ZeroRHS_();
    // void Scheme_U_();
    // void Scheme_B_();
    // void AddSourceToRHS_B();
    // void AddSourceToRHS_Fluid();
    // void ApplyUpdate_Euler_();
    // void calc_Uplus();
    // //=========================================================================

private:
    //     // =============================== Step driver ============================
    //     bool StepOnce();
    //     //-------------------------------------------------------------------------
    //     void Compute_Timestep();
    //     void Calc_Residual();
    //     bool UpdateControlAndOutput();
    //     //=========================================================================

    //     // ======================== Prepare / Sync pipeline =======================
    //     void PrepareStep();
    //     void PrepareSubstep_NoSnapshot();
    //     //-------------------------------------------------------------------------
    //     void SyncPrimaryFaceB();       // B_face_*（ghost）
    //     void ComputeBcellInner();      // 在 inner 域从 B_face_* 重建 cell 磁场。
    //     void SyncDerivedBcell();       // 更新 B_cell 的 ghost
    //     void SyncPrimaryCellU();       // 守恒变量U添加边界条件
    //     void UpdateDerivedPVandDivB(); // 计算原始变量等被动量
    //     void SnapshotOldFields();      // 拷贝保存当前场，用于残差计算
    //     //=========================================================================

    //     // ============================= RHS assembly =============================
    //     void Time_Advance();
    //     //-------------------------------------------------------------------------
    //     void ZeroRHS();
    //     void AssembleRHS_Fluid();     // inv_fluid()
    //     void SyncElectricFace();      // E_face_* BC+halo
    //     void AssembleRHS_Induction(); // inv_induce()
    //     void ApplyTimeUpdate_Euler(); // += dt*RHS
    //     void Update_Physic_Time();    // record physical time
    //     //=========================================================================

    //     // ========================== helper for Fluid ============================
    //     void AssembleOneDirectionFluxAndEMF_(int iblk,
    //                                          int dir,                 // 0 xi, 1 eta, 2 zeta
    //                                          FieldBlock &flux,        // F_xi / F_eta / F_zeta (ncomp=5)
    //                                          FieldBlock &E_face,      // E_face_xi/eta/zeta   (ncomp=3)
    //                                          FieldBlock &B_face,      // B_xi/eta/zeta        (ncomp=1)
    //                                          FieldBlock &B_face_add,  // B_xi/eta/zeta add        (ncomp=1)
    //                                          FieldBlock &metricField, // Xi_/Eta_/Zeta_       (ncomp=3)
    //                                          FieldBlock &PV,
    //                                          FieldBlock &U,
    //                                          FieldBlock &Bcell);
    //     void AssembleCellRHSFromFlux_();
    //     // Reconstruction / flux
    //     void Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k, double *out_flux);
    //     //=========================================================================

    //     // ========================= helper for Induction =========================
    //     void AssembleEdgeEMF_FromFaceE_Ideal_();
    //     void AddExplicitHallToEdgeEMF_(); // 只在 hall_explicit.cpp 实现由宏控制, 对于Ideal Implicit均为空
    //     void ApplyBC_EdgeEMF_();
    //     void AssembleFaceRHS_FromEdgeEMF_Curl_();

    //     // ================================== TOOLS ==============================

    //     void calc_Bcell();
    //     void calc_divB();
    //     void copy_field();
    //     void PrintMinMaxDiagnostics_();
    //     void add_Emag_to_Etotal();
    //     double ComputeLocalMaxRadius_();
    //     //=========================================================================
};