//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"

// #include "MercurySolver.h"
// #include "4_solver/ImplicitHall_Solver.h"
//==============================================================================

//==============================================================================

//==============================================================================

int main(int arg, char **argv)
{
    //=============================================================================================
    // MPI initialization
    int myid;
    PARALLEL::mpi_initial(arg, argv);
    PARALLEL::mpi_rank(&myid);
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // 读入控制参数
    Param *par = new Param;
    par->ReadParam(myid);
    //--------------------------------------------------------------------------
    // 读入网格并作预处理
    Grid *grd = new Grid;
    grd->Grid_Preprocess(par);
    //--------------------------------------------------------------------------
    // 建立topology
    TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    // 建立Field
    Field *fld = new Field(grd, par);
    //-------------------------------------
    // 加入求解物理场
    int ngg = par->GetInt("ngg");
    // 守恒变量、独立变量，用于构建CT方法
    fld->register_field({"U_H", StaggerLocation::Cell, 5, ngg});  // H+
    fld->register_field({"U_Na", StaggerLocation::Cell, 5, ngg}); // Na+
    fld->register_field({"U_b", StaggerLocation::Cell, 3, ngg});  // induced magnetic fields

    // 辅助物理场
    fld->register_field({"Badd", StaggerLocation::Cell, 3, ngg});   // initial applied magnetic fields
    fld->register_field({"B_cell", StaggerLocation::Cell, 3, ngg}); // 总磁场

    fld->register_field(FieldDescriptor{"PV_H", StaggerLocation::Cell, 5, ngg});  // H+的原始变量 u v w p T
    fld->register_field(FieldDescriptor{"PV_Na", StaggerLocation::Cell, 5, ngg}); // Na+的原始变量 u v w p T

    // 辅助通量场
    fld->register_field({"F_xi", StaggerLocation::FaceXi, 5, 0});   // 仅临存流体方程通量，无需虚网格
    fld->register_field({"F_eta", StaggerLocation::FaceEt, 5, 0});  // 仅临存流体方程通量，无需虚网格
    fld->register_field({"F_zeta", StaggerLocation::FaceZe, 5, 0}); // 仅临存流体方程通量，无需虚网格

    // 辅助电场变量
    // fld->register_field({"E_xi", StaggerLocation::EdgeXi, 1, 1});   // 用于CT方法，只需要1层虚网格
    // fld->register_field({"E_eta", StaggerLocation::EdgeEt, 1, 1});  // 用于CT方法，只需要1层虚网格
    // fld->register_field({"E_zeta", StaggerLocation::EdgeZe, 1, 1}); // 用于CT方法，只需要1层虚网格

    // 计算辅助场
    // fld->register_field(FieldDescriptor{"old_U_", StaggerLocation::Cell, 5, 0});
    // fld->register_field(FieldDescriptor{"divB", StaggerLocation::Cell, 1, 1}); // 输出存在插值到node的需要，增加一层虚网格
    // fld->register_field(FieldDescriptor{"old_B_xi", StaggerLocation::FaceXi, 1, 0});
    // fld->register_field(FieldDescriptor{"old_B_eta", StaggerLocation::FaceEt, 1, 0});
    // fld->register_field(FieldDescriptor{"old_B_zeta", StaggerLocation::FaceZe, 1, 0});
    // fld->register_field(FieldDescriptor{"RHS", StaggerLocation::Cell, 5, 0});
    // fld->register_field(FieldDescriptor{"RHS_xi", StaggerLocation::FaceXi, 1, 0});
    // fld->register_field(FieldDescriptor{"RHS_eta", StaggerLocation::FaceEt, 1, 0});
    // fld->register_field(FieldDescriptor{"RHS_zeta", StaggerLocation::FaceZe, 1, 0});
    //--------------------------------------------------------------------------
    // 建立Halo通信
    Halo *hal = new Halo(fld, &topology);
    //=============================================================================================

    //=============================================================================================
    // HallMHDSolver solver(grd, &topology, fld, hal, par, hall_imp); //, {"U_", "B_xi", "B_eta", "B_zeta"});
    // solver.Advance();
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI终止
    PARALLEL::mpi_finalize();
    //--------------------------------------------------------------------------
    // 释放所分配的空间，建议按照创建顺序逆序释放
    delete hal;
    delete fld;
    delete par;
    delete grd;
    //=============================================================================================
    return 0;
}