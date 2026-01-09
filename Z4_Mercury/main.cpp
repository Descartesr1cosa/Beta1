//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

#include "MercurySolver.h"
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
    // 守恒变量、独立变量
    fld->register_field({"U_H", StaggerLocation::Cell, 5, ngg, "Fluid"});  // H+
    fld->register_field({"U_Na", StaggerLocation::Cell, 5, ngg, "Fluid"}); // Na+
    fld->register_field({"U_b", StaggerLocation::Cell, 3, ngg});           // induced magnetic fields

    // 辅助物理场
    fld->register_field({"Na", StaggerLocation::Cell, 1, ngg, "Fluid"});         // Na neutral atom
    fld->register_field({"Photo_rate", StaggerLocation::Cell, 1, ngg, "Fluid"}); // Photoionization rate
    fld->register_field({"Badd", StaggerLocation::Cell, 3, ngg});                // initial applied magnetic fields
    fld->register_field({"B_cell", StaggerLocation::Cell, 3, ngg});              // 总磁场
    fld->register_field({"U_plus", StaggerLocation::Cell, 3, ngg, "Fluid"});     // 按电荷密度加权的平均速度

    fld->register_field(FieldDescriptor{"PV_H", StaggerLocation::Cell, 5, ngg, "Fluid"});  // H+的原始变量 u v w p T
    fld->register_field(FieldDescriptor{"PV_Na", StaggerLocation::Cell, 5, ngg, "Fluid"}); // Na+的原始变量 u v w p T

    // 辅助通量场
    fld->register_field({"F_xi", StaggerLocation::FaceXi, 5, 0, "Fluid"});   // 仅临存流体方程通量，无需虚网格
    fld->register_field({"F_eta", StaggerLocation::FaceEt, 5, 0, "Fluid"});  // 仅临存流体方程通量，无需虚网格
    fld->register_field({"F_zeta", StaggerLocation::FaceZe, 5, 0, "Fluid"}); // 仅临存流体方程通量，无需虚网格

    // 计算辅助场
    // fld->register_field(FieldDescriptor{"old_U_", StaggerLocation::Cell, 5, 0});
    // fld->register_field(FieldDescriptor{"divB", StaggerLocation::Cell, 1, 1}); // 输出存在插值到node的需要，增加一层虚网格
    fld->register_field(FieldDescriptor{"RHS_H", StaggerLocation::Cell, 5, 0, "Fluid"});
    fld->register_field(FieldDescriptor{"RHS_Na", StaggerLocation::Cell, 5, 0, "Fluid"});
    fld->register_field(FieldDescriptor{"RHS_B", StaggerLocation::Cell, 3, 0});
    //--------------------------------------------------------------------------
    // 注册耦合定义（CouplingPairDesc）
    fld->register_coupling_channel("Solid", "Fluid", "U_b", StaggerLocation::Cell, 3, ngg); // Solid -> Fluid
    fld->register_coupling_channel("Fluid", "Solid", "U_b", StaggerLocation::Cell, 3, ngg); // Fluid -> Solid
    // 构建 coupling buffers（一次）
    fld->build_coupling_buffers(topology, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    // 建立Halo通信
    Halo *hal = new Halo(fld, &topology);

    // 注册同物理场之间的halo通信
    std::string fieldname;
    fieldname = "U_H";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "U_Na";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "U_b";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    // 建立同物理场以及多物理场耦合通信的特征
    hal->build_registered_patterns();
    //=============================================================================================

    //=============================================================================================
    MercurySolver solver(grd, &topology, fld, hal, par);
    solver.Advance();
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