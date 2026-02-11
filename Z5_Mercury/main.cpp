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
    // Read control parameters
    Param *par = new Param;
    par->ReadParam(myid);
    //--------------------------------------------------------------------------
    // Read Grid and Preprocess the Grid related info
    Grid *grd = new Grid;
    grd->Grid_Preprocess(par);
    //--------------------------------------------------------------------------
    // Build topology
    TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    int ngg = par->GetInt("ngg");
    // Build Field
    Field *fld = new Field(grd, par, ngg);
    //-------------------------------------
    // Add physical fields for solver
    // Conservative variables, independent
    fld->register_field({"U_H", StaggerLocation::Cell, 5, ngg, "Fluid"});  // H+  (rho rho.u rho.v rho.w rho.e)_H+
    fld->register_field({"U_Na", StaggerLocation::Cell, 5, ngg, "Fluid"}); // Na+ (rho rho.u rho.v rho.w rho.e)_Na+
    fld->register_field({"B_xi", StaggerLocation::FaceXi, 1, ngg});        // induced magnetic flux for Face Xi
    fld->register_field({"B_eta", StaggerLocation::FaceEt, 1, ngg});       // induced magnetic flux for Face Eta
    fld->register_field({"B_zeta", StaggerLocation::FaceZe, 1, ngg});      // induced magnetic flux for Face Zeta

    // Auxiliary physical fields
    fld->register_field({"E_xi", StaggerLocation::EdgeXi, 1, ngg});       // Integration of electric field along Edge Xi
    fld->register_field({"E_eta", StaggerLocation::EdgeEt, 1, ngg});      // Integration of electric field along Edge Eta
    fld->register_field({"E_zeta", StaggerLocation::EdgeZe, 1, ngg});     // Integration of electric field along Edge Zeta
    fld->register_field({"Ehall_xi", StaggerLocation::EdgeXi, 1, ngg});   // Integration of electric field along Edge Xi
    fld->register_field({"Ehall_eta", StaggerLocation::EdgeEt, 1, ngg});  // Integration of electric field along Edge Eta
    fld->register_field({"Ehall_zeta", StaggerLocation::EdgeZe, 1, ngg}); // Integration of electric field along Edge Zeta
    fld->register_field({"Eface_xi", StaggerLocation::FaceXi, 3, ngg});   // Integration of electric field on Face Xi For CT
    fld->register_field({"Eface_eta", StaggerLocation::FaceEt, 3, ngg});  // Integration of electric field on Face Eta For CT
    fld->register_field({"Eface_zeta", StaggerLocation::FaceZe, 3, ngg}); // Integration of electric field on Face Zeta For CT

    fld->register_field({"J_xi", StaggerLocation::EdgeXi, 1, ngg});   // Integration of electric current along Edge Xi
    fld->register_field({"J_eta", StaggerLocation::EdgeEt, 1, ngg});  // Integration of electric current along Edge Eta
    fld->register_field({"J_zeta", StaggerLocation::EdgeZe, 1, ngg}); // Integration of electric current along Edge Zeta

    fld->register_field({"Badd_xi", StaggerLocation::FaceXi, 1, ngg});   // initial applied/added magnetic flux for Face Xi
    fld->register_field({"Badd_eta", StaggerLocation::FaceEt, 1, ngg});  // initial applied/added magnetic flux for Face Eta
    fld->register_field({"Badd_zeta", StaggerLocation::FaceZe, 1, ngg}); // initial applied/added magnetic flux for Face Zeta
    fld->register_field({"B_cell", StaggerLocation::Cell, 3, ngg});      // Total magnetic fields

    fld->register_field({"Na", StaggerLocation::Cell, 1, ngg, "Fluid"});         // Na neutral atom
    fld->register_field({"Photo_rate", StaggerLocation::Cell, 1, ngg, "Fluid"}); // Photoionization rate
    fld->register_field({"U_plus", StaggerLocation::Cell, 3, ngg, "Fluid"});     // Averaged Velocity (electric density weighted), used in induction Eqs

    fld->register_field(FieldDescriptor{"PV_H", StaggerLocation::Cell, 5, ngg, "Fluid"});  // H+  primitive variables: u v w p T
    fld->register_field(FieldDescriptor{"PV_Na", StaggerLocation::Cell, 5, ngg, "Fluid"}); // Na+ primitive variables: u v w p T

    // Auxiliary flux fields, for Solver
    fld->register_field({"F_xi", StaggerLocation::FaceXi, 5, 0, "Fluid"});   // Only stores flux of Fluid Equations, ghost grid is not required
    fld->register_field({"F_eta", StaggerLocation::FaceEt, 5, 0, "Fluid"});  // Only stores flux of Fluid Equations, ghost grid is not required
    fld->register_field({"F_zeta", StaggerLocation::FaceZe, 5, 0, "Fluid"}); // Only stores flux of Fluid Equations, ghost grid is not required

    // Auxiliary fields, for Solver
    fld->register_field(FieldDescriptor{"divB", StaggerLocation::Cell, 1, 1}); // 1 layer of ghost grid for interpolation of Cell to Node (output and visualization)
    fld->register_field(FieldDescriptor{"RHS_H", StaggerLocation::Cell, 5, 0, "Fluid"});
    fld->register_field(FieldDescriptor{"RHS_Na", StaggerLocation::Cell, 5, 0, "Fluid"});
    fld->register_field(FieldDescriptor{"RHS_B_xi", StaggerLocation::FaceXi, 1, 0});
    fld->register_field(FieldDescriptor{"RHS_B_eta", StaggerLocation::FaceEt, 1, 0});
    fld->register_field(FieldDescriptor{"RHS_B_zeta", StaggerLocation::FaceZe, 1, 0});
    //--------------------------------------------------------------------------
    // Register Coupling Pair Description（CouplingPairDesc）
    //   register_coupling_channel("A", "B", "A_field",**):
    //   Let A_field in Block A transfer to coresponding coupling buffer area of Block B
    fld->register_coupling_channel("Solid", "Fluid", "B_xi", StaggerLocation::FaceXi, 1, ngg);       // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "B_eta", StaggerLocation::FaceEt, 1, ngg);      // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "B_zeta", StaggerLocation::FaceZe, 1, ngg);     // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Eface_xi", StaggerLocation::FaceXi, 3, ngg);   // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Eface_eta", StaggerLocation::FaceEt, 3, ngg);  // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Eface_zeta", StaggerLocation::FaceZe, 3, ngg); // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "E_xi", StaggerLocation::EdgeXi, 1, ngg);       // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "E_eta", StaggerLocation::EdgeEt, 1, ngg);      // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "E_zeta", StaggerLocation::EdgeZe, 1, ngg);     // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Ehall_xi", StaggerLocation::EdgeXi, 1, ngg);   // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Ehall_eta", StaggerLocation::EdgeEt, 1, ngg);  // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "Ehall_zeta", StaggerLocation::EdgeZe, 1, ngg); // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "J_xi", StaggerLocation::EdgeXi, 1, ngg);       // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "J_eta", StaggerLocation::EdgeEt, 1, ngg);      // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "J_zeta", StaggerLocation::EdgeZe, 1, ngg);     // Solid -> Fluid
    fld->register_coupling_channel("Solid", "Fluid", "B_cell", StaggerLocation::Cell, 3, ngg);       // Solid -> Fluid

    fld->register_coupling_channel("Fluid", "Solid", "B_xi", StaggerLocation::FaceXi, 1, ngg);       // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "B_eta", StaggerLocation::FaceEt, 1, ngg);      // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "B_zeta", StaggerLocation::FaceZe, 1, ngg);     // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Eface_xi", StaggerLocation::FaceXi, 3, ngg);   // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Eface_eta", StaggerLocation::FaceEt, 3, ngg);  // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Eface_zeta", StaggerLocation::FaceZe, 3, ngg); // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "E_xi", StaggerLocation::EdgeXi, 1, ngg);       // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "E_eta", StaggerLocation::EdgeEt, 1, ngg);      // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "E_zeta", StaggerLocation::EdgeZe, 1, ngg);     // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Ehall_xi", StaggerLocation::EdgeXi, 1, ngg);   // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Ehall_eta", StaggerLocation::EdgeEt, 1, ngg);  // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "Ehall_zeta", StaggerLocation::EdgeZe, 1, ngg); // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "J_xi", StaggerLocation::EdgeXi, 1, ngg);       // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "J_eta", StaggerLocation::EdgeEt, 1, ngg);      // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "J_zeta", StaggerLocation::EdgeZe, 1, ngg);     // Fluid -> Solid
    fld->register_coupling_channel("Fluid", "Solid", "B_cell", StaggerLocation::Cell, 3, ngg);       // Fluid -> Solid
    // Build coupling buffers (YOU CAN ONLY USE it ONCE!)
    fld->build_coupling_buffers(topology, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    // Build Halo Communicator
    Halo *hal = new Halo(fld, &topology);

    // Register halo communicator between blocks with same fields
    std::string fieldname;
    fieldname = "U_H";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "U_Na";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "B_xi";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "B_eta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "B_zeta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "E_xi";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "E_eta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "E_zeta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "J_xi";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "J_eta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    fieldname = "J_zeta";
    hal->register_halo_field(fieldname, HaloLevel::Vertex);
    // Build halo communicator patterns between blocks with same fields and coupling fields
    hal->build_registered_patterns();
    //=============================================================================================

    //=============================================================================================
    MercurySolver solver(grd, &topology, fld, hal, par);
    solver.Advance();
    if (myid == 0)
        std::cout << "Program is finished normally ! !  ^_^\n"
                  << std::flush;
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI finalization
    PARALLEL::mpi_finalize();
    //--------------------------------------------------------------------------
    // Release the allocated memory, recommend to release with inverse order of building
    delete hal;
    delete fld;
    delete par;
    delete grd;
    //=============================================================================================
    return 0;
}