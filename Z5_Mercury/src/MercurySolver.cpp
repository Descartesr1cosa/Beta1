// Core
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

// Z4_Mercury
#include "MercurySolver.h"

MercurySolver::MercurySolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo,
                             Param *par
#ifdef HALL_IMPLICIT
                             ,
                             TOPO::TopologyEquiv *topo_equiv,
                             HALO_OWNER::EdgeOwnerSyncPattern *edge_owner_pat
#endif
                             )
    : grd_(grd),
      topo_(topo),
      fld_(fld),
      halo_(halo),
      par_(par)
#ifdef HALL_IMPLICIT
      ,
      topo_equiv_(topo_equiv),
      edge_owner_pat_(edge_owner_pat)
#endif
{
    // ---- Cache field ids ----
    fid_.Init(fld_);

    // ---- Build IO Module ----
    constexpr int NRES = 13; // 只统计H Na 守恒变量和感应磁场
    io_.Setup(par_, grd_, fld_, NRES);

    {
        std::vector<std::string> bin_name = {"U_H", "U_Na", "B_xi", "B_eta", "B_zeta"};
        io_.ClearRestartFields();
        io_.SetRestartFields(bin_name);

        io_.SetTecplotMode(IOModule::TecplotMode::CellAsNode);
        std::vector<std::string> tec_block_name = {}; // 全部物理块输出
        io_.SetTecplotBlock(tec_block_name);

        std::vector<std::string> plt_name = {"PV_H", "PV_Na", "B_cell", "Na", "U_plus"};
        io_.SetTecplotFields(plt_name);

        std::string fld_name = "PV_H";
        std::vector<std::string> var_name = {"u_H", "v_H", "w_H", "p_H", "T_H"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "PV_Na";
        var_name = {"u_Na", "v_Na", "w_Na", "p_Na", "T_Na"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "B_cell";
        var_name = {"Bx", "By", "Bz"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "Na";
        var_name = {"n_Na"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "U_plus";
        var_name = {"Up_u", "Up_v", "Up_w"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);
    }

    // ---- Calc Constants ----
    calc_physical_constant(par_);

    // ---- Boundary ----
    {
        // 0) 需要添加边界的物理场
        std::vector<std::string> bnd_fields = {
            "U_H",
            "U_Na",
            "B_xi",
            "B_eta",
            "B_zeta",
            "Badd_xi",
            "Badd_eta",
            "Badd_zeta",
            "B_cell",
            "J_xi",
            "J_eta",
            "J_zeta",
            "E_xi",
            "E_eta",
            "E_zeta",
            "Eface_xi",
            "Eface_eta",
            "Eface_zeta",
            "Ehall_xi",
            "Ehall_eta",
            "Ehall_zeta",
            "Bind_cell"};

        // 1) 初始化 Mercury Boundary
        mercury_bound_.Setup(grd_, fld_, topo_, halo_, par_, bnd_fields);
    }

    // ---- Initialization ----
    if (par_->GetBoo("continue_calc"))
    {
        io_.ReadRestartBinFile();
        io_.ReadRunDataFile();
    }
    initial_.Initialization(fld_, fid_);

    // ---- components ----
    control_.Setup(par_);

    run_data_ = &io_.Run();
    runtime_data_ = &io_.Runtime();

    auto count_global_cells = [&]() -> int64_t
    {
        int64_t local = 0;
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &U = fld_->field(fid_.fid_U_H, ib); // 任意 cell-centered 场都行
            if (!U.is_allocated())
                continue;
            const Int3 lo = U.inner_lo();
            const Int3 hi = U.inner_hi();
            local += int64_t(hi.i - lo.i) * int64_t(hi.j - lo.j) * int64_t(hi.k - lo.k);
        }
        double local_d = (double)local;
        double global = 0.0;
        PARALLEL::mpi_sum(&local_d, &global, 1);
        return (int64_t)std::llround(global);
    };

    runtime_data_->Begin(*run_data_, par_, count_global_cells());

#ifdef HALL_IMPLICIT
    if (!topo_equiv_ || !edge_owner_pat_)
        throw std::runtime_error("MercurySolver: hall implicit topology/pattern is null.");

    hall_implicit_.Setup(grd_, topo_, fld_, halo_, par_, &mercury_bound_,
                         fid_, *topo_equiv_, *edge_owner_pat_);

    ImplicitHallSolver::Callbacks cb;
    cb.sync_Bface = [this]()
    {
        mercury_bound_.Sync("Bface");
    };
    cb.sync_Eedge = [this]()
    {
        mercury_bound_.Sync("Eedge");
    };
    cb.calc_PV = [this]()
    {
        calc_PV();
    };
    cb.calc_Uplus = [this]()
    {
        calc_Uplus();
    };
    cb.build_Ehall_from_current_B = [this]()
    {
        Calc_J_Edge();
        calc_Jcell();
        AddHallEdgeEMF_();
    };

    hall_implicit_.SetCallbacks(cb);
    hall_implicit_.SetTheta(1.0); // midpoint
    hall_implicit_.InitializePetsc();
#endif
}
