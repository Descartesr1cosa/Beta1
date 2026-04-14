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

        std::vector<std::string> plt_name = {"PV_H", "PV_Na", "B_cell", "Na", "J_cell"};
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

        fld_name = "J_cell";
        var_name = {"Jx", "Jy", "Jz"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        // fld_name = "U_plus";
        // var_name = {"Up_u", "Up_v", "Up_w"};
        // io_.SetTecplotFieldComponentNames(fld_name, var_name);
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
            "Bind_cell",
            "J_cell",
            "dJ_xi",
            "dJ_eta",
            "dJ_zeta",
            "dE_xi",
            "dE_eta",
            "dE_zeta",
            "dB_xi",
            "dB_eta",
            "dB_zeta",
            "dJ_cell",
            "dEpre_xi",
            "dEpre_eta",
            "dEpre_zeta"};

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
                         fid_, *topo_equiv_, *edge_owner_pat_, &hall_face_scratch_);

    ImplicitHallSolver::Callbacks cb;
    cb.sync_Bface = [this]()
    {
        mercury_bound_.Sync("Bface");
    };
    cb.sync_Ehalledge = [this]()
    {
        mercury_bound_.Sync("Ehall");
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
        calc_Bcell();
        Calc_J_Edge();
        calc_Jcell();
        AddHallEdgeEMF_();
    };
    cb.calc_Bcell_from_current_Bface = [this]()
    {
        calc_Bcell();
    };
    cb.FillFrozenBflatFromCurrentBcell_ = [this]()
    { FillFrozenBflatFromCurrentBcell_(); };
    cb.FillFrozenAlphaFlatCell_ = [this]()
    { FillFrozenAlphaFlatCell_(); };
    cb.sync_dEedge = [this]()
    {
        mercury_bound_.Sync("dE");
    };

    cb.sync_dBface = [this]()
    {
        mercury_bound_.Sync("dB");
    };

    cb.sync_dJedge = [this]()
    {
        mercury_bound_.Sync("dJ");
    };
    cb.sync_dJcell = [this]()
    {
        mercury_bound_.Sync("dJcell");
    };
    cb.sync_dEface = [this]()
    {
        mercury_bound_.Sync("Eface");
    };

    hall_implicit_.SetCallbacks(cb);
    hall_implicit_.SetTheta(1.0); // BE //midpoint
    hall_implicit_.InitializePetsc();

    SetupHallFaceScratch_(); // setup temp block data for Rusanov Scheme
#endif
}

#ifdef HALL_IMPLICIT
void MercurySolver::SetupHallFaceScratch_()
{
    const int nb = fld_->num_blocks();
    hall_face_scratch_.clear();
    hall_face_scratch_.resize(nb);

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Bc = fld_->field(fid_.fid_Bcell, ib);
        if (!Bc.is_allocated())
            continue;

        Int3 clo = Bc.get_lo();
        Int3 chi = Bc.get_hi();

        const int ghost = -clo.i; // 前提：三个方向一致，且 clo = (-g,-g,-g)
        if (clo.i != -ghost || clo.j != -ghost || clo.k != -ghost)
        {
            throw std::runtime_error(
                "SetupHallFaceScratch_: Bcell lo is not compatible with Field_Array ghost indexing.");
        }

        const int dim1 = chi.i - clo.i;
        const int dim2 = chi.j - clo.j;
        const int dim3 = chi.k - clo.k;

        auto &buf = hall_face_scratch_[ib];
        buf.clo = clo;
        buf.chi = chi;
        buf.ni = dim1;
        buf.nj = dim2;
        buf.nk = dim3;

        buf.Ehc.SetSize(dim1, dim2, dim3, ghost, 3);
        buf.beta.SetSize(dim1, dim2, dim3, ghost);

        buf.Bflat.SetSize(dim1, dim2, dim3, ghost, 3);
        buf.alpha_flat.SetSize(dim1, dim2, dim3, ghost);
        buf.dEhc.SetSize(dim1, dim2, dim3, ghost, 3);
        // buf.beta_flat.SetSize(dim1, dim2, dim3, ghost);
        buf.dJcell_w.SetSize(dim1, dim2, dim3, ghost, 36);
        buf.dBcell_w.SetSize(dim1, dim2, dim3, ghost, 18);

        auto setup_like = [](auto &buf, auto &F, int ncomp)
        {
            if (!F.is_allocated())
                return;

            const Int3 lo = F.get_lo();
            const Int3 hi = F.get_hi();

            const int ghost = -lo.i;
            const int dim1 = hi.i - lo.i;
            const int dim2 = hi.j - lo.j;
            const int dim3 = hi.k - lo.k;

            buf.SetSize(dim1, dim2, dim3, ghost, ncomp);
        };
        // face projectors: 6 comps
        setup_like(buf.P_xi, fld_->field(fid_.fid_Eface.xi, ib), 6);
        setup_like(buf.P_eta, fld_->field(fid_.fid_Eface.eta, ib), 6);
        setup_like(buf.P_zeta, fld_->field(fid_.fid_Eface.zeta, ib), 6);

        // edge dr: 3 comps
        setup_like(buf.dr_xi, fld_->field(fid_.fid_Ehall.xi, ib), 3);
        setup_like(buf.dr_eta, fld_->field(fid_.fid_Ehall.eta, ib), 3);
        setup_like(buf.dr_zeta, fld_->field(fid_.fid_Ehall.zeta, ib), 3);
    }

    {
        constexpr double eps = 1e-30;

        // Bcell 比 Jcell 更保守一些：保护前 3 层（layer=0,1,2）
        const int protect_layers = 2;

        auto det3 = [](double a, double b, double c,
                       double d, double e, double f) -> double
        {
            return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
        };

        auto avg4 = [](double a, double b, double c, double d) -> double
        {
            return 0.25 * (a + b + c + d);
        };

        auto avg8 = [](double a0, double a1, double a2, double a3,
                       double a4, double a5, double a6, double a7) -> double
        {
            return 0.125 * (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7);
        };

        auto get_comp = [](const Int3 &a, int ax) -> int
        {
            return (ax == 0 ? a.i : (ax == 1 ? a.j : a.k));
        };

        auto set_comp = [](Int3 &a, int ax, int v)
        {
            if (ax == 0)
                a.i = v;
            else if (ax == 1)
                a.j = v;
            else
                a.k = v;
        };

        // ------------------------------------------------------------
        // 对 cell-centered Bcell:
        // layer=0 就是紧贴 Pole 的 first off-axis shell
        // ------------------------------------------------------------
        auto pole_layer_of_cell =
            [&](int ib, int i, int j, int k, int &layer_out, int &ax_out, int &dir_out) -> bool
        {
            bool found = false;
            int best_layer = 1e9;
            int best_ax = -1;
            int best_dir = 0;

            for (const auto &p : topo_->physical_patches)
            {
                if (p.this_block != ib)
                    continue;
                if (p.bc_name != "Pole")
                    continue;

                const int ax = std::abs(p.direction) - 1;
                const int dir = p.direction;

                const Int3 nlo = p.this_box_node.lo;
                const Int3 nhi = p.this_box_node.hi;

                bool inside_tangent = true;
                for (int d = 0; d < 3; ++d)
                {
                    if (d == ax)
                        continue;

                    const int coord = (d == 0 ? i : (d == 1 ? j : k));
                    const int clo = get_comp(nlo, d);
                    const int chi = get_comp(nhi, d) - 1; // cell-centered half-open
                    if (!(clo <= coord && coord < chi))
                    {
                        inside_tangent = false;
                        break;
                    }
                }
                if (!inside_tangent)
                    continue;

                const int coord_ax = (ax == 0 ? i : (ax == 1 ? j : k));

                int layer = -1;
                if (dir < 0)
                {
                    const int plane0 = get_comp(nlo, ax);
                    layer = coord_ax - plane0;
                }
                else
                {
                    // plus side: first inside cell layer is node_hi - 2
                    const int plane0 = get_comp(nhi, ax) - 2;
                    layer = plane0 - coord_ax;
                }

                if (layer < 0)
                    continue;

                if (layer < best_layer)
                {
                    best_layer = layer;
                    best_ax = ax;
                    best_dir = dir;
                    found = true;
                }
            }

            if (found)
            {
                layer_out = best_layer;
                ax_out = best_ax;
                dir_out = best_dir;
            }
            return found;
        };

        auto copy_W = [&](auto &W, int i_dst, int j_dst, int k_dst,
                          int i_src, int j_src, int k_src)
        {
            for (int m = 0; m < 18; ++m)
                W(i_dst, j_dst, k_dst, m) = W(i_src, j_src, k_src, m);
        };

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Bc = fld_->field(fid_.fid_Bcell, ib);
            if (!Bc.is_allocated())
                continue;

            auto &buf = hall_face_scratch_[ib];
            auto &W = buf.dBcell_w;

            auto &JDxi = fld_->field(fid_.fid_metric.xi, ib);
            auto &JDet = fld_->field(fid_.fid_metric.eta, ib);
            auto &JDze = fld_->field(fid_.fid_metric.zeta, ib);

            auto &x = grd_->grids(ib).x;
            auto &y = grd_->grids(ib).y;
            auto &z = grd_->grids(ib).z;

            Int3 lo = Bc.inner_lo();
            Int3 hi = Bc.inner_hi();

            // ============================================================
            // Pass 1:
            //   regular build with face screening
            // ============================================================
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // ----------------------------------------------------
                        // 1) cell center
                        // ----------------------------------------------------
                        const double xc = avg8(
                            x(i, j, k), x(i + 1, j, k), x(i, j + 1, k), x(i + 1, j + 1, k),
                            x(i, j, k + 1), x(i + 1, j, k + 1), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));

                        const double yc = avg8(
                            y(i, j, k), y(i + 1, j, k), y(i, j + 1, k), y(i + 1, j + 1, k),
                            y(i, j, k + 1), y(i + 1, j, k + 1), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));

                        const double zc = avg8(
                            z(i, j, k), z(i + 1, j, k), z(i, j + 1, k), z(i + 1, j + 1, k),
                            z(i, j, k + 1), z(i + 1, j, k + 1), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

                        // ----------------------------------------------------
                        // 2) 6 face centers
                        // ----------------------------------------------------
                        double fx[6], fy[6], fz[6];

                        fx[0] = avg4(x(i, j, k), x(i, j + 1, k), x(i, j, k + 1), x(i, j + 1, k + 1));
                        fy[0] = avg4(y(i, j, k), y(i, j + 1, k), y(i, j, k + 1), y(i, j + 1, k + 1));
                        fz[0] = avg4(z(i, j, k), z(i, j + 1, k), z(i, j, k + 1), z(i, j + 1, k + 1));

                        fx[1] = avg4(x(i + 1, j, k), x(i + 1, j + 1, k), x(i + 1, j, k + 1), x(i + 1, j + 1, k + 1));
                        fy[1] = avg4(y(i + 1, j, k), y(i + 1, j + 1, k), y(i + 1, j, k + 1), y(i + 1, j + 1, k + 1));
                        fz[1] = avg4(z(i + 1, j, k), z(i + 1, j + 1, k), z(i + 1, j, k + 1), z(i + 1, j + 1, k + 1));

                        fx[2] = avg4(x(i, j, k), x(i + 1, j, k), x(i, j, k + 1), x(i + 1, j, k + 1));
                        fy[2] = avg4(y(i, j, k), y(i + 1, j, k), y(i, j, k + 1), y(i + 1, j, k + 1));
                        fz[2] = avg4(z(i, j, k), z(i + 1, j, k), z(i, j, k + 1), z(i + 1, j, k + 1));

                        fx[3] = avg4(x(i, j + 1, k), x(i + 1, j + 1, k), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));
                        fy[3] = avg4(y(i, j + 1, k), y(i + 1, j + 1, k), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));
                        fz[3] = avg4(z(i, j + 1, k), z(i + 1, j + 1, k), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

                        fx[4] = avg4(x(i, j, k), x(i + 1, j, k), x(i, j + 1, k), x(i + 1, j + 1, k));
                        fy[4] = avg4(y(i, j, k), y(i + 1, j, k), y(i, j + 1, k), y(i + 1, j + 1, k));
                        fz[4] = avg4(z(i, j, k), z(i + 1, j, k), z(i, j + 1, k), z(i + 1, j + 1, k));

                        fx[5] = avg4(x(i, j, k + 1), x(i + 1, j, k + 1), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));
                        fy[5] = avg4(y(i, j, k + 1), y(i + 1, j, k + 1), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));
                        fz[5] = avg4(z(i, j, k + 1), z(i + 1, j, k + 1), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

                        // ----------------------------------------------------
                        // 3) face normals / areas / weights
                        // ----------------------------------------------------
                        double nx[6], ny[6], nz[6], Smag[6], invS[6], ww[6];

                        auto fill_face = [&](int idx,
                                             double Sx, double Sy, double Sz,
                                             double fcx, double fcy, double fcz)
                        {
                            const double S = std::sqrt(Sx * Sx + Sy * Sy + Sz * Sz) + eps;
                            Smag[idx] = S;
                            nx[idx] = Sx / S;
                            ny[idx] = Sy / S;
                            nz[idx] = Sz / S;
                            invS[idx] = 1.0 / S;

                            const double dx = fcx - xc;
                            const double dy = fcy - yc;
                            const double dz = fcz - zc;
                            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                            ww[idx] = 1.0 / std::max(dist, eps);
                        };

                        fill_face(0, -JDxi(i, j, k, 0), -JDxi(i, j, k, 1), -JDxi(i, j, k, 2), fx[0], fy[0], fz[0]);
                        fill_face(1, JDxi(i + 1, j, k, 0), JDxi(i + 1, j, k, 1), JDxi(i + 1, j, k, 2), fx[1], fy[1], fz[1]);

                        fill_face(2, -JDet(i, j, k, 0), -JDet(i, j, k, 1), -JDet(i, j, k, 2), fx[2], fy[2], fz[2]);
                        fill_face(3, JDet(i, j + 1, k, 0), JDet(i, j + 1, k, 1), JDet(i, j + 1, k, 2), fx[3], fy[3], fz[3]);

                        fill_face(4, -JDze(i, j, k, 0), -JDze(i, j, k, 1), -JDze(i, j, k, 2), fx[4], fy[4], fz[4]);
                        fill_face(5, JDze(i, j, k + 1, 0), JDze(i, j, k + 1, 1), JDze(i, j, k + 1, 2), fx[5], fy[5], fz[5]);

                        // local area reference
                        double Sref = 0.0;
                        for (int n = 0; n < 6; ++n)
                            Sref = std::max(Sref, Smag[n]);
                        Sref = std::max(Sref, 1.0e-20);

                        const double Smin = 1.0e-3 * Sref;

                        int nvalid = 0;
                        for (int n = 0; n < 6; ++n)
                        {
                            if (Smag[n] <= Smin)
                                ww[n] = 0.0;
                            if (ww[n] > 0.0)
                                ++nvalid;
                        }

                        // 先清零
                        for (int m = 0; m < 18; ++m)
                            W(i, j, k, m) = 0.0;

                        // 有效 face 太少，留给 Pass 2 fallback
                        if (nvalid < 4)
                            continue;

                        // ----------------------------------------------------
                        // 4) N = A^T W A
                        // ----------------------------------------------------
                        double N00 = 0.0, N01 = 0.0, N02 = 0.0;
                        double N11 = 0.0, N12 = 0.0, N22 = 0.0;

                        for (int n = 0; n < 6; ++n)
                        {
                            const double w = ww[n];
                            if (w == 0.0)
                                continue;

                            N00 += w * nx[n] * nx[n];
                            N01 += w * nx[n] * ny[n];
                            N02 += w * nx[n] * nz[n];
                            N11 += w * ny[n] * ny[n];
                            N12 += w * ny[n] * nz[n];
                            N22 += w * nz[n] * nz[n];
                        }

                        double det = det3(N00, N01, N02, N11, N12, N22);

                        // B 这里正则也加大一点，别太弱
                        const double reg = 1e-10 * (N00 + N11 + N22 + 1.0);

                        if (std::abs(det) < reg)
                        {
                            N00 += reg;
                            N11 += reg;
                            N22 += reg;
                            det = det3(N00, N01, N02, N11, N12, N22);
                        }

                        if (std::abs(det) < 1e-24)
                            continue;

                        const double C00 = (N11 * N22 - N12 * N12);
                        const double C01 = (N02 * N12 - N01 * N22);
                        const double C02 = (N01 * N12 - N02 * N11);
                        const double C11 = (N00 * N22 - N02 * N02);
                        const double C12 = (N01 * N02 - N00 * N12);
                        const double C22 = (N00 * N11 - N01 * N01);

                        const double invdet = 1.0 / det;

                        const double inv00 = C00 * invdet;
                        const double inv01 = C01 * invdet;
                        const double inv02 = C02 * invdet;
                        const double inv11 = C11 * invdet;
                        const double inv12 = C12 * invdet;
                        const double inv22 = C22 * invdet;

                        // ----------------------------------------------------
                        // 5) 3x6 权重
                        // ----------------------------------------------------
                        for (int n = 0; n < 6; ++n)
                        {
                            if (ww[n] == 0.0)
                                continue;

                            const double c0 = ww[n] * nx[n] * invS[n];
                            const double c1 = ww[n] * ny[n] * invS[n];
                            const double c2 = ww[n] * nz[n] * invS[n];

                            W(i, j, k, n) = inv00 * c0 + inv01 * c1 + inv02 * c2;      // Bx
                            W(i, j, k, 6 + n) = inv01 * c0 + inv11 * c1 + inv12 * c2;  // By
                            W(i, j, k, 12 + n) = inv02 * c0 + inv12 * c1 + inv22 * c2; // Bz
                        }
                    }

            // ============================================================
            // Pass 2:
            //   Pole 邻近前 protect_layers+1 层直接从更内层拷贝 W
            // ============================================================
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        int layer, ax, dir;
                        if (!pole_layer_of_cell(ib, i, j, k, layer, ax, dir))
                            continue;

                        if (layer > protect_layers)
                            continue;

                        const int safe_layer = protect_layers + 1;

                        int i_src = i, j_src = j, k_src = k;
                        if (dir < 0)
                        {
                            if (ax == 0)
                                i_src = i + (safe_layer - layer);
                            if (ax == 1)
                                j_src = j + (safe_layer - layer);
                            if (ax == 2)
                                k_src = k + (safe_layer - layer);
                        }
                        else
                        {
                            if (ax == 0)
                                i_src = i - (safe_layer - layer);
                            if (ax == 1)
                                j_src = j - (safe_layer - layer);
                            if (ax == 2)
                                k_src = k - (safe_layer - layer);
                        }

                        i_src = std::max(lo.i, std::min(i_src, hi.i - 1));
                        j_src = std::max(lo.j, std::min(j_src, hi.j - 1));
                        k_src = std::max(lo.k, std::min(k_src, hi.k - 1));

                        copy_W(W, i, j, k, i_src, j_src, k_src);
                    }
        }
    }

    // {
    //     constexpr double eps = 1e-30;

    //     auto det3 = [](double a, double b, double c,
    //                    double d, double e, double f) -> double
    //     {
    //         return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    //     };

    //     auto avg4 = [](double a, double b, double c, double d) -> double
    //     {
    //         return 0.25 * (a + b + c + d);
    //     };

    //     auto avg8 = [](double a0, double a1, double a2, double a3,
    //                    double a4, double a5, double a6, double a7) -> double
    //     {
    //         return 0.125 * (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7);
    //     };

    //     for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    //     {
    //         auto &Bc = fld_->field(fid_.fid_Bcell, ib);
    //         if (!Bc.is_allocated())
    //             continue;

    //         auto &buf = hall_face_scratch_[ib];
    //         auto &W = buf.dBcell_w;

    //         auto &JDxi = fld_->field(fid_.fid_metric.xi, ib);
    //         auto &JDet = fld_->field(fid_.fid_metric.eta, ib);
    //         auto &JDze = fld_->field(fid_.fid_metric.zeta, ib);

    //         auto &x = grd_->grids(ib).x;
    //         auto &y = grd_->grids(ib).y;
    //         auto &z = grd_->grids(ib).z;

    //         Int3 lo = Bc.inner_lo();
    //         Int3 hi = Bc.inner_hi();

    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     // ----------------------------------------------------
    //                     // 1) cell center (8 nodes average)
    //                     // ----------------------------------------------------
    //                     const double xc = avg8(
    //                         x(i, j, k), x(i + 1, j, k), x(i, j + 1, k), x(i + 1, j + 1, k),
    //                         x(i, j, k + 1), x(i + 1, j, k + 1), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));

    //                     const double yc = avg8(
    //                         y(i, j, k), y(i + 1, j, k), y(i, j + 1, k), y(i + 1, j + 1, k),
    //                         y(i, j, k + 1), y(i + 1, j, k + 1), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));

    //                     const double zc = avg8(
    //                         z(i, j, k), z(i + 1, j, k), z(i, j + 1, k), z(i + 1, j + 1, k),
    //                         z(i, j, k + 1), z(i + 1, j, k + 1), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

    //                     // ----------------------------------------------------
    //                     // 2) 6 face centers
    //                     //    顺序: xi-, xi+, eta-, eta+, zeta-, zeta+
    //                     // ----------------------------------------------------
    //                     double fx[6], fy[6], fz[6];

    //                     // xi-
    //                     fx[0] = avg4(x(i, j, k), x(i, j + 1, k), x(i, j, k + 1), x(i, j + 1, k + 1));
    //                     fy[0] = avg4(y(i, j, k), y(i, j + 1, k), y(i, j, k + 1), y(i, j + 1, k + 1));
    //                     fz[0] = avg4(z(i, j, k), z(i, j + 1, k), z(i, j, k + 1), z(i, j + 1, k + 1));

    //                     // xi+
    //                     fx[1] = avg4(x(i + 1, j, k), x(i + 1, j + 1, k), x(i + 1, j, k + 1), x(i + 1, j + 1, k + 1));
    //                     fy[1] = avg4(y(i + 1, j, k), y(i + 1, j + 1, k), y(i + 1, j, k + 1), y(i + 1, j + 1, k + 1));
    //                     fz[1] = avg4(z(i + 1, j, k), z(i + 1, j + 1, k), z(i + 1, j, k + 1), z(i + 1, j + 1, k + 1));

    //                     // eta-
    //                     fx[2] = avg4(x(i, j, k), x(i + 1, j, k), x(i, j, k + 1), x(i + 1, j, k + 1));
    //                     fy[2] = avg4(y(i, j, k), y(i + 1, j, k), y(i, j, k + 1), y(i + 1, j, k + 1));
    //                     fz[2] = avg4(z(i, j, k), z(i + 1, j, k), z(i, j, k + 1), z(i + 1, j, k + 1));

    //                     // eta+
    //                     fx[3] = avg4(x(i, j + 1, k), x(i + 1, j + 1, k), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));
    //                     fy[3] = avg4(y(i, j + 1, k), y(i + 1, j + 1, k), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));
    //                     fz[3] = avg4(z(i, j + 1, k), z(i + 1, j + 1, k), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

    //                     // zeta-
    //                     fx[4] = avg4(x(i, j, k), x(i + 1, j, k), x(i, j + 1, k), x(i + 1, j + 1, k));
    //                     fy[4] = avg4(y(i, j, k), y(i + 1, j, k), y(i, j + 1, k), y(i + 1, j + 1, k));
    //                     fz[4] = avg4(z(i, j, k), z(i + 1, j, k), z(i, j + 1, k), z(i + 1, j + 1, k));

    //                     // zeta+
    //                     fx[5] = avg4(x(i, j, k + 1), x(i + 1, j, k + 1), x(i, j + 1, k + 1), x(i + 1, j + 1, k + 1));
    //                     fy[5] = avg4(y(i, j, k + 1), y(i + 1, j, k + 1), y(i, j + 1, k + 1), y(i + 1, j + 1, k + 1));
    //                     fz[5] = avg4(z(i, j, k + 1), z(i + 1, j, k + 1), z(i, j + 1, k + 1), z(i + 1, j + 1, k + 1));

    //                     // ----------------------------------------------------
    //                     // 3) 6 face normals / areas / weights
    //                     // ----------------------------------------------------
    //                     double nx[6], ny[6], nz[6], invS[6], ww[6];

    //                     auto fill_face = [&](int idx,
    //                                          double Sx, double Sy, double Sz,
    //                                          double fcx, double fcy, double fcz)
    //                     {
    //                         const double Smag = std::sqrt(Sx * Sx + Sy * Sy + Sz * Sz) + eps;
    //                         nx[idx] = Sx / Smag;
    //                         ny[idx] = Sy / Smag;
    //                         nz[idx] = Sz / Smag;
    //                         invS[idx] = 1.0 / Smag;

    //                         const double dx = fcx - xc;
    //                         const double dy = fcy - yc;
    //                         const double dz = fcz - zc;
    //                         const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    //                         // 这里就沿用你 calc_Bcell 里的距离权定义
    //                         ww[idx] = 1.0 / std::max(dist, eps);
    //                     };

    //                     // xi- / xi+
    //                     fill_face(0,
    //                               -JDxi(i, j, k, 0), -JDxi(i, j, k, 1), -JDxi(i, j, k, 2),
    //                               fx[0], fy[0], fz[0]);

    //                     fill_face(1,
    //                               JDxi(i + 1, j, k, 0), JDxi(i + 1, j, k, 1), JDxi(i + 1, j, k, 2),
    //                               fx[1], fy[1], fz[1]);

    //                     // eta- / eta+
    //                     fill_face(2,
    //                               -JDet(i, j, k, 0), -JDet(i, j, k, 1), -JDet(i, j, k, 2),
    //                               fx[2], fy[2], fz[2]);

    //                     fill_face(3,
    //                               JDet(i, j + 1, k, 0), JDet(i, j + 1, k, 1), JDet(i, j + 1, k, 2),
    //                               fx[3], fy[3], fz[3]);

    //                     // zeta- / zeta+
    //                     fill_face(4,
    //                               -JDze(i, j, k, 0), -JDze(i, j, k, 1), -JDze(i, j, k, 2),
    //                               fx[4], fy[4], fz[4]);

    //                     fill_face(5,
    //                               JDze(i, j, k + 1, 0), JDze(i, j, k + 1, 1), JDze(i, j, k + 1, 2),
    //                               fx[5], fy[5], fz[5]);

    //                     // ----------------------------------------------------
    //                     // 4) N = A^T W A
    //                     // ----------------------------------------------------
    //                     double N00 = 0.0, N01 = 0.0, N02 = 0.0;
    //                     double N11 = 0.0, N12 = 0.0, N22 = 0.0;

    //                     for (int n = 0; n < 6; ++n)
    //                     {
    //                         const double w = ww[n];
    //                         N00 += w * nx[n] * nx[n];
    //                         N01 += w * nx[n] * ny[n];
    //                         N02 += w * nx[n] * nz[n];
    //                         N11 += w * ny[n] * ny[n];
    //                         N12 += w * ny[n] * nz[n];
    //                         N22 += w * nz[n] * nz[n];
    //                     }

    //                     double det = det3(N00, N01, N02, N11, N12, N22);
    //                     const double reg = 1e-14 * (N00 + N11 + N22 + 1.0);

    //                     if (std::abs(det) < reg)
    //                     {
    //                         N00 += reg;
    //                         N11 += reg;
    //                         N22 += reg;
    //                         det = det3(N00, N01, N02, N11, N12, N22);
    //                     }

    //                     const double C00 = (N11 * N22 - N12 * N12);
    //                     const double C01 = (N02 * N12 - N01 * N22);
    //                     const double C02 = (N01 * N12 - N02 * N11);
    //                     const double C11 = (N00 * N22 - N02 * N02);
    //                     const double C12 = (N01 * N02 - N00 * N12);
    //                     const double C22 = (N00 * N11 - N01 * N01);

    //                     const double invdet = 1.0 / det;

    //                     const double inv00 = C00 * invdet;
    //                     const double inv01 = C01 * invdet;
    //                     const double inv02 = C02 * invdet;
    //                     const double inv11 = C11 * invdet;
    //                     const double inv12 = C12 * invdet;
    //                     const double inv22 = C22 * invdet;

    //                     // ----------------------------------------------------
    //                     // 5) 预计算 3x6 权重，直接吸收 invS
    //                     //    B = Wpre * Phi
    //                     // ----------------------------------------------------
    //                     for (int n = 0; n < 6; ++n)
    //                     {
    //                         const double c0 = ww[n] * nx[n] * invS[n];
    //                         const double c1 = ww[n] * ny[n] * invS[n];
    //                         const double c2 = ww[n] * nz[n] * invS[n];

    //                         W(i, j, k, n) = inv00 * c0 + inv01 * c1 + inv02 * c2;      // Bx
    //                         W(i, j, k, 6 + n) = inv01 * c0 + inv11 * c1 + inv12 * c2;  // By
    //                         W(i, j, k, 12 + n) = inv02 * c0 + inv12 * c1 + inv22 * c2; // Bz
    //                     }
    //                 }
    //     }
    // }

    {
        constexpr double eps = 1e-25;

        // near-axis 保护层数：
        // protect_layers = 1  => layer 0,1 都不用本地病态W，改为从更内层拷贝
        const int protect_layers = 1;

        auto det3 = [](double a, double b, double c,
                       double d, double e, double f) -> double
        {
            return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
        };

        auto get_comp = [](const Int3 &a, int ax) -> int
        {
            return (ax == 0 ? a.i : (ax == 1 ? a.j : a.k));
        };

        auto set_comp = [](Int3 &a, int ax, int v)
        {
            if (ax == 0)
                a.i = v;
            else if (ax == 1)
                a.j = v;
            else
                a.k = v;
        };

        // ------------------------------------------------------------
        // 判断某个 cell 是否在 Pole 邻近层
        // 返回:
        //   true  -> 找到了 Pole patch，layer/ax/dir_sign 有效
        //   false -> 普通区域
        //
        // 注意:
        // 对 cell-centered 量，
        // layer=0 就是紧贴 Pole 边界的第一层 cell，也就是 first off-axis shell
        // ------------------------------------------------------------
        auto pole_layer_of_cell =
            [&](int ib, int i, int j, int k, int &layer_out, int &ax_out, int &dir_out) -> bool
        {
            bool found = false;
            int best_layer = 1e9;
            int best_ax = -1;
            int best_dir = 0;

            for (const auto &p : topo_->physical_patches)
            {
                if (p.this_block != ib)
                    continue;
                if (p.bc_name != "Pole")
                    continue;

                const int ax = std::abs(p.direction) - 1;
                const int dir = p.direction;

                const Int3 nlo = p.this_box_node.lo;
                const Int3 nhi = p.this_box_node.hi;

                // tangential cell range: [node_lo, node_hi-1)
                bool inside_tangent = true;
                for (int d = 0; d < 3; ++d)
                {
                    if (d == ax)
                        continue;

                    const int coord = (d == 0 ? i : (d == 1 ? j : k));
                    const int clo = get_comp(nlo, d);
                    const int chi = get_comp(nhi, d) - 1; // half-open for cell-centered
                    if (!(clo <= coord && coord < chi))
                    {
                        inside_tangent = false;
                        break;
                    }
                }
                if (!inside_tangent)
                    continue;

                const int coord_ax = (ax == 0 ? i : (ax == 1 ? j : k));

                int layer = -1;
                if (dir < 0)
                {
                    // minus side: first cell layer is node_lo + 0
                    const int plane0 = get_comp(nlo, ax);
                    layer = coord_ax - plane0;
                }
                else
                {
                    // plus side: first cell layer is node_hi - 2
                    const int plane0 = get_comp(nhi, ax) - 2;
                    layer = plane0 - coord_ax;
                }

                if (layer < 0)
                    continue;

                if (layer < best_layer)
                {
                    best_layer = layer;
                    best_ax = ax;
                    best_dir = dir;
                    found = true;
                }
            }

            if (found)
            {
                layer_out = best_layer;
                ax_out = best_ax;
                dir_out = best_dir;
            }
            return found;
        };

        // ------------------------------------------------------------
        // 从更内层 cell 复制一整套 3x12 权重
        // ------------------------------------------------------------
        auto copy_W = [&](auto &W, int i_dst, int j_dst, int k_dst,
                          int i_src, int j_src, int k_src)
        {
            for (int m = 0; m < 36; ++m)
                W(i_dst, j_dst, k_dst, m) = W(i_src, j_src, k_src, m);
        };

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Bc = fld_->field(fid_.fid_Bcell, ib);
            if (!Bc.is_allocated())
                continue;

            auto &buf = hall_face_scratch_[ib];
            auto &W = buf.dJcell_w;

            auto &dl_xi = fld_->field("dl_xi", ib);
            auto &dl_eta = fld_->field("dl_eta", ib);
            auto &dl_zeta = fld_->field("dl_zeta", ib);

            auto &alpha_xi = fld_->field(fid_.Edge_alpha.xi, ib);
            auto &alpha_eta = fld_->field(fid_.Edge_alpha.eta, ib);
            auto &alpha_zeta = fld_->field(fid_.Edge_alpha.zeta, ib);

            auto &x = grd_->grids(ib).x;
            auto &y = grd_->grids(ib).y;
            auto &z = grd_->grids(ib).z;

            auto edge_xi = [&](int i, int j, int k,
                               double &tx, double &ty, double &tz,
                               double &Lraw, double &alpha)
            {
                Lraw = dl_xi(i, j, k, 0);
                alpha = alpha_xi(i, j, k, 0);

                const double Linv = 1.0 / std::max(Lraw, eps);
                tx = (x(i + 1, j, k) - x(i, j, k)) * Linv;
                ty = (y(i + 1, j, k) - y(i, j, k)) * Linv;
                tz = (z(i + 1, j, k) - z(i, j, k)) * Linv;
            };

            auto edge_eta = [&](int i, int j, int k,
                                double &tx, double &ty, double &tz,
                                double &Lraw, double &alpha)
            {
                Lraw = dl_eta(i, j, k, 0);
                alpha = alpha_eta(i, j, k, 0);

                const double Linv = 1.0 / std::max(Lraw, eps);
                tx = (x(i, j + 1, k) - x(i, j, k)) * Linv;
                ty = (y(i, j + 1, k) - y(i, j, k)) * Linv;
                tz = (z(i, j + 1, k) - z(i, j, k)) * Linv;
            };

            auto edge_zeta = [&](int i, int j, int k,
                                 double &tx, double &ty, double &tz,
                                 double &Lraw, double &alpha)
            {
                Lraw = dl_zeta(i, j, k, 0);
                alpha = alpha_zeta(i, j, k, 0);

                const double Linv = 1.0 / std::max(Lraw, eps);
                tx = (x(i, j, k + 1) - x(i, j, k)) * Linv;
                ty = (y(i, j, k + 1) - y(i, j, k)) * Linv;
                tz = (z(i, j, k + 1) - z(i, j, k)) * Linv;
            };

            Int3 lo = Bc.inner_lo();
            Int3 hi = Bc.inner_hi();

            // ============================================================
            // Pass 1:
            //   先对所有 cell 构造 regular W
            //   但屏蔽退化 edge；如果矩阵太差，则先置零
            // ============================================================
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double tx[12], ty[12], tz[12];
                        double Lraw[12], alpha_e[12];
                        double ww[12];

                        // 0..3 xi
                        edge_xi(i, j, k, tx[0], ty[0], tz[0], Lraw[0], alpha_e[0]);
                        edge_xi(i, j + 1, k, tx[1], ty[1], tz[1], Lraw[1], alpha_e[1]);
                        edge_xi(i, j, k + 1, tx[2], ty[2], tz[2], Lraw[2], alpha_e[2]);
                        edge_xi(i, j + 1, k + 1, tx[3], ty[3], tz[3], Lraw[3], alpha_e[3]);

                        // 4..7 eta
                        edge_eta(i, j, k, tx[4], ty[4], tz[4], Lraw[4], alpha_e[4]);
                        edge_eta(i + 1, j, k, tx[5], ty[5], tz[5], Lraw[5], alpha_e[5]);
                        edge_eta(i, j, k + 1, tx[6], ty[6], tz[6], Lraw[6], alpha_e[6]);
                        edge_eta(i + 1, j, k + 1, tx[7], ty[7], tz[7], Lraw[7], alpha_e[7]);

                        // 8..11 zeta
                        edge_zeta(i, j, k, tx[8], ty[8], tz[8], Lraw[8], alpha_e[8]);
                        edge_zeta(i + 1, j, k, tx[9], ty[9], tz[9], Lraw[9], alpha_e[9]);
                        edge_zeta(i, j + 1, k, tx[10], ty[10], tz[10], Lraw[10], alpha_e[10]);
                        edge_zeta(i + 1, j + 1, k, tx[11], ty[11], tz[11], Lraw[11], alpha_e[11]);

                        // local reference length
                        double href = 0.0;
                        for (int n = 0; n < 12; ++n)
                            href = std::max(href, Lraw[n]);
                        href = std::max(href, 1.0e-12);

                        const double Lmin = 1.0e-3 * href;

                        int nvalid = 0;
                        for (int n = 0; n < 12; ++n)
                        {
                            const bool ok = (alpha_e[n] > 0.0) && (Lraw[n] > Lmin);
                            ww[n] = ok ? 1.0 : 0.0;
                            if (ok)
                                ++nvalid;
                        }

                        // 默认先清零
                        for (int m = 0; m < 36; ++m)
                            W(i, j, k, m) = 0.0;

                        // 有效样本太少，直接留待 Pass 2 fallback
                        if (nvalid < 6)
                            continue;

                        double N00 = 0.0, N01 = 0.0, N02 = 0.0;
                        double N11 = 0.0, N12 = 0.0, N22 = 0.0;

                        for (int n = 0; n < 12; ++n)
                        {
                            const double w = ww[n];
                            if (w == 0.0)
                                continue;

                            N00 += w * tx[n] * tx[n];
                            N01 += w * tx[n] * ty[n];
                            N02 += w * tx[n] * tz[n];
                            N11 += w * ty[n] * ty[n];
                            N12 += w * ty[n] * tz[n];
                            N22 += w * tz[n] * tz[n];
                        }

                        double det = det3(N00, N01, N02, N11, N12, N22);

                        // 正则稍微放大一些，别再用过于弱的 1e-14
                        const double reg = 1e-10 * (N00 + N11 + N22 + 1.0);

                        if (std::abs(det) < reg)
                        {
                            N00 += reg;
                            N11 += reg;
                            N22 += reg;
                            det = det3(N00, N01, N02, N11, N12, N22);
                        }

                        // 还是太差，留给 Pass 2 fallback
                        if (std::abs(det) < 1e-20)
                            continue;

                        const double C00 = (N11 * N22 - N12 * N12);
                        const double C01 = (N02 * N12 - N01 * N22);
                        const double C02 = (N01 * N12 - N02 * N11);
                        const double C11 = (N00 * N22 - N02 * N02);
                        const double C12 = (N01 * N02 - N00 * N12);
                        const double C22 = (N00 * N11 - N01 * N01);

                        const double invdet = 1.0 / det;

                        const double inv00 = C00 * invdet;
                        const double inv01 = C01 * invdet;
                        const double inv02 = C02 * invdet;
                        const double inv11 = C11 * invdet;
                        const double inv12 = C12 * invdet;
                        const double inv22 = C22 * invdet;

                        for (int n = 0; n < 12; ++n)
                        {
                            if (ww[n] == 0.0)
                                continue;

                            const double invL = 1.0 / std::max(Lraw[n], eps);

                            const double c0 = ww[n] * tx[n] * invL;
                            const double c1 = ww[n] * ty[n] * invL;
                            const double c2 = ww[n] * tz[n] * invL;

                            W(i, j, k, n) = inv00 * c0 + inv01 * c1 + inv02 * c2;      // Jx
                            W(i, j, k, 12 + n) = inv01 * c0 + inv11 * c1 + inv12 * c2; // Jy
                            W(i, j, k, 24 + n) = inv02 * c0 + inv12 * c1 + inv22 * c2; // Jz
                        }
                    }

            // ============================================================
            // Pass 2:
            //   Pole 邻近层(layer <= protect_layers) 直接从更内层 safe layer 拷贝 W
            //   这是最保守、最容易先跑稳的做法
            // ============================================================
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        int layer, ax, dir;
                        if (!pole_layer_of_cell(ib, i, j, k, layer, ax, dir))
                            continue;

                        if (layer > protect_layers)
                            continue;

                        // safe source plane: first unprotected layer
                        const int safe_layer = protect_layers + 1;

                        int i_src = i, j_src = j, k_src = k;
                        if (dir < 0)
                        {
                            if (ax == 0)
                                i_src = i + (safe_layer - layer);
                            else if (ax == 1)
                                j_src = j + (safe_layer - layer);
                            else
                                k_src = k + (safe_layer - layer);
                        }
                        else
                        {
                            if (ax == 0)
                                i_src = i - (safe_layer - layer);
                            else if (ax == 1)
                                j_src = j - (safe_layer - layer);
                            else
                                k_src = k - (safe_layer - layer);
                        }

                        // clamp
                        i_src = std::max(lo.i, std::min(i_src, hi.i - 1));
                        j_src = std::max(lo.j, std::min(j_src, hi.j - 1));
                        k_src = std::max(lo.k, std::min(k_src, hi.k - 1));

                        copy_W(W, i, j, k, i_src, j_src, k_src);
                    }
        }
    }

    // {
    //     constexpr double eps = 1e-25;

    //     auto det3 = [](double a, double b, double c,
    //                    double d, double e, double f) -> double
    //     {
    //         return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    //     };

    //     for (int ib = 0; ib < nb; ++ib)
    //     {
    //         auto &Bc = fld_->field(fid_.fid_Bcell, ib);
    //         if (!Bc.is_allocated())
    //             continue;

    //         auto &buf = hall_face_scratch_[ib];
    //         auto &W = buf.dJcell_w;

    //         auto &dl_xi = fld_->field("dl_xi", ib);
    //         auto &dl_eta = fld_->field("dl_eta", ib);
    //         auto &dl_zeta = fld_->field("dl_zeta", ib);

    //         auto &x = grd_->grids(ib).x;
    //         auto &y = grd_->grids(ib).y;
    //         auto &z = grd_->grids(ib).z;

    //         auto unit_t_xi = [&](int i, int j, int k,
    //                              double &tx, double &ty, double &tz, double &invL)
    //         {
    //             const double L = std::max(dl_xi(i, j, k, 0), eps);
    //             invL = 1.0 / L;
    //             tx = (x(i + 1, j, k) - x(i, j, k)) * invL;
    //             ty = (y(i + 1, j, k) - y(i, j, k)) * invL;
    //             tz = (z(i + 1, j, k) - z(i, j, k)) * invL;
    //         };

    //         auto unit_t_eta = [&](int i, int j, int k,
    //                               double &tx, double &ty, double &tz, double &invL)
    //         {
    //             const double L = std::max(dl_eta(i, j, k, 0), eps);
    //             invL = 1.0 / L;
    //             tx = (x(i, j + 1, k) - x(i, j, k)) * invL;
    //             ty = (y(i, j + 1, k) - y(i, j, k)) * invL;
    //             tz = (z(i, j + 1, k) - z(i, j, k)) * invL;
    //         };

    //         auto unit_t_zeta = [&](int i, int j, int k,
    //                                double &tx, double &ty, double &tz, double &invL)
    //         {
    //             const double L = std::max(dl_zeta(i, j, k, 0), eps);
    //             invL = 1.0 / L;
    //             tx = (x(i, j, k + 1) - x(i, j, k)) * invL;
    //             ty = (y(i, j, k + 1) - y(i, j, k)) * invL;
    //             tz = (z(i, j, k + 1) - z(i, j, k)) * invL;
    //         };

    //         Int3 lo = Bc.inner_lo();
    //         Int3 hi = Bc.inner_hi();

    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     double tx[12], ty[12], tz[12], invL[12];
    //                     double ww[12];
    //                     for (int n = 0; n < 12; ++n)
    //                         ww[n] = 1.0;

    //                     // 0..3 xi
    //                     unit_t_xi(i, j, k, tx[0], ty[0], tz[0], invL[0]);
    //                     unit_t_xi(i, j + 1, k, tx[1], ty[1], tz[1], invL[1]);
    //                     unit_t_xi(i, j, k + 1, tx[2], ty[2], tz[2], invL[2]);
    //                     unit_t_xi(i, j + 1, k + 1, tx[3], ty[3], tz[3], invL[3]);

    //                     // 4..7 eta
    //                     unit_t_eta(i, j, k, tx[4], ty[4], tz[4], invL[4]);
    //                     unit_t_eta(i + 1, j, k, tx[5], ty[5], tz[5], invL[5]);
    //                     unit_t_eta(i, j, k + 1, tx[6], ty[6], tz[6], invL[6]);
    //                     unit_t_eta(i + 1, j, k + 1, tx[7], ty[7], tz[7], invL[7]);

    //                     // 8..11 zeta
    //                     unit_t_zeta(i, j, k, tx[8], ty[8], tz[8], invL[8]);
    //                     unit_t_zeta(i + 1, j, k, tx[9], ty[9], tz[9], invL[9]);
    //                     unit_t_zeta(i, j + 1, k, tx[10], ty[10], tz[10], invL[10]);
    //                     unit_t_zeta(i + 1, j + 1, k, tx[11], ty[11], tz[11], invL[11]);

    //                     // N = A^T W A
    //                     double N00 = 0.0, N01 = 0.0, N02 = 0.0;
    //                     double N11 = 0.0, N12 = 0.0, N22 = 0.0;

    //                     for (int n = 0; n < 12; ++n)
    //                     {
    //                         const double w = ww[n];
    //                         N00 += w * tx[n] * tx[n];
    //                         N01 += w * tx[n] * ty[n];
    //                         N02 += w * tx[n] * tz[n];
    //                         N11 += w * ty[n] * ty[n];
    //                         N12 += w * ty[n] * tz[n];
    //                         N22 += w * tz[n] * tz[n];
    //                     }

    //                     double det = det3(N00, N01, N02, N11, N12, N22);
    //                     const double reg = 1e-14 * (N00 + N11 + N22 + 1.0);

    //                     if (std::abs(det) < reg)
    //                     {
    //                         N00 += reg;
    //                         N11 += reg;
    //                         N22 += reg;
    //                         det = det3(N00, N01, N02, N11, N12, N22);
    //                     }

    //                     const double C00 = (N11 * N22 - N12 * N12);
    //                     const double C01 = (N02 * N12 - N01 * N22);
    //                     const double C02 = (N01 * N12 - N02 * N11);
    //                     const double C11 = (N00 * N22 - N02 * N02);
    //                     const double C12 = (N01 * N02 - N00 * N12);
    //                     const double C22 = (N00 * N11 - N01 * N01);

    //                     const double invdet = 1.0 / det;

    //                     const double inv00 = C00 * invdet;
    //                     const double inv01 = C01 * invdet;
    //                     const double inv02 = C02 * invdet;
    //                     const double inv11 = C11 * invdet;
    //                     const double inv12 = C12 * invdet;
    //                     const double inv22 = C22 * invdet;

    //                     // 预计算 3x12 权重，直接吸收 invL
    //                     for (int n = 0; n < 12; ++n)
    //                     {
    //                         const double c0 = ww[n] * tx[n] * invL[n];
    //                         const double c1 = ww[n] * ty[n] * invL[n];
    //                         const double c2 = ww[n] * tz[n] * invL[n];

    //                         W(i, j, k, n) = inv00 * c0 + inv01 * c1 + inv02 * c2;      // Jx
    //                         W(i, j, k, 12 + n) = inv01 * c0 + inv11 * c1 + inv12 * c2; // Jy
    //                         W(i, j, k, 24 + n) = inv02 * c0 + inv12 * c1 + inv22 * c2; // Jz
    //                     }
    //                 }
    //     }
    // }

    {
        constexpr double eps = 1e-30;

        auto fill_projector = [&](auto &Pbuf, auto &Sfield)
        {
            if (!Sfield.is_allocated())
                return;

            Int3 lo = Sfield.get_lo();
            Int3 hi = Sfield.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double Sx = Sfield(i, j, k, 0);
                        const double Sy = Sfield(i, j, k, 1);
                        const double Sz = Sfield(i, j, k, 2);

                        const double Smag = std::sqrt(Sx * Sx + Sy * Sy + Sz * Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        const double Pxx = 1.0 - nx * nx;
                        const double Pxy = -nx * ny;
                        const double Pxz = -nx * nz;
                        const double Pyy = 1.0 - ny * ny;
                        const double Pyz = -ny * nz;
                        const double Pzz = 1.0 - nz * nz;

                        Pbuf(i, j, k, 0) = Pxx;
                        Pbuf(i, j, k, 1) = Pxy;
                        Pbuf(i, j, k, 2) = Pxz;
                        Pbuf(i, j, k, 3) = Pyy;
                        Pbuf(i, j, k, 4) = Pyz;
                        Pbuf(i, j, k, 5) = Pzz;
                    }
        };

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &buf = hall_face_scratch_[ib];
            fill_projector(buf.P_xi, fld_->field(fid_.fid_metric.xi, ib));
            fill_projector(buf.P_eta, fld_->field(fid_.fid_metric.eta, ib));
            fill_projector(buf.P_zeta, fld_->field(fid_.fid_metric.zeta, ib));
        }
    }
}
#endif
