#include "1_Boundary.h"
#include "0_basic/Error.h"
#include "00_Lunar_Const.h"

void LunarBoundary::InstallHandlers()
{
    if (!par_)
        ERROR::Abort("InstallHandlers: call Setup first");

    bool is_Lunar_interior_resis = par_->GetBoo("is_Lunar_resistance");

    auto copy = [](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto nop = [](FieldBlock &, Field *, const BOUND::PhysicalRegion &, int) {};

    // -------------------------------------------------------------------------
    // 0. Default physical handlers
    //
    // Install default physical behavior, then override fields that need a
    // dedicated lunar surface rule.
    //
    // 后面只覆盖真正需要特殊处理的字段。
    // -------------------------------------------------------------------------
    for (const auto &fn : boundary_fields_)
    {
        RegisterPhysical_(fn, "Outflow", copy);
        RegisterPhysical_(fn, "Farfield", copy);
        RegisterPhysical_(fn, "Pole", copy);
        RegisterPhysical_(fn, "Solid_Surface", copy);
    }

    RegisterPhysical_("U_plus", "Outflow", nop);
    RegisterPhysical_("U_plus", "Farfield", nop);
    RegisterPhysical_("U_plus", "Pole", nop);
    RegisterPhysical_("U_plus", "Solid_Surface", nop);

    // ============================================================================================
    // Fluid BCs
    // -------------------------------------------------------------------------
    //  Fluid farfield
    // -------------------------------------------------------------------------
    RegisterPhysical_("U_H", "Farfield",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_UH_Farfield_H_(U, fld, r, ngh);
                      });

    // -------------------------------------------------------------------------
    //  Solid wall: fluid variables only
    // -------------------------------------------------------------------------
    RegisterPhysical_("U_H", "Solid_Surface",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Solid_Surface_(U, fld, r, ngh);
                      });

    // -------------------------------------------------------------------------
    //  Pole: fluid cells
    // -------------------------------------------------------------------------
    RegisterPhysical_("U_H", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Cell_(U, fld, r, ngh);
                      });

    // ============================================================================================

    // ============================================================================================
    // Electromagnetic BCs
    // -------------------------------------------------------------------------
    // Cell: B_cell J_cell
    // -------------------------------------------------------------------------
    RegisterPhysical_("B_cell", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh);
                      });

    RegisterPhysical_("Bind_cell", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh);
                      });
    RegisterPhysical_("B_cell", "Farfield", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Farfield_Bface(U, fld, r, ngh); });
    RegisterPhysical_("Bind_cell", "Farfield", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Farfield_Bface(U, fld, r, ngh); });
    // J_cell Pole无需特别处理，在计算的时候就已经处理过了
    // J_cell 壁面层 cell 复制域内相邻一层，虚网格再从该壁面层拷贝
    if (!is_Lunar_interior_resis)
    {
        RegisterPhysical_("J_cell", "Solid_Surface",
                          [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                          {
                              this->BC_Solid_Surface_Jcell(U, fld, r, ngh);
                          });
    }
    // -------------------------------------------------------------------------
    // Face: B_xi/eta/zeta Eface_xi/eta/zeta
    // -------------------------------------------------------------------------
    // B_xi/eta/zeta应该严格按照CT求解的，因此绝大部分不需要施加边界条件
    // Pole的正则性要求，需要对Pole上的进行处理
    // ------------------------------------------
    // Pole: face magnetic flux
    // B_xi / B_eta 使用 collapse。
    // B_zeta 目前显式保持 copy。不要在 Bcell handler 里顺手改 B_zeta。
    // ------------------------------------------
    auto Bface_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Bface_Collapse_(U, fld, r, ngh);
    };
    RegisterPhysical_("B_xi", "Pole", Bface_pole);
    RegisterPhysical_("B_eta", "Pole", Bface_pole);
    // 当前策略：B_zeta copy。
    RegisterPhysical_("B_zeta", "Pole", copy);

    RegisterPhysical_("B_xi", "Farfield", nop);
    RegisterPhysical_("B_eta", "Farfield", nop);
    RegisterPhysical_("B_zeta", "Farfield", nop);
    // ------------------------------------------
    // E_face
    // ------------------------------------------
    if (!is_Lunar_interior_resis)
    {
        auto Eface_Wall = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
        {
            this->BC_Solid_Surface_Eface_(U, fld, r, ngh);
        };
        RegisterPhysical_("Eface_xi", "Solid_Surface", Eface_Wall);
        RegisterPhysical_("Eface_eta", "Solid_Surface", Eface_Wall);
        RegisterPhysical_("Eface_zeta", "Solid_Surface", Eface_Wall);
    }
    else
    {
        auto Eface_Wall = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
        {
            this->BC_Solid_Surface_Eface_ghots_zero(U, fld, r, ngh);
        };
        RegisterPhysical_("Eface_xi", "Solid_Surface", Eface_Wall);
        RegisterPhysical_("Eface_eta", "Solid_Surface", Eface_Wall);
        RegisterPhysical_("Eface_zeta", "Solid_Surface", Eface_Wall);
    }

    RegisterPhysical_("Eface_xi", "Farfield", nop);
    RegisterPhysical_("Eface_eta", "Farfield", nop);
    RegisterPhysical_("Eface_zeta", "Farfield", nop);
    // -------------------------------------------------------------------------
    // Edge: E_xi/eta/zeta
    // -------------------------------------------------------------------------
    auto Eedge_1stlayer_to_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        // no operators
        // this->BC_Solid_Surface_Eedge_(U, fld, r, ngh);
    };
    if (!is_Lunar_interior_resis)
    {
        RegisterPhysical_("E_xi", "Solid_Surface", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_eta", "Solid_Surface", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_zeta", "Solid_Surface", Eedge_1stlayer_to_zero);
    }

    auto Eedge_boundandghost_to_constant = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Farfield_Eedge_set_zerocurl(U, fld, r, ngh);
    };
    RegisterPhysical_("E_xi", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("E_eta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("E_zeta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_xi", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_eta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_zeta", "Farfield", Eedge_boundandghost_to_constant);
    // ------------------------------------------
    // Pole
    //   E_xi / E_eta 保留你原来的 Pole 处理；
    //   E_zeta 显式保持 copy。
    // ------------------------------------------
    auto E_xi_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 1)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh); // xi norm
        else if (dir == 2)
            this->BC_Pole_Eedge_Axis(U, fld, r, ngh); // eta norm, xi using Pole -u\timesB\codt dr
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto E_eta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 2)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh); // eta norm
        else if (dir == 1)
            this->BC_Pole_Eedge_Axis(U, fld, r, ngh); // xi norm, eta using Pole -u\timesB\codt dr
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto E_zeta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_xi_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 1)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh);
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_eta_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 2)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh);
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_zeta_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    RegisterPhysical_("E_xi", "Pole", E_xi_pole);
    RegisterPhysical_("E_eta", "Pole", E_eta_pole);
    RegisterPhysical_("E_zeta", "Pole", E_zeta_pole);
    RegisterPhysical_("Eres_xi", "Pole", Eres_xi_pole);
    RegisterPhysical_("Eres_eta", "Pole", Eres_eta_pole);
    RegisterPhysical_("Eres_zeta", "Pole", Eres_zeta_pole);

    // ------------------------------------------
    //   J_xi / E_eta / J_zeta
    // ------------------------------------------
    // Pole
    auto J_xi_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 1)
            this->BC_Pole_Jedge_RegulateK_Norm(U, fld, r, ngh); // xi norm
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    auto J_eta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 2)
            this->BC_Pole_Jedge_RegulateK_Norm(U, fld, r, ngh); // eta norm
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    auto J_zeta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    RegisterPhysical_("J_xi", "Pole", J_xi_pole);
    RegisterPhysical_("J_eta", "Pole", J_eta_pole);
    RegisterPhysical_("J_zeta", "Pole", J_zeta_pole);

}

void LunarBoundary::InstallDefaultGroups()
{
    BoundGroup gU;
    gU.name = "Ucell";
    gU.fields = {"U_H"};
    gU.do_physical = true;
    gU.do_halo = true;
    gU.halo_level = HaloLevel::Vertex;
    AddGroup(gU);

    BoundGroup gJ;
    gJ.name = "Jedge";
    gJ.fields = {"J_xi", "J_eta", "J_zeta"};
    gJ.do_physical = true;
    gJ.do_halo = true;
    gJ.halo_level = HaloLevel::Vertex;
    AddGroup(gJ);

    BoundGroup gE;
    gE.name = "Eedge";
    gE.fields = {"E_xi", "E_eta", "E_zeta"};
    gE.do_physical = true;
    gE.do_halo = true;
    gE.halo_level = HaloLevel::Vertex;
    AddGroup(gE);

    BoundGroup gEhall;
    gEhall.name = "Ehall";
    gEhall.fields = {"Ehall_xi", "Ehall_eta", "Ehall_zeta"};
    gEhall.do_physical = true;
    gEhall.do_halo = true;
    gEhall.halo_level = HaloLevel::Vertex;
    AddGroup(gEhall);

    BoundGroup gEres;
    gEres.name = "Eres";
    gEres.fields = {"Eres_xi", "Eres_eta", "Eres_zeta"};
    gEres.do_physical = true;
    gEres.do_halo = true;
    gEres.halo_level = HaloLevel::Vertex;
    AddGroup(gEres);

    BoundGroup gEres1form;
    gEres1form.name = "Eres1form";
    gEres1form.fields = {"Eres_xi", "Eres_eta", "Eres_zeta"};
    gEres1form.do_physical = true;
    gEres1form.do_halo = true;
    gEres1form.halo_level = HaloLevel::Vertex;
    AddGroup(gEres1form);

    BoundGroup gB;
    gB.name = "Bface";
    gB.fields = {"B_xi", "B_eta", "B_zeta"};
    gB.do_physical = true;
    gB.do_halo = true;
    gB.halo_level = HaloLevel::Vertex;
    AddGroup(gB);

    BoundGroup EfaceB;
    EfaceB.name = "Eface";
    EfaceB.fields = {"Eface_xi", "Eface_eta", "Eface_zeta"};
    EfaceB.do_physical = true;
    EfaceB.do_halo = true;
    EfaceB.halo_level = HaloLevel::Vertex;
    AddGroup(EfaceB);

    BoundGroup gBc;
    gBc.name = "B_cell";
    gBc.fields = {"B_cell", "Bind_cell"};
    gBc.do_physical = true;
    gBc.do_halo = true;
    gBc.halo_level = HaloLevel::Vertex;
    AddGroup(gBc);

    BoundGroup gJc;
    gJc.name = "J_cell";
    gJc.fields = {"J_cell"};
    gJc.do_physical = true;
    gJc.do_halo = true;
    gJc.halo_level = HaloLevel::Vertex;
    AddGroup(gJc);

    BoundGroup gUplus;
    gUplus.name = "Uplus";
    gUplus.fields = {"U_plus"};
    gUplus.do_physical = false;
    gUplus.do_halo = true;
    gUplus.halo_level = HaloLevel::Vertex;
    AddGroup(gUplus);

    BoundGroup gdE;
    gdE.name = "dE";
    gdE.fields = {"dE_xi", "dE_eta", "dE_zeta"};
    gdE.do_physical = true;
    gdE.do_halo = true;
    gdE.halo_level = HaloLevel::Vertex;
    AddGroup(gdE);

    BoundGroup gdB;
    gdB.name = "dB";
    gdB.fields = {"dB_xi", "dB_eta", "dB_zeta"};
    gdB.do_physical = true;
    gdB.do_halo = true;
    gdB.halo_level = HaloLevel::Vertex;
    AddGroup(gdB);

    BoundGroup gdJ;
    gdJ.name = "dJ";
    gdJ.fields = {"dJ_xi", "dJ_eta", "dJ_zeta"};
    gdJ.do_physical = true;
    gdJ.do_halo = true;
    gdJ.halo_level = HaloLevel::Vertex;
    AddGroup(gdJ);

    BoundGroup gdJcell;
    gdJcell.name = "dJcell";
    gdJcell.fields = {"dJ_cell"};
    gdJcell.do_physical = true;
    gdJcell.do_halo = true;
    gdJcell.halo_level = HaloLevel::Vertex;
    AddGroup(gdJcell);

    BoundGroup gBadd;
    gBadd.name = "Badd";
    gBadd.fields = {"Badd_xi", "Badd_eta", "Badd_zeta"};
    gBadd.do_physical = true;
    gBadd.do_halo = true;
    gBadd.halo_level = HaloLevel::Vertex;
    AddGroup(gBadd);
}

void LunarBoundary::Build(bool strict_check)
{
    if (!halo_)
        ERROR::Abort("Build: call Setup first");

    // build halo patterns once
    // halo_->build_registered_patterns();

    if (strict_check)
        bound_.CheckPhysicalHandlers(boundary_fields_);

    built_ = true;
}

void LunarBoundary::InitBCStateFromParam_()
{
    bc_state_.SetBackground(BuildLunarBackgroundState(par_));
}
