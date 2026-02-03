// SolverFields.h
#pragma once

#include <array>
#include <cstdio>
#include <cstdlib>

// Core
#include "3_field/2_MPCNS_Field.h"

// 用于缓存“字段 id 的三分量”，这里存 int
struct IdTriplet
{
    int xi = -1;
    int eta = -1;
    int zeta = -1;

    bool all_valid() const { return (xi >= 0) && (eta >= 0) && (zeta >= 0); }

    void require_all(const char *what) const
    {
        if (!all_valid())
        {
            std::fprintf(stderr, "[IdTriplet] %s not fully bound (xi=%d eta=%d zeta=%d)\n",
                         what ? what : "(null)", xi, eta, zeta);
            std::abort();
        }
    }

    const int &at(int dir) const
    {
        switch (dir)
        {
        case 0:
            return xi;
        case 1:
            return eta;
        case 2:
            return zeta;
        default:
            std::fprintf(stderr, "[IdTriplet] invalid dir=%d (expect 0/1/2)\n", dir);
            std::abort();
        }
    }
};

struct SolverFields
{
    Field *fld = nullptr;

    // ---- geometry ----
    int fid_Jac = -1;
    IdTriplet fid_metric; // (xi,eta,zeta) <- (JDxi,JDet,JDze)
    IdTriplet fid_pinvGT; // (pinvGT_xi, pinvGT_eta, pinvGT_zeta)  ncomp=9
    IdTriplet fid_pinvAT; // (pinvAT_xi, pinvAT_eta, pinvAT_zeta)  ncomp=9

    // ---- field ids ----
    int fid_U_H = -1;
    int fid_U_Na = -1;
    IdTriplet fid_B;     // (xi,eta,zeta) <- (B_xi,B_eta,B_zeta)
    IdTriplet fid_E;     // (xi,eta,zeta) <- (E_xi,E_eta,E_zeta)
    IdTriplet fid_Eface; // (xi,eta,zeta) <- (E_xi,E_eta,E_zeta)
    IdTriplet fid_J;     // (xi,eta,zeta) <- (J_xi,J_eta,J_zeta)

    // ---- auxiliary ----
    int fid_PV_H = -1;
    int fid_PV_Na = -1;
    int fid_Bcell = -1;
    IdTriplet fid_Badd;
    int fid_Na = -1;
    int fid_Photo = -1;

    // ---- fluid flux and face E ----
    IdTriplet fid_F; // F_xi/F_eta/F_zeta

    // -------- buffers / time advance ----------
    int fid_RHS_H = -1;
    int fid_RHS_Na = -1;
    IdTriplet fid_RHS_b;
    int fid_U_plus = -1;
    // int fid_old_U = -1;
    int fid_divB = -1;
    // IdTriplet fid_old_Bface; // old_B_xi/eta/zeta

    // int fid_RHS_U = -1;          // RHS (cell,5)
    // IdTriplet fid_RHS_Bface;     // RHS_xi/eta/zeta (face,1)
    // IdTriplet fid_RHShall_Bface; // RHShall_xi/eta/zeta

    void Init(Field *fld_in)
    {
        fld = fld_in;

        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        // ---- geometry ----
        fid_Jac = fld->field_id("Jac");
        fid_metric.xi = fld->field_id("JDxi");
        fid_metric.eta = fld->field_id("JDet");
        fid_metric.zeta = fld->field_id("JDze");
        fid_pinvGT.xi = fld->field_id("pinvGT_xi");
        fid_pinvGT.eta = fld->field_id("pinvGT_eta");
        fid_pinvGT.zeta = fld->field_id("pinvGT_zeta");
        fid_pinvAT.xi = fld->field_id("pinvAT_xi");
        fid_pinvAT.eta = fld->field_id("pinvAT_eta");
        fid_pinvAT.zeta = fld->field_id("pinvAT_zeta");

        // ---- field ids ----
        fid_U_H = fld->field_id("U_H");
        fid_U_Na = fld->field_id("U_Na");
        fid_B.xi = fld->field_id("B_xi");
        fid_B.eta = fld->field_id("B_eta");
        fid_B.zeta = fld->field_id("B_zeta");

        fid_E.xi = fld->field_id("E_xi");
        fid_E.eta = fld->field_id("E_eta");
        fid_E.zeta = fld->field_id("E_zeta");

        fid_Eface.xi = fld->field_id("Eface_xi");
        fid_Eface.eta = fld->field_id("Eface_eta");
        fid_Eface.zeta = fld->field_id("Eface_zeta");

        fid_J.xi = fld->field_id("J_xi");
        fid_J.eta = fld->field_id("J_eta");
        fid_J.zeta = fld->field_id("J_zeta");

        // ---- auxiliary ----
        fid_PV_H = fld->field_id("PV_H");
        fid_PV_Na = fld->field_id("PV_Na");
        fid_Bcell = fld->field_id("B_cell");
        fid_Badd.xi = fld->field_id("Badd_xi");
        fid_Badd.eta = fld->field_id("Badd_eta");
        fid_Badd.zeta = fld->field_id("Badd_zeta");
        fid_Na = fld->field_id("Na");
        fid_Photo = fld->field_id("Photo_rate");

        fid_F.xi = fld->field_id("F_xi");
        fid_F.eta = fld->field_id("F_eta");
        fid_F.zeta = fld->field_id("F_zeta");

        // -------- buffers / time advance ----------
        fid_RHS_H = fld->field_id("RHS_H");
        fid_RHS_Na = fld->field_id("RHS_Na");
        fid_RHS_b.xi = fld->field_id("RHS_B_xi");
        fid_RHS_b.eta = fld->field_id("RHS_B_eta");
        fid_RHS_b.zeta = fld->field_id("RHS_B_zeta");
        fid_U_plus = fld->field_id("U_plus");
        // fid_old_U = fld->field_id("old_U_");
        fid_divB = fld->field_id("divB");
        // fid_old_Bface.xi = fld->field_id("old_B_xi");
        // fid_old_Bface.eta = fld->field_id("old_B_eta");
        // fid_old_Bface.zeta = fld->field_id("old_B_zeta");

        // fid_RHS_U = fld->field_id("RHS");
        // fid_RHS_Bface.xi = fld->field_id("RHS_xi");
        // fid_RHS_Bface.eta = fld->field_id("RHS_eta");
        // fid_RHS_Bface.zeta = fld->field_id("RHS_zeta");

        // -------- Check --------
        Validate();
    }

    void Validate() const
    {
        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        auto require_id = [](int id, const char *name)
        {
            if (id < 0)
            {
                std::fprintf(stderr, "[SolverFields] missing field id: %s\n", name ? name : "(null)");
                std::abort();
            }
        };

        // ---- geometry ----
        require_id(fid_Jac, "Jac");
        fid_metric.require_all("metric(JDxi/JDet/JDze)");
        fid_pinvGT.require_all("pinvGT(edge)");
        fid_pinvAT.require_all("pinvAT(edge)");

        // ---- field ids ----
        require_id(fid_U_H, "U_H");
        require_id(fid_U_Na, "U_Na");
        fid_B.require_all("B_xi/B_eta/B_zeta");
        fid_E.require_all("E_xi/E_eta/E_zeta");
        fid_Eface.require_all("Eface_xi/Eface_eta/Eface_zeta");
        fid_J.require_all("J_xi/J_eta/J_zeta");

        // ---- primary / auxiliary ----
        require_id(fid_PV_H, "PV_H");
        require_id(fid_PV_Na, "PV_Na");
        require_id(fid_Bcell, "B_cell");
        fid_Badd.require_all("Badd_xi/Badd_eta/Badd_zeta");
        require_id(fid_Na, "Na");
        require_id(fid_Photo, "Photo_rate");
        fid_F.require_all("Flux(F_xi/F_eta/F_zeta)");

        // ---- buffers / time advance ----
        require_id(fid_RHS_H, "RHS_H");
        require_id(fid_RHS_Na, "RHS_Na");
        fid_RHS_b.require_all("Flux(RHS_B_xi/RHS_B_eta/RHS_B_zeta)");
        require_id(fid_U_plus, "U_plus");
        // require_id(fid_old_U, "old_U_");
        require_id(fid_divB, "divB");
        // fid_old_Bface.require_all("old_B_face(old_B_xi/eta/zeta)");

        // require_id(fid_RHS_U, "RHS");
        // fid_RHS_Bface.require_all("RHS_B_face(RHS_xi/eta/zeta)");
    }
};