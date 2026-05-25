#pragma once

#include "3_field/2_MPCNS_Field.h"

#include "00_Lunar_Const.h"
#include "0_BackgroundState.h"
#include "0_SolverFields.h"
#include "operators/Dipole.h"

class Lunar_Initial
{
public:
    Lunar_Initial() = default;
    ~Lunar_Initial() = default;

    void Initialization(Field *fld, SolverFields &fid)
    {
        build_background_state(fld->par);

        if (!fld->par->GetBoo("continue_calc"))
            Common_Initial(fld, fid);

        Initial_Badd(fld, fid);

        if (fld->par->GetInt("myid") == 0)
        {
            std::cout << "********************Finish the Initial Process! !*********************\n\n\n";
            std::cout << "====================================="
                      << "===============================================================\n";
            std::cout << "\t \t Multi-Physics Coupling Numerical Simulation Solver Begins NOW!\n";
            std::cout << "====================================="
                      << "===============================================================\n\n\n";
        }
    }

private:
    double qinf[5]{};
    double B_imf[3]{};
    DipoleField dipB;

    void build_background_state(Param *par)
    {
        const LunarBackgroundState state = BuildLunarBackgroundState(par);
        for (int i = 0; i < 5; ++i)
            qinf[i] = state.qinf[i];
        for (int i = 0; i < 3; ++i)
            B_imf[i] = state.B_imf[i];
        dipB.load_from_param(par);
    }

    void Common_Initial(Field *fld, SolverFields &fid)
    {
        if (fld->par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        for (int32_t iblock = 0; iblock < fld->num_blocks(); ++iblock)
        {
            FieldBlock &U = fld->field(fid.fid_U_H, iblock);
            if (!U.is_allocated())
                continue;
            const Int3 &lo = U.get_lo();
            const Int3 &hi = U.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int32_t m = 0; m < U.descriptor().ncomp; ++m)
                            U(i, j, k, m) = qinf[m];
        }

        for (int32_t iblock = 0; iblock < fld->num_blocks(); ++iblock)
        {
            FieldBlock &Bxi = fld->field(fid.fid_B.xi, iblock);
            FieldBlock &Beta = fld->field(fid.fid_B.eta, iblock);
            FieldBlock &Bzeta = fld->field(fid.fid_B.zeta, iblock);
            if (!Bxi.is_allocated() || !Beta.is_allocated() || !Bzeta.is_allocated())
                continue;
            zero_field(Bxi);
            zero_field(Beta);
            zero_field(Bzeta);
        }
    }

    static void zero_field(FieldBlock &field)
    {
        const Int3 &lo = field.get_lo();
        const Int3 &hi = field.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    for (int32_t m = 0; m < field.descriptor().ncomp; ++m)
                        field(i, j, k, m) = 0.0;
    }

    void Initial_Badd(Field *fld, const SolverFields &fid)
    {
        dipB.Build_Badd_FaceFlux(fld->grd, fld, fld->par,
                                 fid.fid_Badd.xi, fid.fid_Badd.eta, fid.fid_Badd.zeta);
    }
};
