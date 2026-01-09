#include <cmath>
#include <algorithm>
#include <iostream>

#include "MercurySolver.h"

void MercurySolver::Advance()
{
    // 初始先算一遍派生量（用于输出/诊断）
    Boundary_Condition();
    calc_Bcell();
    calc_PV();
    calc_Uplus();

    // step=0 也允许输出一次
    control_.UpdateSwitches(*run_data_);
    UpdateControlAndOutput();

    while (!control_.if_stop)
    {
        StepOnce();
    }
}

bool MercurySolver::StepOnce()
{
    // Calculate time step From CFL
    Compute_Timestep();

    // Time Advance
    Time_Advance();

    // Record and Update Runtime DATA
    {
        run_data_->dt = dt;
        run_data_->time += dt;
        run_data_->step += 1;
    }

    // Add Boundary Condition And Prepare for Next Step
    {
        Boundary_Condition();
        calc_Bcell();
        calc_PV();
        calc_Uplus();
    }

    // When Stop/Output/Print
    control_.UpdateSwitches(*run_data_);
    return UpdateControlAndOutput();
}

void MercurySolver::Boundary_Condition()
{
    // 1) 物理边界
    bound_.ApplyPhysical(std::vector<std::string>{"U_H", "U_Na", "U_b"});

    // 2) halo
    {
        std::string fn = "U_H";
        halo_->data_trans_1DCorner(fn);
        fn = "U_Na";
        halo_->data_trans_1DCorner(fn);
        fn = "U_b";
        halo_->data_trans_1DCorner(fn);
    }

    // 3) coupling（填 buffer -> 写入 dst ghost）
    {
        std::string src = "Solid", dst = "Fluid";
        halo_->coupling_trans_1DCorner(src, dst);
        bound_.ApplyCouplingPair(src, dst);

        src = "Fluid";
        dst = "Solid";
        halo_->coupling_trans_1DCorner(src, dst);
        bound_.ApplyCouplingPair(src, dst);
    }
}

bool MercurySolver::UpdateControlAndOutput()
{
    RunData &run = io_.Run();

    if (control_.if_outres)
    {
        if (par_->GetInt("myid") == 0)
        {
            std::cout << "[Mercury] step=" << run.step
                      << " time=" << run.time
                      << " dt=" << run.dt << "\n";
        }
    }

    if (control_.if_outfile)
    {
        io_.WriteTecplotBinFile(run.step, run.time);
        io_.WriteRestartBinFile(run.step, run.time);
        io_.WriteRunDataFile();
    }

    return control_.if_stop;
}