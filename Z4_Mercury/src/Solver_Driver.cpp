#include <cmath>
#include <algorithm>
#include <iostream>

#include "MercurySolver.h"

void MercurySolver::Advance()
{
    RunData &run = io_.Run();

    // 初始先算一遍派生量（用于输出/诊断）
    PrepareStep();
    calc_Bcell();
    calc_PV();

    // step=0 也允许输出一次
    control_.UpdateSwitches(run);
    UpdateControlAndOutput();

    while (!control_.if_stop)
    {
        StepOnce();
    }
}

bool MercurySolver::StepOnce()
{
    RunData &run = io_.Run();

    PrepareStep();
    calc_Bcell();
    calc_PV();
    Compute_Timestep();

    // 这里先只推进“时间与步数”，场量暂不更新（下一步再加 Time_Advance）
    run.dt = dt;
    run.time += dt;
    run.step += 1;

    control_.UpdateSwitches(run);
    return UpdateControlAndOutput();
}

void MercurySolver::PrepareStep()
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