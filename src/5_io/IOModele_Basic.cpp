#include "5_io/IOModule.h"

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void IOModule::Setup(Param *par, Grid *grd, Field *fld)
{
    par_ = par;
    grd_ = grd;
    fld_ = fld;
    if (!par_ || !fld_ || !grd_)
        Fail_("[IOModule] Setup: null par/fld/grd");

    std::system("mkdir ./DATA");

    // 获取bin文件的输出路径
    const int myid = par_->GetInt("myid");
    auto rank4 = [](int id) -> std::string
    {
        char buf[8];
        std::snprintf(buf, sizeof(buf), "%04d", id);
        return std::string(buf);
    };
    restart_path_ = "./DATA/flow_field" + rank4(myid) + ".bin";
    tecplot_path_ = "./DATA/flow_field" + rank4(myid) + ".plt";
}

void IOModule::Fail_(const std::string &msg)
{
    std::fprintf(stderr, "%s\n", msg.c_str());
    std::abort();
}