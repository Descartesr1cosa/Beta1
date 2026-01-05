#pragma once

#include "3_field/2_MPCNS_Field.h"

#include "00_Mercury_Const.h"

class Mercury_Initial
{
public:
    double qinf[5], q_pv_inf[5], B_imf[3]; // SW: conservative quanities, primitive quantities, IMF
    double qinfs[5], q_pv_infs[5];         // Na+: seed initial state
    void Initialization(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;
        par->AddParam("Nstep", 0);

        build_background_state(par);

        if (par->GetBoo("continue_calc"))
        {
            Read_bin_Initial(fld);
        }
        else
        {
            Common_Initial(fld, fid);
        }

        if (par->GetInt("myid") == 0)
        {
            std::cout << "********************Finish the Initial Process! !*********************\n\n\n";
            std::cout << "====================================="
                      << "===============================================================\n";
            std::cout << "\t \t Multi-Physics Coupling Numerical Simulation Solver Begins NOW!\n";
            std::cout << "====================================="
                      << "===============================================================\n\n\n";
        }
    }

    void build_background_state(Param *par)
    {
        // ---- 常数 ----
        double gamma = par->GetDou_List("constant").data["gamma"];
        double NA = par->GetDou_List("constant").data["NA"];
        double R_uni = par->GetDou_List("constant").data["R_uni"];
        double q_e = par->GetDou_List("constant").data["q_e"];
        double k_Boltz = R_uni / NA;
        double mu_mag = par->GetDou_List("constant").data["mu_mag"];

        List<double> ini = par->GetDou_List("INITIAL");
        List<double> ref = par->GetDou_List("REF");

        // ---- 背景磁场（物理单位）----
        double Bx_phy = ini.data["Bx"];
        double By_phy = ini.data["By"];
        double Bz_phy = ini.data["Bz"];

        double B_ref = ref.data["B_ref"];

        // 无量纲IMF磁场
        double Bx = Bx_phy / B_ref;
        double By = By_phy / B_ref;
        double Bz = Bz_phy / B_ref;

        // 流体参考量
        double L_ref = ref.data["L_ref"];
        double U_ref = ref.data["U"];
        double n_ref = ref.data["n"];
        double T_ref = ref.data["T"];
        double Molecular_mass_ref = ref.data["Molecular_mass"];
        double rho_ref = Molecular_mass_ref / NA * n_ref;

        // ---- 无量纲原始量 ----
        double rho0 = (ini.data["n"] / n_ref) * (ini.data["Molecular_mass"] / Molecular_mass_ref);

        double c_y = ini.data["c_y"];
        double c_z = ini.data["c_z"];
        double c_x = std::sqrt(1.0 - c_y * c_y - c_z * c_z);
        double u0 = c_x * ini.data["U"] / U_ref;
        double v0 = c_y * ini.data["U"] / U_ref;
        double w0 = c_z * ini.data["U"] / U_ref;

        double p_ini = ini.data["n"] * k_Boltz * ini.data["T"];
        double p0 = p_ini / (rho_ref * U_ref * U_ref);

        double T0 = ini.data["T"] / T_ref;

        // SW state
        q_pv_inf[0] = u0;
        q_pv_inf[1] = v0;
        q_pv_inf[2] = w0;
        q_pv_inf[3] = p0;
        q_pv_inf[4] = T0;

        qinf[0] = rho0;
        qinf[1] = rho0 * u0;
        qinf[2] = rho0 * v0;
        qinf[3] = rho0 * w0;
        qinf[4] = 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) // 动能
                  + p0 / (gamma - 1.0);                      // 内能

        // IMF
        B_imf[0] = Bx;
        B_imf[1] = By;
        B_imf[2] = Bz;

        // qinfs[5], q_pv_infs[5];// Na+: seed initial state

        q_pv_infs[0] = 0.0;       // 静止，速度为零
        q_pv_infs[1] = 0.0;       // 静止，速度为零
        q_pv_infs[2] = 0.0;       // 静止，速度为零
        q_pv_infs[3] = p0 * 1E-6; // 极低压力背景场
        q_pv_infs[4] = T0;        // 温度与来流一致

        qinfs[0] = rho_small * rho0;
        qinfs[1] = 0.0;
        qinfs[2] = 0.0;
        qinfs[3] = 0.0;
        qinfs[4] = q_pv_infs[3] / (gamma - 1.0) + 0.5 * qinfs[0] * (q_pv_infs[0] * q_pv_infs[0] + q_pv_infs[1] * q_pv_infs[1] + q_pv_infs[2] * q_pv_infs[2]);
    }

public:
    Mercury_Initial() {};
    ~Mercury_Initial() = default;

private:
    void Compute_Fluid_Energy(Field *fld_)
    {
        double gamma = fld_->par->GetDou_List("constant").data["gamma"];

        const int nblock = fld_->num_blocks();

        for (int ib = 0; ib < nblock; ++ib)
        {
            auto &U = fld_->field("U_", ib);
            auto &PV = fld_->field("PV_", ib);

            Int3 lo = U.get_lo();
            Int3 hi = U.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        double rho = U(i, j, k, 0);
                        double u = U(i, j, k, 1) / rho;
                        double v = U(i, j, k, 2) / rho;
                        double w = U(i, j, k, 3) / rho;
                        double kin = 0.5 * rho * (u * u + v * v + w * w);
                        double E_inner = PV(i, j, k, 3) / (gamma - 1.0);

                        U(i, j, k, 4) = E_inner + kin; // 这里只存储流体能量
                    }
        }
    }

    void Common_Initial(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        int ngg = par->GetInt("ngg");

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &UH = fld->field(fid.fid_U_H, iblock);
            if (!UH.is_allocated())
                continue;
            const Int3 &sub = UH.get_lo();
            const Int3 &sup = UH.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < UH.descriptor().ncomp; ll++)
                            UH(i, j, k, ll) = qinf[ll];
                    }
        }

        par->AddParam("Physic_Time", 0.0);
    }

    void Read_bin_Initial(Field *fld)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        //----------------------------------------------------------------------
        // 打开二进制文件  read the "flow_field   *.bin" binary files
        std::string _my_id_s;
        int my_id = par->GetInt("myid");
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }
        std::ifstream inFile_bin("./DATA/flow_field" + _my_id_s + ".bin",
                                 std::ios_base::in | std::ios::binary);
        if (!inFile_bin)
        {
            std::cerr << "Error: cannot open ./DATA/flow_field"
                      << _my_id_s << ".bin" << std::endl;
            return;
        }

        //----------------------------------------------------------------------
        // 一些小工具：读名字（带长度）
        auto bin_read_name = [](std::ifstream &file) -> std::string
        {
            int32_t temp_i;
            //  1、读入名称长度
            file.read((char *)&temp_i, sizeof(int32_t));
            //  2、分配空间并读入名称，存入string后释放
            char *temp_char = new char[temp_i];
            file.read(temp_char, temp_i);
            std::string physical_name;
            physical_name.reserve(temp_i);
            physical_name = temp_char;
            delete[] temp_char;
            return physical_name;
        };

        // 读 int 参数： (length + name + value)
        auto bin_read_int_param = [&](std::ifstream &file)
        {
            std::string pname = bin_read_name(file);
            int value = 0;
            file.read(reinterpret_cast<char *>(&value), sizeof(int));
            // 根据你的 Param 接口写回去，这里假定有 SetInt / SetDou
            par->AddParam(pname, value);
        };

        // 读 double 参数： (length + name + value)
        auto bin_read_double_param = [&](std::ifstream &file)
        {
            std::string pname = bin_read_name(file);
            double value = 0.0;
            file.read(reinterpret_cast<char *>(&value), sizeof(double));
            par->AddParam(pname, value);
        };

        // kind(int) -> StaggerLocation 的反向映射
        auto map_int2loc = [](int kind) -> StaggerLocation
        {
            switch (kind)
            {
            case 0:
                return StaggerLocation::Cell;
            case 1:
                return StaggerLocation::FaceXi;
            case 2:
                return StaggerLocation::FaceEt;
            case 3:
                return StaggerLocation::FaceZe;
            case 4:
                return StaggerLocation::EdgeXi;
            case 5:
                return StaggerLocation::EdgeEt;
            case 6:
                return StaggerLocation::EdgeZe;
            case 7:
                return StaggerLocation::Node;
            default:
                return StaggerLocation::Cell; // fallback
            }
        };

        //----------------------------------------------------------------------
        // 1. 读入 int 类型记录变量
        int32_t num_int = 0;
        inFile_bin.read(reinterpret_cast<char *>(&num_int), sizeof(int32_t));
        for (int32_t n = 0; n < num_int; ++n)
        {
            bin_read_int_param(inFile_bin);
        }

        //----------------------------------------------------------------------
        // 2. 读入 double 类型记录变量
        int32_t num_double = 0;
        inFile_bin.read(reinterpret_cast<char *>(&num_double), sizeof(int32_t));
        for (int32_t n = 0; n < num_double; ++n)
        {
            bin_read_double_param(inFile_bin);
        }

        //----------------------------------------------------------------------
        // 3. 读入 data 信息
        // 3.1 网格块数
        int32_t nblock_file = 0;
        inFile_bin.read(reinterpret_cast<char *>(&nblock_file), sizeof(int32_t));

        int32_t nblock = fld->num_blocks();
        if (nblock_file != nblock)
        {
            std::cerr << "Error: nblock in file (" << nblock_file
                      << ") != fld->num_blocks() (" << nblock << ")\n";
            exit(-1);
        }

        // 3.2 物理场个数（文件中实际输出的）
        int32_t n_field_file = 0;
        inFile_bin.read(reinterpret_cast<char *>(&n_field_file), sizeof(int32_t));

        // 注意：文件里只包含你写出时没有被 skip 的物理场，
        // 我们按名字匹配到当前 Field 中的 fid 再回填数据。
        for (int32_t ifld = 0; ifld < n_field_file; ++ifld)
        {
            // 3.2.1 读 FieldDescriptor：name + kind + ncomp + nghost
            std::string fname = bin_read_name(inFile_bin);

            int32_t kind = 0;
            int32_t ncomp_file = 0;
            int32_t nghost_file = 0;
            inFile_bin.read(reinterpret_cast<char *>(&kind), sizeof(int32_t));
            inFile_bin.read(reinterpret_cast<char *>(&ncomp_file), sizeof(int32_t));
            inFile_bin.read(reinterpret_cast<char *>(&nghost_file), sizeof(int32_t));

            // 在当前 Field 中查这个名字的 field ID
            int fid = -1;
            try
            {
                fid = fld->field_id(fname); // Field::field_id(name):contentReference[oaicite:1]{index=1}
            }
            catch (const std::out_of_range &)
            {
                std::cerr << "Error: field name \"" << fname
                          << "\" not found in current Field. Skip this field.\n";
                // 如果找不到，只能把这块数据从文件中 “读掉丢弃”
                // 以保证后面能对齐。
                // 这里简单实现一个丢弃过程：
                for (int ib = 0; ib < nblock_file; ++ib)
                {
                    int32_t lo[3], hi[3];
                    inFile_bin.read(reinterpret_cast<char *>(&lo[0]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&lo[1]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&lo[2]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[0]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[1]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[2]), sizeof(int32_t));

                    const int Ni = hi[0] - lo[0];
                    const int Nj = hi[1] - lo[1];
                    const int Nk = hi[2] - lo[2];
                    const std::size_t nval =
                        static_cast<std::size_t>(Ni) *
                        static_cast<std::size_t>(Nj) *
                        static_cast<std::size_t>(Nk) *
                        static_cast<std::size_t>(ncomp_file);

                    // 直接跳过这些 double
                    inFile_bin.seekg(static_cast<std::streamoff>(
                                         nval * sizeof(double)),
                                     std::ios_base::cur);
                }
                continue; // 进入下一个 field
            }

            const FieldDescriptor &desc = fld->descriptor(fid); //: contentReference[oaicite:2]{index=2}

            // 一些健壮性检查（可选）
            if (desc.ncomp != ncomp_file)
            {
                std::cerr << "Error: field \"" << fname << "\" ncomp mismatch: file="
                          << ncomp_file << ", current=" << desc.ncomp << std::endl;
                exit(-1);
            }
            if (desc.nghost != nghost_file)
            {
                std::cerr << "Error: field \"" << fname << "\" nghost mismatch: file="
                          << nghost_file << ", current=" << desc.nghost << std::endl;
                exit(-1);
            }
            // 如果你希望严格一致，也可以这里直接 return 或 throw

            // 3.2.2 对每个 block 读 lo/hi 和数据
            for (int ib = 0; ib < nblock_file; ++ib)
            {
                int32_t lo[3], hi[3];
                inFile_bin.read(reinterpret_cast<char *>(&lo[0]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&lo[1]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&lo[2]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[0]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[1]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[2]), sizeof(int32_t));

                FieldBlock &fb = fld->field(fid, ib); //: contentReference[oaicite:3]{index=3}

                // 可选：检查 lo/hi 是否和当前 fb 一致
                // const Int3 &cur_lo = fb.get_lo();
                // const Int3 &cur_hi = fb.get_hi();

                // 实际数据读取：完全按写出时的循环顺序读回来
                for (int i = lo[0]; i < hi[0]; ++i)
                    for (int j = lo[1]; j < hi[1]; ++j)
                        for (int k = lo[2]; k < hi[2]; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                            {
                                double val = 0.0;
                                inFile_bin.read(reinterpret_cast<char *>(&val), sizeof(double));
                                fb(i, j, k, m) = val;
                            }
            }
        }

        inFile_bin.close();
    }
};
