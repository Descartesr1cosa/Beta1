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
        double c_x = -std::sqrt(1.0 - c_y * c_y - c_z * c_z);
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

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &UNa = fld->field(fid.fid_U_Na, iblock);
            if (!UNa.is_allocated())
                continue;
            const Int3 &sub = UNa.get_lo();
            const Int3 &sup = UNa.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < UNa.descriptor().ncomp; ll++)
                            UNa(i, j, k, ll) = qinfs[ll];
                    }
        }

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &Ub = fld->field(fid.fid_U_b, iblock);
            if (!Ub.is_allocated())
                continue;
            const Int3 &sub = Ub.get_lo();
            const Int3 &sup = Ub.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < Ub.descriptor().ncomp; ll++)
                            Ub(i, j, k, ll) = 0.0; // induced Magnetic Field
                    }
        }

        par->AddParam("Physic_Time", 0.0);

        Initial_Badd_And_Bcell(fld, fid); // Badd Bcell

        Initial_Na_And_PhotoRate_(fld, fid); // Na Photo_ratec
    }

    // 构造初始外加磁场：IMF + 多个局部偶极子
    void Initial_Badd_And_Bcell(Field *fld, const SolverFields &fid)
    {

        const int nblock = fld->num_blocks();

        for (int ib = 0; ib < nblock; ++ib)
        {
            FieldBlock &Badd = fld->field(fid.fid_Badd, ib);
            if (!Badd.is_allocated())
                continue;

            Block &blk = fld->grd->grids(ib); // cell-center coords via dual_x/y/z

            const Int3 lo = Badd.get_lo();
            const Int3 hi = Badd.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // 1) IMF
                        double bx = B_imf[0];
                        double by = B_imf[1];
                        double bz = B_imf[2];

                        // 2) Local dipoles (pseudo)
                        const double x = blk.dual_x(i, j, k);
                        const double y = blk.dual_y(i, j, k);
                        const double z = blk.dual_z(i, j, k);

                        // for (auto& dip : dipoles_) {
                        //    add_dipole_field(dip, x,y,z, bx,by,bz);
                        // }

                        Badd(i, j, k, 0) = bx;
                        Badd(i, j, k, 1) = by;
                        Badd(i, j, k, 2) = bz;
                    }
        }

        // Calc B_cell
        for (int ib = 0; ib < nblock; ++ib)
        {
            FieldBlock &Badd = fld->field(fid.fid_Badd, ib);
            FieldBlock &Binduce = fld->field(fid.fid_U_b, ib);
            FieldBlock &Bcell = fld->field(fid.fid_Bcell, ib);
            if (!Bcell.is_allocated())
                continue;

            const Int3 lo = Bcell.get_lo();
            const Int3 hi = Bcell.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        Bcell(i, j, k, 0) = 0.0;
                        Bcell(i, j, k, 1) = 0.0;
                        Bcell(i, j, k, 2) = 0.0;
                    }

            if (Binduce.is_allocated())
            {
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            Bcell(i, j, k, 0) += Binduce(i, j, k, 0);
                            Bcell(i, j, k, 1) += Binduce(i, j, k, 1);
                            Bcell(i, j, k, 2) += Binduce(i, j, k, 2);
                        }
            }

            if (Badd.is_allocated())
            {
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            Bcell(i, j, k, 0) += Badd(i, j, k, 0);
                            Bcell(i, j, k, 1) += Badd(i, j, k, 1);
                            Bcell(i, j, k, 2) += Badd(i, j, k, 2);
                        }
            }
        }
    }

    // 构造中性 Na 原子密度与光致电离率
    void Initial_Na_And_PhotoRate_(Field *fld, const SolverFields &fid)
    {
        Param *par = fld->par;

        // ---- Fortran common1.f: sk1, sk2 ----
        // 建议在参数文件中提供 sk1/sk2；若还没有，就先在 ReadParam 后 AddParam 默认值
        const double sk1 = 5e-5; // par->GetDou("sk1"); // e.g. 5e-5
        const double sk2 = 1e-5; // par->GetDou("sk2"); // e.g. 1e-5

        // ---- Fortran common1.f: sl8 = R_M (m) ----
        // 用 REF.L_ref 作为水星半径（m）
        const double sl8 = par->GetDou_List("REF").data["L_ref"];

        const double pi = std::acos(-1.0);

        // ---- Na model constants (from 03_initiation.for) ----
        const double n_miv0 = 7.84e6;
        const double H_miv = 431.0e3;

        const double n_td0 = 8.86e9;
        const double H_td = 100.0e3;
        const double H_theta_td = 15.0 * pi / 180.0;

        const double n_psd0 = 4.06e10;
        const double H_psd = 232.0e3;
        const double H_theta_psd = 20.0 * pi / 180.0;
        const double theta_psd_n = 50.0 * pi / 180.0;
        const double theta_psd_s = -50.0 * pi / 180.0;

        const double n_sp0 = 5.67e6;
        const double H_sp = 748.0e3;
        const double H_theta_sp = 15.0 * pi / 180.0;
        const double H_theta_sp_night = 10.0 * pi / 180.0;
        const double theta_sp_day = 80.0 * pi / 180.0;
        const double theta_sp_night = 15.0 * pi / 180.0;

        auto gauss2 = [](double d, double H) -> double
        {
            if (H <= 0.0)
                return 0.0;
            const double a = d / H;
            return std::exp(-(a * a));
        };

        const int nblock = fld->num_blocks();

        for (int ib = 0; ib < nblock; ++ib)
        {
            FieldBlock &Na = fld->field(fid.fid_Na, ib);       // (Cell,1) neutral Na density
            FieldBlock &Photo = fld->field(fid.fid_Photo, ib); // (Cell,1) photo-ionization rate

            // 只对分配了这些场的 block（通常是 Fluid）计算
            if (!Na.is_allocated() || !Photo.is_allocated())
                continue;

            // cell-center coordinates (assumed nondimensional in R_M units)
            Block &blk = fld->grd->grids(ib);

            const Int3 lo = Na.get_lo();
            const Int3 hi = Na.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double x = blk.dual_x(i, j, k);
                        const double y = blk.dual_y(i, j, k);
                        const double z = blk.dual_z(i, j, k);

                        // hxz0 = r/R_M  (dimensionless)
                        const double hxz0 = std::sqrt(x * x + y * y + z * z);

                        // inside/on surface: set to zero (Fortran: hxz0 <= 1)
                        if (!(hxz0 > 1.0))
                        {
                            Na(i, j, k, 0) = 0.0;
                            Photo(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // physical radius r_norm (m)
                        const double r_norm = sl8 * hxz0;

                        // r_factor = (R_M/r)^2 = (sl8/r_norm)^2 = 1/hxz0^2
                        const double r_factor = 1.0 / (hxz0 * hxz0);

                        // theta = atan2(sqrt(y^2+z^2), x)  (0..pi)
                        const double r_perp = std::sqrt(y * y + z * z);
                        const double theta = std::atan2(r_perp, x);

                        // ------------------------
                        // 1) Neutral Na: roenum_neu
                        // ------------------------
                        const double dr = (r_norm - sl8); // (m), altitude above surface in meters

                        // MIV
                        const double n_miv = n_miv0 * std::exp(-dr / H_miv);

                        // TD
                        const double n_td = n_td0 * std::exp(-dr / H_td) * gauss2(std::fabs(theta), H_theta_td);

                        // PSD (north + south)
                        const double radial_psd = n_psd0 * std::exp(-dr / H_psd);
                        const double n_psd = radial_psd *
                                             (gauss2(std::fabs(theta - theta_psd_n), H_theta_psd) + gauss2(std::fabs(theta - theta_psd_s), H_theta_psd));

                        // SP day (two lobes: +/- theta_sp_day)
                        const double radial_sp = n_sp0 * std::exp(-dr / H_sp);
                        const double n_sp_day = radial_sp *
                                                (gauss2(std::fabs(theta - theta_sp_day), H_theta_sp) + gauss2(std::fabs(theta + theta_sp_day), H_theta_sp));

                        // SP night (only x < 0)
                        const double n_sp_night = (x < 0.0)
                                                      ? (radial_sp * gauss2(std::fabs(theta - theta_sp_night), H_theta_sp_night))
                                                      : 0.0;

                        // roenum_neu in m^-3
                        double roenum_m3 = r_factor * (n_miv + n_td + n_psd + n_sp_day + n_sp_night);
                        if (roenum_m3 < 0.0)
                            roenum_m3 = 0.0;

                        // convert to cm^-3 (Fortran: *1e-6)
                        const double roenum_cm3 = roenum_m3 * 1.0e-6;
                        Na(i, j, k, 0) = roenum_cm3;

                        // ------------------------
                        // 2) Photo-ionization rate: qmm
                        // ------------------------
                        // Fortran's cosxz = |x| / sqrt(x^2+y^2+z^2) = |cos(theta)|
                        const double cosxz = std::fabs(x) / hxz0; // in [0,1]

                        double qmm = 0.0;
                        if (x > 0.0)
                        {
                            // day side
                            qmm = (((sk1 - sk2) * cosxz) + sk2) * roenum_cm3;
                        }
                        else
                        {
                            // night side
                            qmm = sk2 * (1.0 - cosxz) * roenum_cm3;
                        }

                        if (qmm < 0.0)
                            qmm = 0.0;
                        Photo(i, j, k, 0) = qmm;
                    }
        }
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
