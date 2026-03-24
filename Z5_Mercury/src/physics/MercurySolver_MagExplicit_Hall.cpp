#include "MercurySolver.h"
void MercurySolver::AddHallEdgeEMF_()
{
    // 0) 清零 Ehall_edge / Ehall_face
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto zero_scalar_edge = [&](FieldBlock &F)
        {
            if (!F.is_allocated())
                return;
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        F(i, j, k, 0) = 0.0;
        };

        auto zero_vec_face = [&](FieldBlock &F)
        {
            if (!F.is_allocated())
                return;
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 3; ++m)
                            F(i, j, k, m) = 0.0;
        };

        zero_scalar_edge(fld_->field(fid_.fid_Ehall.xi, ib));
        zero_scalar_edge(fld_->field(fid_.fid_Ehall.eta, ib));
        zero_scalar_edge(fld_->field(fid_.fid_Ehall.zeta, ib));

        zero_vec_face(fld_->field(fid_.fid_Eface.xi, ib));
        zero_vec_face(fld_->field(fid_.fid_Eface.eta, ib));
        zero_vec_face(fld_->field(fid_.fid_Eface.zeta, ib));
    }

    // 2) face 上做 Hall-Rusanov
    BuildHallFaceEMF_Rusanov_();

    // 3) sync Hall face field
    mercury_bound_.Sync("Eface");

    // 4) face -> edge
    AssembleEdgeEMF_FromFaceE_Hall_();

    // 5) sync edge Hall EMF
    mercury_bound_.Sync("Ehall");

    // 6) 加到总 E 上
}

void MercurySolver::BuildHallFaceEMF_Rusanov_()
{
    constexpr double eps = 1e-14;
    const double Cwh = 0.5;
    const double hall_coef_abs = std::abs(hall_coef);

    auto norm3 = [](double x, double y, double z) -> double
    {
        return std::sqrt(x * x + y * y + z * z);
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Bc = fld_->field(fid_.fid_Bcell, ib); // total B on cells
        auto &Binduce = fld_->field(fid_.fid_Bindcell, ib);
        auto &Jc = fld_->field(fid_.fid_Jcell, ib);

        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);

        auto &Efxi = fld_->field(fid_.fid_Eface.xi, ib);
        auto &Efet = fld_->field(fid_.fid_Eface.eta, ib);
        auto &Efze = fld_->field(fid_.fid_Eface.zeta, ib);

        auto &JDxi = fld_->field(fid_.fid_metric.xi, ib);
        auto &JDet = fld_->field(fid_.fid_metric.eta, ib);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, ib);

        auto &dlst_xi = fld_->field("dlstar_xi", ib);
        auto &dlst_et = fld_->field("dlstar_eta", ib);
        auto &dlst_ze = fld_->field("dlstar_zeta", ib);

        if (!Bc.is_allocated() || !Binduce.is_allocated() || !Jc.is_allocated() ||
            !UH.is_allocated() || !UNa.is_allocated() ||
            !Efxi.is_allocated() || !Efet.is_allocated() || !Efze.is_allocated())
        {
            continue;
        }

        auto &buf = hall_face_scratch_[ib];
        auto &Ehc = buf.Ehc;        // Ehc(i,j,k,comp)
        auto &beta_hall = buf.beta; // beta_hall(i,j,k)

        const Int3 clo = Bc.get_lo();
        const Int3 chi = Bc.get_hi();

        // ============================================================
        // 1) cell 上预计算 Ehc = alpha * (J x B), beta_hall
        // ============================================================
        for (int i = clo.i; i < chi.i; ++i)
            for (int j = clo.j; j < chi.j; ++j)
                for (int k = clo.k; k < chi.k; ++k)
                {
                    double num[3];
                    Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0), num);
                    const double ne = num[2];
                    const double alpha = hall_coef / ne;

                    const double Jx = Jc(i, j, k, 0);
                    const double Jy = Jc(i, j, k, 1);
                    const double Jz = Jc(i, j, k, 2);

                    const double Bx = Bc(i, j, k, 0);
                    const double By = Bc(i, j, k, 1);
                    const double Bz = Bc(i, j, k, 2);

                    Ehc(i, j, k, 0) = alpha * (Jy * Bz - Jz * By);
                    Ehc(i, j, k, 1) = alpha * (Jz * Bx - Jx * Bz);
                    Ehc(i, j, k, 2) = alpha * (Jx * By - Jy * Bx);

                    beta_hall(i, j, k) = hall_coef_abs * norm3(Bx, By, Bz) / ne;
                }

        // ============================================================
        // 2) xi-face
        // ============================================================
        {
            Int3 lo = Efxi.inner_lo();
            Int3 hi = Efxi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int iL = i - 1;
                        const int iR = i;

                        double ELx = Ehc(iL, j, k, 0);
                        double ELy = Ehc(iL, j, k, 1);
                        double ELz = Ehc(iL, j, k, 2);

                        double ERx = Ehc(iR, j, k, 0);
                        double ERy = Ehc(iR, j, k, 1);
                        double ERz = Ehc(iR, j, k, 2);

                        const double Sx = JDxi(i, j, k, 0);
                        const double Sy = JDxi(i, j, k, 1);
                        const double Sz = JDxi(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_xi(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(iL, j, k), beta_hall(iR, j, k)) / h_n;

                        const double BLx = Binduce(iL, j, k, 0);
                        const double BLy = Binduce(iL, j, k, 1);
                        const double BLz = Binduce(iL, j, k, 2);

                        const double BRx = Binduce(iR, j, k, 0);
                        const double BRy = Binduce(iR, j, k, 1);
                        const double BRz = Binduce(iR, j, k, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efxi(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efxi(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efxi(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }

        // ============================================================
        // 3) eta-face
        // ============================================================
        {
            Int3 lo = Efet.inner_lo();
            Int3 hi = Efet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int jL = j - 1;
                        const int jR = j;

                        double ELx = Ehc(i, jL, k, 0);
                        double ELy = Ehc(i, jL, k, 1);
                        double ELz = Ehc(i, jL, k, 2);

                        double ERx = Ehc(i, jR, k, 0);
                        double ERy = Ehc(i, jR, k, 1);
                        double ERz = Ehc(i, jR, k, 2);

                        const double Sx = JDet(i, j, k, 0);
                        const double Sy = JDet(i, j, k, 1);
                        const double Sz = JDet(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_et(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(i, jL, k), beta_hall(i, jR, k)) / h_n;

                        const double BLx = Binduce(i, jL, k, 0);
                        const double BLy = Binduce(i, jL, k, 1);
                        const double BLz = Binduce(i, jL, k, 2);

                        const double BRx = Binduce(i, jR, k, 0);
                        const double BRy = Binduce(i, jR, k, 1);
                        const double BRz = Binduce(i, jR, k, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efet(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efet(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efet(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }

        // ============================================================
        // 4) zeta-face
        // ============================================================
        {
            Int3 lo = Efze.inner_lo();
            Int3 hi = Efze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int kL = k - 1;
                        const int kR = k;

                        double ELx = Ehc(i, j, kL, 0);
                        double ELy = Ehc(i, j, kL, 1);
                        double ELz = Ehc(i, j, kL, 2);

                        double ERx = Ehc(i, j, kR, 0);
                        double ERy = Ehc(i, j, kR, 1);
                        double ERz = Ehc(i, j, kR, 2);

                        const double Sx = JDze(i, j, k, 0);
                        const double Sy = JDze(i, j, k, 1);
                        const double Sz = JDze(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_ze(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(i, j, kL), beta_hall(i, j, kR)) / h_n;

                        const double BLx = Binduce(i, j, kL, 0);
                        const double BLy = Binduce(i, j, kL, 1);
                        const double BLz = Binduce(i, j, kL, 2);

                        const double BRx = Binduce(i, j, kR, 0);
                        const double BRy = Binduce(i, j, kR, 1);
                        const double BRz = Binduce(i, j, kR, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efze(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efze(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efze(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }
    }
}

void MercurySolver::AssembleEdgeEMF_FromFaceE_Hall_()
{
    Int3 sub, sup;

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        double3D &x = fld_->grd->grids(iblk).x;
        double3D &y = fld_->grd->grids(iblk).y;
        double3D &z = fld_->grd->grids(iblk).z;

        Vec3 E, dr;

        // -----------------------------------------
        // Exi(edge) from face-eta and face-zeta
        // -----------------------------------------
        {
            auto &Exi = fld_->field(fid_.fid_Ehall.xi, iblk);
            auto &Ef_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_eta(i, j, k, 0) + Ef_eta(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i, j - 1, k, 0));
                        E.vec[1] = 0.25 * (Ef_eta(i, j, k, 1) + Ef_eta(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i, j - 1, k, 1));
                        E.vec[2] = 0.25 * (Ef_eta(i, j, k, 2) + Ef_eta(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i, j - 1, k, 2));

                        dr = {
                            x(i + 1, j, k) - x(i, j, k),
                            y(i + 1, j, k) - y(i, j, k),
                            z(i + 1, j, k) - z(i, j, k)};

                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Eeta(edge) from face-xi and face-zeta
        // -----------------------------------------
        {
            auto &Eeta = fld_->field(fid_.fid_Ehall.eta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i - 1, j, k, 2));

                        dr = {
                            x(i, j + 1, k) - x(i, j, k),
                            y(i, j + 1, k) - y(i, j, k),
                            z(i, j + 1, k) - z(i, j, k)};

                        Eeta(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Ezeta(edge) from face-xi and face-eta
        // -----------------------------------------
        {
            auto &Eze = fld_->field(fid_.fid_Ehall.zeta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_et = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Eze.inner_lo();
            sup = Eze.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j - 1, k, 0) + Ef_et(i, j, k, 0) + Ef_et(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j - 1, k, 1) + Ef_et(i, j, k, 1) + Ef_et(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j - 1, k, 2) + Ef_et(i, j, k, 2) + Ef_et(i - 1, j, k, 2));

                        dr = {
                            x(i, j, k + 1) - x(i, j, k),
                            y(i, j, k + 1) - y(i, j, k),
                            z(i, j, k + 1) - z(i, j, k)};

                        Eze(i, j, k, 0) = E * dr;
                    }
        }
    }
}
