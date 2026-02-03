
#include "MercurySolver.h"

void MercurySolver::AddIdealEdgeEMF_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Uplus = fld_->field(fid_.fid_U_plus, iblk);
        if (!Uplus.is_allocated())
            continue;
        auto &Bcell = fld_->field(fid_.fid_Bcell, iblk);
        //  三方向：通量 + ideal face EMF
        for (int dir = 0; dir < 3; ++dir)
        {
            auto &E_face = fld_->field(fid_.fid_Eface.at(dir), iblk);
            auto &B_face = fld_->field(fid_.fid_B.at(dir), iblk);
            auto &B_face_add = fld_->field(fid_.fid_Badd.at(dir), iblk);
            auto &metric = fld_->field(fid_.fid_metric.at(dir), iblk); // Xi_[iblk]/Eta_[iblk]/Zeta_[iblk]

            AssembleOneDirectionEMF_(iblk, dir, E_face, B_face, B_face_add, Bcell, metric, Uplus);
        }
    }

    AssembleEdgeEMF_FromFaceE_Ideal_();
}

//=========================================================================
// Eelectric face  to  edge
void MercurySolver::AssembleEdgeEMF_FromFaceE_Ideal_()
{
    Int3 sub, sup;
    // 插值计算电场E=-u\times B, 存储的E_xi eta zeta均为E\cdot dr的线积分量
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Bcell = fld_->field(fid_.fid_Bcell, iblk);
        Vec3 E, dr;
        double3D &x = fld_->grd->grids(iblk).x;
        double3D &y = fld_->grd->grids(iblk).y;
        double3D &z = fld_->grd->grids(iblk).z;

        {
            auto &Exi = fld_->field(fid_.fid_E.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
            sub = Exi.inner_lo();
            sup = Exi.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {

                        E.vec[0] = 0.25 * (E_face_eta(i, j, k, 0) + E_face_eta(i, j, k - 1, 0) + E_face_zeta(i, j, k, 0) + E_face_zeta(i, j - 1, k, 0));
                        E.vec[1] = 0.25 * (E_face_eta(i, j, k, 1) + E_face_eta(i, j, k - 1, 1) + E_face_zeta(i, j, k, 1) + E_face_zeta(i, j - 1, k, 1));
                        E.vec[2] = 0.25 * (E_face_eta(i, j, k, 2) + E_face_eta(i, j, k - 1, 2) + E_face_zeta(i, j, k, 2) + E_face_zeta(i, j - 1, k, 2));

                        dr = {x(i + 1, j, k) - x(i, j, k),
                              y(i + 1, j, k) - y(i, j, k),
                              z(i + 1, j, k) - z(i, j, k)};
                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        {
            auto &Eeta = fld_->field(fid_.fid_E.eta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        E.vec[0] = 0.25 * (E_face_xi(i, j, k, 0) + E_face_xi(i, j, k - 1, 0) + E_face_zeta(i, j, k, 0) + E_face_zeta(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (E_face_xi(i, j, k, 1) + E_face_xi(i, j, k - 1, 1) + E_face_zeta(i, j, k, 1) + E_face_zeta(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (E_face_xi(i, j, k, 2) + E_face_xi(i, j, k - 1, 2) + E_face_zeta(i, j, k, 2) + E_face_zeta(i - 1, j, k, 2));

                        dr = {x(i, j + 1, k) - x(i, j, k),
                              y(i, j + 1, k) - y(i, j, k),
                              z(i, j + 1, k) - z(i, j, k)};
                        Eeta(i, j, k, 0) = E * dr;
                    }
        }

        {
            auto &Ezeta = fld_->field(fid_.fid_E.zeta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            sub = Ezeta.inner_lo();
            sup = Ezeta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        E.vec[0] = 0.25 * (E_face_xi(i, j, k, 0) + E_face_xi(i, j - 1, k, 0) + E_face_eta(i, j, k, 0) + E_face_eta(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (E_face_xi(i, j, k, 1) + E_face_xi(i, j - 1, k, 1) + E_face_eta(i, j, k, 1) + E_face_eta(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (E_face_xi(i, j, k, 2) + E_face_xi(i, j - 1, k, 2) + E_face_eta(i, j, k, 2) + E_face_eta(i - 1, j, k, 2));

                        dr = {x(i, j, k + 1) - x(i, j, k),
                              y(i, j, k + 1) - y(i, j, k),
                              z(i, j, k + 1) - z(i, j, k)};
                        Ezeta(i, j, k, 0) = E * dr;
                    }
        }
    }
}

//=========================================================================
void MercurySolver::AssembleOneDirectionEMF_(
    int iblk,
    int dir,                // 0 xi, 1 eta, 2 zeta
    FieldBlock &E_face,     // E_face_xi/eta/zeta   (ncomp=3)
    FieldBlock &B_face,     // B_xi/eta/zeta        (ncomp=1)
    FieldBlock &B_face_add, // B_xi/eta/zeta add    (ncomp=1)
    FieldBlock &Bcell,
    FieldBlock &metricField, // Xi_/Eta_/Zeta_      (ncomp=3)
    FieldBlock &Uplus)
{
    Int3 sub = E_face.inner_lo();
    Int3 sup = E_face.inner_hi();

    double metric[3];
    double flux3[3]; // 注意：Reconstruction 输出 8 个，其中[5..7] 已被旋转成 E_face

    for (int i = sub.i; i < sup.i; ++i)
        for (int j = sub.j; j < sup.j; ++j)
            for (int k = sub.k; k < sup.k; ++k)
            {
                metric[0] = metricField(i, j, k, 0);
                metric[1] = metricField(i, j, k, 1);
                metric[2] = metricField(i, j, k, 2);

                ReconstructionEMF_(metric, dir, Uplus, Bcell,
                                   B_face(i, j, k, 0) + B_face_add(i, j, k, 0), iblk, i, j, k, flux3);

                // 2) Ideal MHD electric field on face
                E_face(i, j, k, 0) = flux3[0];
                E_face(i, j, k, 1) = flux3[1];
                E_face(i, j, k, 2) = flux3[2];
            }
}

void MercurySolver::ReconstructionEMF_(double *metric, int32_t direction,
                                       FieldBlock &Uplus, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k,
                                       double *out_flux)
{
    auto calc_Jac_radius_GCL = [&](double &out, double *pv, double *B, double *metric)
    {
        double u, v, w;
        double BB2, Bx, By, Bz;
        double K1, K2, K3;
        u = pv[0];
        v = pv[1];
        w = pv[2];

        Bx = B[0];
        By = B[1];
        Bz = B[2];

        BB2 = Bx * Bx + By * By + Bz * Bz;

        K1 = metric[0];
        K2 = metric[1];
        K3 = metric[2];

        double uvw = K1 * u + K2 * v + K3 * w;
        // double cc1 = sqrt((gamma_ * p / rho) * (K1 * K1 + K2 * K2 + K3 * K3));
        // double cc = sqrt((gamma_ * p / rho + BB2 / rho * inver_MA2) * (K1 * K1 + K2 * K2 + K3 * K3));
        // constexpr double C_hall_safe = 1.5;
        // double c_hall = C_hall_safe * ion_inertial_len * sqrt(BB2 / rho * inver_MA2 * (K1 * K1 + K2 * K2 + K3 * K3)); // ≈ Jac * v_A * d_i * |k|
        // out = fabs(uvw) + cc + c_hall;
        out = 1.5 * fabs(uvw); //+ cc;
        return;
    };

    // 注意这里的B长度为4，最后一个为B_Jac_nabla\xi eta zeta
    auto calc_Jac_Flux_GCL = [&](double *flux, double *pv, double *B, double *metric)
    {
        double k1, k2, k3; // GCL 这里为Jac *k1, Jac *k2, Jac *k3
        k1 = metric[0];
        k2 = metric[1];
        k3 = metric[2];

        double u, v, w, uvw;
        double Bx, By, Bz, B_Jac_nabla;

        u = pv[0];
        v = pv[1];
        w = pv[2];

        Bx = B[0];
        By = B[1];
        Bz = B[2];
        B_Jac_nabla = B[3]; // Bx * k1 + By * k2 + Bz * k3;
        uvw = k1 * u + k2 * v + k3 * w;

        flux[0] = uvw * Bx - B_Jac_nabla * u;
        flux[1] = uvw * By - B_Jac_nabla * v;
        flux[2] = uvw * Bz - B_Jac_nabla * w;
    };

    int i = index_i;
    int j = index_j;
    int k = index_k;

    double ppvvL[3], ppvvR[3], BL[4], BR[4];
    double radius[2];

    auto fill_state = [&](int ic, int jc, int kc, double *pv, double *B)
    {
        double u = Uplus(ic, jc, kc, 0);
        double v = Uplus(ic, jc, kc, 1);
        double w = Uplus(ic, jc, kc, 2);

        double Bx = B_cell(ic, jc, kc, 0); // including B_add
        double By = B_cell(ic, jc, kc, 1); // including B_add
        double Bz = B_cell(ic, jc, kc, 2); // including B_add

        double inner_product = Bx * metric[0] + By * metric[1] + Bz * metric[2];
        // double inner_product_add = B_add_x * metric[0] + B_add_y * metric[1] + B_add_z * metric[2];
        double inver_norm2 = 1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);

        double B_jac_total = B_jac_nabla; // 法向通量仅仅为induced部分

        pv[0] = u;
        pv[1] = v;
        pv[2] = w;

        B[0] = Bx;
        B[1] = By;
        B[2] = Bz;
        B[3] = B_jac_total;
    };

    if (direction == 0)
    {
        int iL = index_i - 1;
        int iR = index_i;
        fill_state(iL, j, k, ppvvL, BL);
        fill_state(iR, j, k, ppvvR, BR);
    }
    else if (direction == 1)
    {
        int jL = index_j - 1;
        int jR = index_j;
        fill_state(i, jL, k, ppvvL, BL);
        fill_state(i, jR, k, ppvvR, BR);
    }
    else
    { // direction == 2
        int kL = index_k - 1;
        int kR = index_k;
        fill_state(i, j, kL, ppvvL, BL);
        fill_state(i, j, kR, ppvvR, BR);
    }

    calc_Jac_radius_GCL(radius[0], ppvvL, BL, metric);
    calc_Jac_radius_GCL(radius[1], ppvvR, BR, metric);

    double radius_max = std::max(radius[0], radius[1]);

    double FL[8], FR[8];
    calc_Jac_Flux_GCL(FL, ppvvL, BL, metric);
    calc_Jac_Flux_GCL(FR, ppvvR, BR, metric);

    for (int m = 0; m < 3; ++m)
        out_flux[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * radius_max * (BR[m] - BL[m]);

    double Elec_flux[3] = {out_flux[0], out_flux[1], out_flux[2]};
    double norm2 = -1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);
    out_flux[0] = norm2 * (metric[1] * Elec_flux[2] - metric[2] * Elec_flux[1]); // Averaged Electric in Face xi eta zeta
    out_flux[1] = norm2 * (metric[2] * Elec_flux[0] - metric[0] * Elec_flux[2]); // Averaged Electric in Face xi eta zeta
    out_flux[2] = norm2 * (metric[0] * Elec_flux[1] - metric[1] * Elec_flux[0]); // Averaged Electric in Face xi eta zeta
}
