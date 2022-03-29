#ifndef __MPM_h__
#define __MPM_h__
#include "Common.h"
#include "Grid.h"
#include "Particles.h"

template <int d>
class GridBasis
{
    using VectorD = Vector<real, d>;
    using VectorDi = Vector<int, d>;
    real h;

public:
    const real BSPLINEEPSILON = 5e-10;

    void init(real _h)
    {
        h = _h;
    }

    real CubicBsplines(real x)
    {
        real y = abs(x);
        if (y < 1)
        {
            return y * y * (0.5 * y - 1) + 2. / 3.;
        }
        else if (y < 2)
        {
            return y * (y * (-y / 6 + 1) - 2) + 4. / 3.;
        }
        else
        {
            return 0;
        }
    }

    real D_CubicBsplines(real x)
    {
        real y = abs(x);
        if (y < 1)
        {
            return 1.5 * x * y - 2 * x;
        }
        else if (y < 2)
        {
            return -0.5 * x * y + 2 * x - 2 * x / y;
        }
        else
        {
            return 0;
        }
    }

    real DyadicCubicBsplines(VectorD &x, VectorDi &node)
    {
        VectorD x_incell = (x - node.template cast<real>() * h) / h;

        real dyadic_cb = 1;
        for (int j = 0; j < d; j++)
        {
            dyadic_cb *= CubicBsplines(x_incell[j]);
        }
        return dyadic_cb;
    }

    VectorD gridientDCB(VectorD &x, VectorDi &node)
    {
        VectorD x_incell = (x - node.template cast<real>() * h) / h;

        VectorD gridient_dyadic_cb = VectorD::Ones();
        for (int j = 0; j < d; j++)
        {
            for (int k = 0; k < d; k++)
            {
                if (k != j)
                {
                    gridient_dyadic_cb[j] *= CubicBsplines(x_incell[k]);
                }
                else
                {
                    gridient_dyadic_cb[j] *= D_CubicBsplines(x_incell[k]);
                }
            }
        }
        return gridient_dyadic_cb / h;
    }
};

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template <int d>
class MPM
{
    using VectorD = Vector<real, d>;
    using Vector2 = Vector<real, 2>;
    using Vector3 = Vector<real, 3>;

    using VectorDi = Vector<int, d>;
    using Vector2i = Vector<int, 2>;
    using Vector3i = Vector<int, 3>;

    using MatrixD = Matrix<real, d>;
    using Matrix2 = Matrix<real, 2>;
    using Matrix3 = Matrix<real, 3>;

private:
    // Constants defined here
    const real EPS_MASS = 1e-100;

    const real HARDENING = 20.0,
               YOUNGS_MODULUS = 1.4e5,
               POISSONS_RATIO = 0.2, //0.2 soft  0.499 fluid
               STICKY = 1.,
               CRIT_COMPRESS = 1 - 2.5e-2,
               CRIT_STRETCH = 1 + 7.5e-3;

    const real PIC = 0.05,
               FLIP = 0.95;

    const real LAMBDA = YOUNGS_MODULUS * POISSONS_RATIO / ((1 + POISSONS_RATIO) * (1 - 2 * POISSONS_RATIO)),
               MU = YOUNGS_MODULUS / (2 + 2 * POISSONS_RATIO);

    const VectorD GRAVITY = -200. * VectorD::Unit(1);

    const real CGEPSILON = 1e-10;

    const int CGMAXITERATION = 50;

    real THICKNESS;

    const bool enable_implicit = false;

public:
    // Grid quantities
    GridBasis<d> gridBasis;
    Grid<d> grid;
    Array<VectorD> force_grid_node;
    Array<VectorD> delta_force_grid_node;
    Array<VectorD> v_grid_node;      // v_{n+1}
    Array<VectorD> v_grid_node_old;  // v_n
    Array<VectorD> v_grid_node_next; // v_n

    // Array<VectorD> pos_grid_node;
    Array<real> mass_grid_node;
    Array<real> mass_inv_grid_node;
    Array<bool> active_nodes;
    // CG values
    Array<VectorD> r_cg;
    Array<VectorD> s_cg;
    Array<VectorD> p_cg;
    Array<VectorD> q_cg;
    Array<real> gamma_cg;
    Array<real> alpha_cg;
    Array<real> beta_cg;
    Array<bool> implicit_active_nodes;

    Particles<d> particles;

    int node_num;
    real unit_grid_vol;
    real unit_grid_vol_inv;
    int n;

    virtual void Initialize()
    {
        n = 64;
        VectorDi cell_counts = VectorDi::Ones() * n;
        real dx = (real)1. / n;
        THICKNESS = 0.5 * dx;
        gridBasis.init(dx);

        unit_grid_vol = 1;

        for (int i = 0; i < d; i++)
        {
            unit_grid_vol *= dx;
        }

        unit_grid_vol_inv = 1. / unit_grid_vol;

        VectorD domain_min = VectorD::Zero();

        // grid domain: [0, 1]*[0, 1]*[0, 1]
        grid.Initialize(cell_counts, dx, domain_min);
        node_num = grid.node_counts.prod();

        // pos_grid_node.resize(node_num, VectorD::Unit(0) * (real).01);
        mass_grid_node.resize(node_num, 0.);
        mass_inv_grid_node.resize(node_num, 0.);

        v_grid_node.resize(node_num, VectorD::Zero());
        v_grid_node_old.resize(node_num, VectorD::Zero());

        force_grid_node.resize(node_num, VectorD::Zero());
        active_nodes.resize(node_num, false);

        if (enable_implicit)
        {
            InitializeCG();
        }

        Rasterize();

// Compute particle volume
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {
            real den = 0;
            for (int i = 0; i < node_num; i++)
            {
                if (active_nodes[i])
                {
                    VectorDi node = Coord(i);
                    real w_ip = gridBasis.DyadicCubicBsplines(particles.X(p), node);
                    den += mass_grid_node[i] * w_ip;
                }
            }
            particles.D(p) = den * unit_grid_vol_inv;
            particles.Vol(p) = particles.M(p) / particles.D(p);
        }
    }

    void InitializeCG()
    {
        r_cg.resize(node_num, VectorD::Zero());
        s_cg.resize(node_num, VectorD::Zero());
        p_cg.resize(node_num, VectorD::Zero());
        q_cg.resize(node_num, VectorD::Zero());
        delta_force_grid_node.resize(node_num, VectorD::Zero());

        gamma_cg.resize(node_num, (real)0.);
        beta_cg.resize(node_num, (real)0.);
        alpha_cg.resize(node_num, (real)0.);

        v_grid_node_next.resize(node_num, VectorD::Zero());
        implicit_active_nodes.resize(node_num, false);
    }

    // P2G
    virtual void Rasterize()
    {
#pragma omp parallel for
        for (int i = 0; i < node_num; i++)
        {
            VectorDi node = Coord(i);

            active_nodes[i] = false;
            v_grid_node[i] = VectorD::Zero();
            mass_grid_node[i] = 0.;
            VectorD vel_node = VectorD::Zero();
            real mass_node = 0.;
            force_grid_node[i] = VectorD::Zero();

            for (int j = 0; j < particles.Size(); j++)
            {
                real w_ip = gridBasis.DyadicCubicBsplines(particles.X(j), node);
                if (w_ip > gridBasis.BSPLINEEPSILON)
                {
                    mass_node += particles.M(j) * w_ip;
                    vel_node += particles.M(j) * w_ip * particles.V(j);
                    active_nodes[i] = true;
                }
            }

            if (active_nodes[i])
            {
                mass_grid_node[i] = mass_node;
                mass_inv_grid_node[i] = 1. / (mass_node + EPS_MASS);
                v_grid_node[i] = vel_node * mass_inv_grid_node[i];
                v_grid_node_old[i] = v_grid_node[i];
            }
        }
    }

    virtual void Update_Grid_Force()
    {
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {

            MatrixD &FE = particles.F_E(p);
            MatrixD &FP = particles.F_P(p);
            VectorD &Xp = particles.X(p);
            real Volp = particles.Vol(p);

            real JP = FP.determinant();
            real JE = FE.determinant();

            real harden = exp(HARDENING * (1 - JP));
            real mu = MU * harden;
            real lambda = LAMBDA * harden;

            Eigen::JacobiSVD<MatrixD> svd(FE, Eigen::ComputeFullU | Eigen::ComputeFullV);
            MatrixD RE = svd.matrixU() * svd.matrixV().transpose();

            for (int i = 0; i < node_num; i++)
            {
                if (active_nodes[i])
                {
                    VectorDi node = Coord(i);
                    VectorD grad_wip = gridBasis.gridientDCB(Xp, node);
                    force_grid_node[i] -=
                        Volp * (2 * mu * (FE - RE) * FE.transpose() + lambda * (JE - 1) * JE * MatrixD::Identity()) *
                        grad_wip;
                }
            }
        }
    }

    //// For explicit method, this is the grid velocity at n+1-th time step
    virtual void Update_RHS_Grid_Vel(real dt)
    {
#pragma omp parallel for
        for (int i = 0; i < node_num; i++)
        {
            if (active_nodes[i])
            {
                v_grid_node[i] += dt * (force_grid_node[i] * mass_inv_grid_node[i] + GRAVITY);
            }
        }
    }

    //// TODO: modify to arbitary env objects (e.g. bunny.txt)
    virtual void Grid_Collision(real dt)
    {
#pragma omp parallel for
        for (int i = 0; i < node_num; i++)
        {
            if (active_nodes[i])
            {
                VectorD &v = v_grid_node[i];
                // VectorD pos_new = Pos(i) + v * dt;
                VectorD pos_new = Pos(i);

                for (int j = 0; j < d; j++)
                {
                    // the y = 0 plane is smooth
                    if (j == 1 && pos_new[j] < THICKNESS)
                    {
                        v[1] = std::max((float)0.0, v[1]);
                    }
                    // stick to other boundaries
                    else if (pos_new[j] < THICKNESS || pos_new[j] > 1 - THICKNESS)
                    {
                        v = VectorD::Zero();
                        active_nodes[i] = false;
                    }
                    // if (pos_new[j] < THICKNESS || pos_new[j] > 1 - THICKNESS) {
                    //     v[j] *= -STICKY;
                    //     // v[j] = std::max(0.0, v[j]);
                    //     //active_nodes[i] = false;
                    // }
                }
            }
        }
    }

    const inline Matrix2 cofactor(Matrix2 m) const
    {
        Matrix2 M;
        M << m(1, 1), -m(1, 0),
            -m(0, 1), m(0, 0);
        return M;
    }

    const inline Matrix3 cofactor(Matrix3 m) const
    {
        Matrix3 M;
        M << m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1), -m(1, 0) * m(2, 2) + m(2, 0) * m(1, 2), m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1),
            -m(0, 1) * m(2, 2) + m(0, 2) * m(2, 1), m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0), -m(0, 0) * m(2, 1) + m(0, 1) * m(2, 0),
            m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1), -m(0, 0) * m(1, 2) + m(0, 2) * m(1, 0), m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
        return M;
    }

    // TODO: complete the semi implicit solver to solve for v_{n+1}
    // This function can be merged into grid force update to speed up the code
    // 2D case
    void Update_Particle_A(int p, Vector2 &u, real dt)
    {
        // 1. calculate delta_FE
        auto &FE = particles.F_E(p);
        auto &FP = particles.F_P(p);
        auto &Xp = particles.X(p);
        Matrix2 delta_FE = Matrix2::Zero();
        real Volp = particles.Vol(p);

        real JP = FP.determinant();
        real JE = FE.determinant();

        real harden = exp(HARDENING * (1 - JP));
        real mu = MU * harden;
        real lambda = LAMBDA * harden;

        for (int i = 0; i < node_num; i++)
        {
            if (active_nodes[i])
            {
                auto node = Coord(i);
                auto grad_wip = gridBasis.gridientDCB(Xp, node);
                delta_FE += dt * u * grad_wip.transpose();
            }
        }
        delta_FE *= FE;

        // 2. calculate delta_RE
        // x = (R^T delta_F - delta_F^T R)[0][1]
        Eigen::JacobiSVD<Matrix2> svd(FE, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Matrix2 RE = svd.matrixU() * svd.matrixV().transpose();

        // Matrix2 tempM = RE.transpose()* delta_FE - delta_FE.transpose() * RE;

        real x = RE(0, 0) * delta_FE(0, 1) + RE(1, 0) * delta_FE(1, 1) -
                 delta_FE(0, 0) * RE(0, 1) - delta_FE(1, 0) * RE(1, 1);

        x /= (RE(0, 0) + RE(1, 1));

        // 3. calculate R (R^T delta_R) = delta_R

        Matrix2 delta_RE;
        delta_RE << -x * RE(0, 1), x * RE(0, 0),
            -x * RE(1, 1), x * RE(1, 0);

        // 4. calculate JF^-T = cofactor(F)
        auto cofactor_FE = cofactor(FE);

        // 5. calculate delta(JF^-T)
        // In 2d, this is equal to cofactor(delta_F) https://github.com/Azmisov/snow
        // In 3d, see https://berkeley.mintkit.net/cs284b-projects/mpm-snow/assets/files/docs.pdf
        auto delta_cofactor_FE = cofactor(delta_FE);

        // 6. calculate A_p
        auto Ap = 2 * mu * (delta_FE - delta_RE) + lambda * cofactor_FE * (cofactor_FE.cwiseProduct(delta_FE)).sum() + lambda * (JE - 1) * delta_cofactor_FE;
        particles.A(p) = Ap;
    }
    //// TODO: 3D case
    ////
    VectorD Delta_Grid_Force(int i, VectorD &u, real dt)
    {
        VectorD delta_f = VectorD::Zero();
        if (active_nodes[i])
        {
            auto node = Coord(i);
#pragma omp parallel for
            for (int p = 0; p < particles.Size(); p++)
            {
                VectorD Xp = particles.X(p);
                auto &FE = particles.F_E(p);
                real Volp = particles.Vol(p);
                auto grad_wip = gridBasis.gridientDCB(Xp, node);
                // if (grad_wip.norm() > gridBasis.BSPLINEEPSILON) {
                Update_Particle_A(p, u, dt);
                delta_f -= Volp * (particles.A(p) * (FE.transpose() * grad_wip));
                // }
            }
        }
        return delta_f;
    }

    //// Conjugate residual solver
    //// https://en.wikipedia.org/wiki/Conjugate_residual_method
    virtual void Semi_Implicit_Update_Grid_Vel(real dt, real implicit_ratio)
    {
        real res_sum = 0;
        Array<real> residuals(node_num, (real)0.);
// initialize conjugate residual solver
#pragma omp parallel for reduction(+ \
                                   : res_sum)
        for (int i = 0; i < node_num; i++)
        {
            implicit_active_nodes[i] = active_nodes[i];
            if (implicit_active_nodes[i])
            {
                VectorD u = v_grid_node[i];
                v_grid_node_next[i] = u;
                r_cg[i] = dt * mass_inv_grid_node[i] * implicit_ratio * Delta_Grid_Force(i, u, dt);
                s_cg[i] = r_cg[i] - dt * mass_inv_grid_node[i] * implicit_ratio * Delta_Grid_Force(i, r_cg[i], dt);
                p_cg[i] = r_cg[i];
                q_cg[i] = s_cg[i];

                gamma_cg[i] = r_cg[i].dot(s_cg[i]);
                alpha_cg[i] = gamma_cg[i] / (q_cg[i].dot(q_cg[i]) + EPS_MASS);
                residuals[i] = alpha_cg[i] * alpha_cg[i] * p_cg[i].dot(p_cg[i]);
            }
            res_sum += residuals[i];
        }

        // CG solver
        int k = 0;
        for (; k < CGMAXITERATION && res_sum > CGEPSILON; k++)
        {
            res_sum = 0;
#pragma omp parallel for reduction(+ \
                                   : res_sum)
            for (int i = 0; i < node_num; i++)
            {
                if (implicit_active_nodes[i])
                {
                    alpha_cg[i] = gamma_cg[i] / (q_cg[i].dot(q_cg[i]) + EPS_MASS);
                    residuals[i] = alpha_cg[i] * alpha_cg[i] * p_cg[i].dot(p_cg[i]);
                    v_grid_node_next[i] += alpha_cg[i] * p_cg[i];
                    r_cg[i] -= alpha_cg[i] * q_cg[i];
                    s_cg[i] = r_cg[i] - dt * mass_inv_grid_node[i] * implicit_ratio * Delta_Grid_Force(i, r_cg[i], dt);
                    beta_cg[i] = r_cg[i].dot(s_cg[i]) / (gamma_cg[i] + EPS_MASS);
                    gamma_cg[i] *= beta_cg[i];
                    p_cg[i] = beta_cg[i] * p_cg[i] + r_cg[i];
                    q_cg[i] = beta_cg[i] * q_cg[i] + s_cg[i];
                }
                res_sum += residuals[i];
                if (!(res_sum > 1 || isnan(res_sum)))
                {
                    v_grid_node[i] = v_grid_node_next[i];
                }
                else
                {
                    residuals[i] = 0;
                    implicit_active_nodes[i] = false;
                }
            }
        }

        std::cout << k << ": " << res_sum << std::endl;
    }

    template <typename T>
    inline T clamp(const T &a, const T &min, const T &max)
    {
        if (a < min)
            return min;
        if (a > max)
            return max;
        return a;
    }

    virtual void Update_Deformation_Gradient(real dt)
    {
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {
            MatrixD &FE = particles.F_E(p);
            MatrixD &FP = particles.F_P(p);
            MatrixD F_total;
            VectorD &X = particles.X(p);

            MatrixD grad_vel_p = MatrixD::Zero();
            for (int i = 0; i < node_num; i++)
            {
                VectorDi node = Coord(i);
                grad_vel_p += v_grid_node[i] * gridBasis.gridientDCB(X, node).transpose();
            }

            FE += dt * grad_vel_p * FE;
            F_total = FE * FP;
            Eigen::JacobiSVD<MatrixD> svd(FE, Eigen::ComputeFullU | Eigen::ComputeFullV);

            MatrixD U = svd.matrixU();
            MatrixD Sigma = svd.singularValues().asDiagonal();
            MatrixD V = svd.matrixV();
            MatrixD Sigma_Inv = MatrixD::Zero();
            // Do clamp
            for (int k = 0; k < d; k++)
            {
                Sigma(k, k) = clamp(Sigma(k, k), CRIT_COMPRESS, CRIT_STRETCH);
                Sigma_Inv(k, k) = 1. / Sigma(k, k);
            }
            FE = U * Sigma * V.transpose();
            FP = V * Sigma_Inv * U.transpose() * F_total;
        }
    }

    virtual void Update_Particle_Velocity()
    {
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {
            VectorD pic = VectorD::Zero();
            VectorD flip = particles.V(p);
            VectorD xp = particles.X(p);
            for (int i = 0; i < node_num; i++)
            {
                if (active_nodes[i])
                {
                    VectorDi node = Coord(i);
                    real w_ip = gridBasis.DyadicCubicBsplines(xp, node);
                    pic += v_grid_node[i] * w_ip;
                    flip += (v_grid_node[i] - v_grid_node_old[i]) * w_ip;
                }
            }
            particles.V(p) = PIC * pic + FLIP * flip;
        }
    }

    // Particles collision
    virtual void Particle_Collision()
    {
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {
            VectorD &v = particles.V(p);
            VectorD &pos_new = particles.X(p);
            for (int j = 0; j < d; j++)
            {
                // the y = 0 plane is smooth
                // if (j == 1 && pos_new[j] < THICKNESS)
                // {
                //     v[1] = std::max(0.0, v[1]);
                // }
                // // stick to other boundaries
                // else if (pos_new[j] < THICKNESS || pos_new[j] > 1 - THICKNESS)
                // {
                //     v = VectorD::Zero();
                // }
                if (pos_new[j] < THICKNESS || pos_new[j] > 1 - THICKNESS)
                {
                    v[j] *= -STICKY;
                    // v[j] = std::max(0.0, v[j]);
                    // active_nodes[i] = false;
                }
            }
        }
    }

    virtual void Update_Particle_Position(real dt)
    {
#pragma omp parallel for
        for (int p = 0; p < particles.Size(); p++)
        {
            particles.X(p) += dt * particles.V(p);
        }
    }

    virtual void Advance(const real dt)
    {
        Update_Grid_Force();
        Update_RHS_Grid_Vel(dt);
        Grid_Collision(dt);
        // if (enable_implicit) {
        //     Semi_Implicit_Update_Grid_Vel(dt, 1.);
        // }
        Update_Deformation_Gradient(dt);
        Update_Particle_Velocity();
        Particle_Collision();
        Update_Particle_Position(dt);
        Rasterize();
    }

    ////Helper functions
protected:
    ////return the node index given its coordinate
    int Idx(const VectorDi &node_coord) const
    {
        return grid.Node_Index(node_coord);
    }

    ////return the coordinate given its index
    VectorDi Coord(const int node_index) const
    {
        return grid.Node_Coord(node_index);
    }

    ////return the node position given its index
    VectorD Pos(const int node_index) const
    {
        return grid.Node(node_index);
    }
};

#endif
