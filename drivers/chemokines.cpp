#include <finite_difference_problem.h>
#include <config_file.h>

#include <iostream>
#include <fstream>

using namespace mjrfd;

struct ChemokinesParams
{
    double p_u;
    double alpha;
    double beta;
    double gamma_u;
    double gamma_b;
    double D_su;
    double D_ju;
    double nu;
    double lambda;
};

class ChemokinesProblem1D : public FiniteDifferenceProblem
{
    enum class Variable : unsigned
    {
        c_u = 0,
        c_b = 1,
        c_s = 2,
        phi = 3
    };

    typedef Eigen::Triplet<double> T;

public:
    // dofs:
    // C_u - 0--n_node-1
    // C_b - n_node--2*n_node-1
    // C_s - 2*n_node--3*n_node-1
    // phi - 3*n_node--4*n_node-1
    ChemokinesProblem1D(const unsigned n_node, const double dt) :
        FiniteDifferenceProblem(4, n_node, dt),
        n_node_(n_node),
        dx_(1.0/(n_node-1)),
        cn_theta_(1.0),
        time_factor_(1.0),
        c_u_offset_(0),
        c_b_offset_(n_node_),
        c_s_offset_(2*n_node_),
        phi_offset_(3*n_node_)
    {
        std::cout << "n_node = " << n_node_ << '\n';
        std::cout << "n_dof  = " << n_dof_ << '\n';
        std::cout << "dx     = " << dx_ << '\n';
        std::cout << "dt     = " << dt_ << "\n\n";
    }

    ~ChemokinesProblem1D()
    {
    }

    void set_initial_conditions()
    {
        // set zero initial conditions
        for(unsigned i = 0; i < n_node_; ++i)
        {
            c_u(0,i) = 0.0;
            c_b(0,i) = 0.0;
            c_s(0,i) = 0.0;
            phi(0,i) = 0.0;
        }

        c_u(0,0) = 1.0;
        phi(0,n_node_-1) = 1.0;
    }

    ChemokinesParams p;

    const double c_u(const unsigned t, const unsigned i) const
    {
        return u(t, Variable::c_u, i);
    }

    const double c_b(const unsigned t, const unsigned i) const
    {
        return u(t, Variable::c_b, i);
    }

    const double c_s(const unsigned t, const unsigned i) const
    {
        return u(t, Variable::c_s, i);
    }

    const double phi(const unsigned t, const unsigned i) const
    {
        return u(t, Variable::phi, i);
    }

    double& c_u(const unsigned t, const unsigned i)
    {
        return u(t, Variable::c_u, i);
    }

    double& c_b(const unsigned t, const unsigned i)
    {
        return u(t, Variable::c_b, i);
    }

    double& c_s(const unsigned t, const unsigned i)
    {
        return u(t, Variable::c_s, i);
    }

    double& phi(const unsigned t, const unsigned i)
    {
        return u(t, Variable::phi, i);
    }

    const double d_c_u_dx(const unsigned t, const unsigned i) const
    {
        return d_u_dx_helper(t, Variable::c_u, i);
    }

    const double d_c_b_dx(const unsigned t, const unsigned i) const
    {
        return d_u_dx_helper(t, Variable::c_b, i);
    }

    const double d_c_s_dx(const unsigned t, const unsigned i) const
    {
        return d_u_dx_helper(t, Variable::c_s, i);
    }

    const double d_phi_dx(const unsigned t, const unsigned i) const
    {
        return d_u_dx_helper(t, Variable::phi, i);
    }

    void make_steady() override
    {
        Problem::make_steady();
        time_factor_ = 0.0;
        old_cn_theta_ = cn_theta_;
        cn_theta_ = 1.0;
    }

    void make_unsteady() override
    {
        Problem::make_unsteady();
        time_factor_ = 1.0;
        cn_theta_ = old_cn_theta_;
    }

    void output(std::ostream &out) const override
    {
        // Output headers
        out << "x c_u c_b c_s phi dc_udx dc_bdx dc_sdx dphidx\n";

        for(unsigned i = 0; i < n_node_; ++i)
        {
            // Output values
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' '
                << c_u(0,i) << ' '
                << c_b(0,i) << ' '
                << c_s(0,i) << ' '
                << phi(0,i) << ' ';

            // Output derivatives
            out << d_c_u_dx(0,i) << ' '
                << d_c_b_dx(0,i) << ' '
                << d_c_s_dx(0,i) << ' '
                << d_phi_dx(0,i) << '\n';
        }
    }

private:
    const unsigned n_node_;
    const double dx_;
    double cn_theta_;
    double time_factor_;
    double old_cn_theta_;

    const unsigned c_u_offset_;
    const unsigned c_b_offset_;
    const unsigned c_s_offset_;
    const unsigned phi_offset_;

    const double d_u_dx_helper(const unsigned t,
                               const Variable variable,
                               const unsigned i) const
    {
        if(i == 0)
        {
            // Left boundary
            return stencil_1_forward(t,variable,i)/dx_;
        }
        else if(i == (n_node_-1))
        {
            // Right boundary
            return stencil_1_backward(t,variable,i)/dx_;
        }
        else
        {
            // Sanity check (debug only)
            assert(i > 0 && i < (n_node_-1));

            // Bulk
            return stencil_1_central(t,variable,i)/dx_;
        }
    }

    const double upwind_d_u_helper(const unsigned t, const Variable variable, const unsigned i, const double v) const
    {
        // TODO rewrite as an if statement? Possibly avoid calculating fw+bw diffs every time?
        //return std::max(0.0,v)*upwind_d_u_backward_helper(t, variable, i) + std::min(0.0,v)*upwind_d_u_forward_helper(t, variable, i);

        if(v > 0.0)
        {
            return upwind_d_u_backward_helper(t, variable, i);
        }
        else
        {
            return upwind_d_u_forward_helper(t, variable, i);
        }
    }

    const double upwind_v_d_u_helper(const unsigned t, const Variable variable, const unsigned i, const double v) const
    {
        // TODO rewrite as an if statement? Possibly avoid calculating fw+bw diffs every time?
        //return std::max(0.0,v)*upwind_d_u_backward_helper(t, variable, i) + std::min(0.0,v)*upwind_d_u_forward_helper(t, variable, i);

        if(v > 0.0)
        {
            return v*upwind_d_u_backward_helper(t, variable, i);
        }
        else
        {
            return v*upwind_d_u_forward_helper(t, variable, i);
        }
    }

    const double upwind_d_u_backward_helper(const unsigned t, const Variable variable, const unsigned i) const
    {
        if(i == 0)
        {
            // Left boundary (forward difference)
            return stencil_1_forward(t, variable, i);
        }
        else if(i == 1)
        {
            // Left boundary+1 (central difference)
            return stencil_1_central(t, variable, i);
        }
        else
        {
            // Bulk and right boundary (backward difference)
            return stencil_1_backward(t, variable, i);
        }
    }

    const double upwind_d_u_forward_helper(const unsigned t, const Variable variable, const unsigned i) const
    {
        if(i == (n_node_-1))
        {
            // Right boundary (backward difference)
            return stencil_1_backward(t, variable, i);
        }
        else if(i == (n_node_-2))
        {
            // Right boundary-1 (central difference)
            return stencil_1_central(t, variable, i);
        }
        else
        {
            // Bulk and left boundary (forward difference)
            return stencil_1_forward(t, variable, i);
        }
    }

    const double upwind_v_d_c_u(const unsigned t, const unsigned i, const double v) const
    {
        return upwind_v_d_u_helper(t, Variable::c_u, i, v);
    }

    const double upwind_v_d_c_b(const unsigned t, const unsigned i, const double v) const
    {
        return upwind_v_d_u_helper(t, Variable::c_b, i, v);
    }

    const double upwind_v_d_c_s(const unsigned t, const unsigned i, const double v) const
    {
        return upwind_v_d_u_helper(t, Variable::c_s, i, v);
    }

    const double upwind_v_d_phi(const unsigned t, const unsigned i, const double v) const
    {
        return upwind_v_d_u_helper(t, Variable::phi, i, v);
    }

    const double upwind_d_phi(const unsigned t, const unsigned i, const double v) const
    {
        return upwind_d_u_helper(t, Variable::phi, i, v);
    }

    void calculate_residual()
    {
        // set all entries to zero
        // not neccessary here since all entries are set explicitly
        residual_.setZero();

        // LHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u
        residual_(c_u_offset_+0) += (c_u(0,0) - 1.0)*dx_*dx_;

        // C_b (full equation, not just BC)
        // time derivative
        residual_(c_b_offset_+0) += -time_factor_*(c_b(0,0) - c_b(1,0))*dx_*dx_/dt_;

        // binding/unbinding (Crank-Nicolson)
        residual_(c_b_offset_+0) += cn_theta_*(p.alpha*c_u(0,0) - p.beta*c_b(0,0) - p.gamma_b*c_b(0,0)*phi(0,0))*dx_*dx_;
        residual_(c_b_offset_+0) += (1.0-cn_theta_)*(p.alpha*c_u(1,0) - p.beta*c_b(1,0) - p.gamma_b*c_b(1,0)*phi(1,0))*dx_*dx_;

        // C_s
        residual_(c_s_offset_+0) += c_s(0,0)*dx_*dx_;

        // phi
        residual_(phi_offset_+0) += stencil_1_forward(0,Variable::phi,0)*dx_ - p.lambda*phi(0,0)*dx_*dx_;

        // Bulk equations
        // --------------------------------------------------------------------

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // C_u
            // time derivative
            residual_(c_u_offset_+i) += -time_factor_*(c_u(0,i) - c_u(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(c_u_offset_+i) += cn_theta_*stencil_2_central(0,Variable::c_u,i);
            residual_(c_u_offset_+i) += (1.0-cn_theta_)*stencil_2_central(1,Variable::c_u,i);

            // advection (Crank-Nicolson)
            residual_(c_u_offset_+i) += -cn_theta_*upwind_v_d_c_u(0,i,p.p_u)*dx_;
            residual_(c_u_offset_+i) += -(1.0-cn_theta_)*upwind_v_d_c_u(1,i,p.p_u)*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_u_offset_+i) += cn_theta_*(-p.alpha*c_u(0,i) + p.beta*c_b(0,i) - p.gamma_u*c_u(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_u_offset_+i) += (1.0-cn_theta_)*(-p.alpha*c_u(1,i) + p.beta*c_b(1,i) - p.gamma_u*c_u(1,i)*phi(1,i))*dx_*dx_;

            // C_b
            // time derivative
            residual_(c_b_offset_+i) += -time_factor_*(c_b(0,i) - c_b(1,i))*dx_*dx_/dt_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_b_offset_+i) += cn_theta_*(p.alpha*c_u(0,i) - p.beta*c_b(0,i) - p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_b_offset_+i) += (1.0-cn_theta_)*(p.alpha*c_u(1,i) - p.beta*c_b(1,i) - p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;

            // C_s
            // time derivative
            residual_(c_s_offset_+i) += -time_factor_*(c_s(0,i) - c_s(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(c_s_offset_+i) += cn_theta_*p.D_su*stencil_2_central(0,Variable::c_s,i);
            residual_(c_s_offset_+i) += (1.0-cn_theta_)*p.D_su*stencil_2_central(1,Variable::c_s,i);

            // advection (Crank-Nicolson)
            residual_(c_s_offset_+i) += -cn_theta_*upwind_v_d_c_s(0,i,p.p_u)*dx_;
            residual_(c_s_offset_+i) += -(1.0-cn_theta_)*upwind_v_d_c_s(1,i,p.p_u)*dx_;

            // binding/unbinding (Crank-Nicolson)
            residual_(c_s_offset_+i) += cn_theta_*(p.gamma_u*c_u(0,i)*phi(0,i) + p.gamma_b*c_b(0,i)*phi(0,i))*dx_*dx_;
            residual_(c_s_offset_+i) += (1.0-cn_theta_)*(p.gamma_u*c_u(1,i)*phi(1,i) + p.gamma_b*c_b(1,i)*phi(1,i))*dx_*dx_;

            // phi
            // time derivative
            residual_(phi_offset_+i) += -time_factor_*(phi(0,i) - phi(1,i))*dx_*dx_/dt_;

            // diffusion (Crank-Nicolson)
            residual_(phi_offset_+i) += cn_theta_*p.D_ju*stencil_2_central(0,Variable::phi,i);
            residual_(phi_offset_+i) += (1.0-cn_theta_)*p.D_ju*stencil_2_central(1,Variable::phi,i);

            // advection (Chemotaxis?) (Crank-Nicolson)
            residual_(phi_offset_+i) += -cn_theta_*p.nu*upwind_v_d_phi(0,i,stencil_1_central(0,Variable::c_b,i));
            residual_(phi_offset_+i) += -(1.0-cn_theta_)*p.nu*upwind_v_d_phi(1,i,stencil_1_central(1,Variable::c_b,i));

            residual_(phi_offset_+i) += -cn_theta_*p.nu*phi(0,i)*stencil_2_central(0,Variable::c_b,i);
            residual_(phi_offset_+i) += -(1.0-cn_theta_)*p.nu*phi(1,i)*stencil_2_central(1,Variable::c_b,i);
        }

        // RHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u
        residual_(c_u_offset_+n_node_-1) += (c_u(0,n_node_-1) - 0.0)*dx_*dx_;

        // C_b (full equation, not just BC)
        // time derivative
        residual_(c_b_offset_+n_node_-1) += -time_factor_*(c_b(0,n_node_-1) - c_b(1,n_node_-1))*dx_*dx_/dt_;

        // binding/unbinding (Crank-Nicolson)
        residual_(c_b_offset_+n_node_-1) += cn_theta_*(p.alpha*c_u(0,n_node_-1) - p.beta*c_b(0,n_node_-1) - p.gamma_b*c_b(0,n_node_-1)*phi(0,n_node_-1))*dx_*dx_;
        residual_(c_b_offset_+n_node_-1) += (1.0-cn_theta_)*(p.alpha*c_u(1,n_node_-1) - p.beta*c_b(1,n_node_-1) - p.gamma_b*c_b(1,n_node_-1)*phi(1,n_node_-1))*dx_*dx_;

        // C_s
        residual_(c_s_offset_+n_node_-1) += c_s(0,n_node_-1)*dx_*dx_;

        // phi
        residual_(phi_offset_+n_node_-1) += (phi(0,n_node_-1) - 1.0)*dx_*dx_;
    }

    void calculate_jacobian()
    {
        // set all entries to zero
        //jacobian_.setZero();

        std::vector<T> triplet_list;
        // TODO update the number of entries after completion of upwinding implementation
        //triplet_list.reserve(10 + 20*(n_node_-2) + 8);
        triplet_list.reserve(8 + 36*(n_node_-2) + 6);

        // LHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u / C_u
        triplet_list.push_back( T(c_u_offset_+0, c_u_offset_+0, 1.0*dx_*dx_) );

        // C_b / C_u
        triplet_list.push_back( T(c_b_offset_+0, c_u_offset_+0, cn_theta_*p.alpha*dx_*dx_) );

        // C_b / C_b
        triplet_list.push_back( T(c_b_offset_+0, c_b_offset_+0, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,0))*dx_*dx_) );

        // C_b / C_s
        // (no entries)

        // C_b / phi
        triplet_list.push_back( T(c_b_offset_+0, phi_offset_+0, -cn_theta_*(p.gamma_b*c_b(0,0))*dx_*dx_) );

        // C_s / C_s
        triplet_list.push_back( T(c_s_offset_+0, c_s_offset_+0, 1.0*dx_*dx_) );

        // phi / phi
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+0, -1.5*dx_ - p.lambda*dx_*dx_) );
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+1, 2.0*dx_) );
        triplet_list.push_back( T(phi_offset_+0, phi_offset_+2, -0.5*dx_) );

        // Bulk equations
        // --------------------------------------------------------------------

        for(unsigned i = 1; i < n_node_-1; ++i)
        {
            // C_u / C_u
            // time derivative
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i, -time_factor_*dx_*dx_/dt_) );

            // diffusion
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-1,  cn_theta_*1.0) );
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i,   -cn_theta_*2.0) );
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+1,  cn_theta_*1.0) );

            // advection (upwinding)
            if(p.p_u > 0)
            {
                // if peclet is positive
                if(i == 1)
                {
                    // use central difference at left boundary+1
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-1,  cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+1, -cn_theta_*p.p_u*0.5*dx_) );
                }
                else
                {
                    // use backward difference away from left boundary
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-2, -cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-1,  cn_theta_*p.p_u*2.0*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i,   -cn_theta_*p.p_u*1.5*dx_) );
                }
            }
            else
            {
                // if peclet is negative
                if(i == (n_node_-2))
                {
                    // use central difference at right boundary-1
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i-1,  cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+1, -cn_theta_*p.p_u*0.5*dx_) );
                }
                else
                {
                    // use forward difference away from right boundary
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i  ,  cn_theta_*p.p_u*1.5*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+1, -cn_theta_*p.p_u*2.0*dx_) );
                    triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i+2,  cn_theta_*p.p_u*0.5*dx_) );
                }
            }

            // binding/unbinding
            triplet_list.push_back( T(c_u_offset_+i, c_u_offset_+i, cn_theta_*((p.alpha + p.gamma_u*phi(0,i))*dx_*dx_)) );

            // C_u / C_b
            triplet_list.push_back( T(c_u_offset_+i, c_b_offset_+i, cn_theta_*p.beta*dx_*dx_) );

            // C_u / C_s
            // (no entries)

            // C_u / phi
            triplet_list.push_back( T(c_u_offset_+i, phi_offset_+i, -cn_theta_*p.gamma_u*c_u(0,i)*dx_*dx_) );

            // C_b / C_u
            triplet_list.push_back( T(c_b_offset_+i, c_u_offset_+i, cn_theta_*p.alpha*dx_*dx_) );

            // C_b / C_b
            triplet_list.push_back( T(c_b_offset_+i, c_b_offset_+i, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,i))*dx_*dx_) );

            // C_b / C_s
            // (no entries)

            // C_b / phi
            triplet_list.push_back( T(c_b_offset_+i, phi_offset_+i, -cn_theta_*(p.gamma_b*c_b(0,i))*dx_*dx_) );

            // C_s / C_u
            triplet_list.push_back( T(c_s_offset_+i, c_u_offset_+i, cn_theta_*(p.gamma_u*phi(0,i))*dx_*dx_) );

            // C_s / C_b
            triplet_list.push_back( T(c_s_offset_+i, c_b_offset_+i, cn_theta_*(p.gamma_b*phi(0,i)*dx_*dx_)) );

            // C_s / C_s

            // time derivative
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i, -time_factor_*dx_*dx_/dt_) );

            // diffusion
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-1,  cn_theta_*p.D_su) );
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i,   -cn_theta_*2.0*p.D_su) );
            triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+1,  cn_theta_*p.D_su) );

            // advection (upwinding)
            if(p.p_u > 0)
            {
                // if peclet is positive
                if(i == 1)
                {
                    // use central difference at left boundary+1
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-1,  cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+1, -cn_theta_*p.p_u*0.5*dx_) );
                }
                else
                {
                    // use backward difference away from left boundary
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-2, -cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-1,  cn_theta_*p.p_u*2.0*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i,   -cn_theta_*p.p_u*1.5*dx_) );
                }
            }
            else
            {
                // if peclet is negative
                if(i == (n_node_-2))
                {
                    // use central difference at right boundary-1
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i-1,  cn_theta_*p.p_u*0.5*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+1, -cn_theta_*p.p_u*0.5*dx_) );
                }
                else
                {
                    // use forward difference away from right boundary
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i,    cn_theta_*p.p_u*1.5*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+1, -cn_theta_*p.p_u*2.0*dx_) );
                    triplet_list.push_back( T(c_s_offset_+i, c_s_offset_+i+2,  cn_theta_*p.p_u*0.5*dx_) );
                }
            }

            // C_s / phi
            triplet_list.push_back( T(c_s_offset_+i, phi_offset_+i, cn_theta_*(p.gamma_u*c_u(0,i) + p.gamma_b*c_b(0,i))*dx_*dx_) );

            // phi / C_u
            // (no entries)

            // phi / C_b

            const double d_c_b = stencil_1_central(0,Variable::c_b,i);
            const double nu_d_c_b = p.nu*d_c_b;
            const double d_phi = upwind_d_phi(0,i,nu_d_c_b);

            // first chemotaxis term
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i-1,  cn_theta_*0.5*p.nu*d_phi) );
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i+1, -cn_theta_*0.5*p.nu*d_phi) );

            // second chemotaxis term
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i-1, -cn_theta_*1.0*p.nu*phi(0,i)) );
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i,    cn_theta_*2.0*p.nu*phi(0,i)) );
            triplet_list.push_back( T(phi_offset_+i, c_b_offset_+i+1, -cn_theta_*1.0*p.nu*phi(0,i)) );

            // phi / C_s
            // (no entries)

            // phi / phi

            // time derivative
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i, -time_factor_*dx_*dx_/dt_) );

            // diffusion
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-1,  cn_theta_*p.D_ju) );
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i,   -cn_theta_*2.0*p.D_ju) );
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+1,  cn_theta_*p.D_ju) );

            // chemotaxis (upwinding)
            if(nu_d_c_b > 0)
            {
                // if nu_d_c_b is positive
                if(i == 1)
                {
                    // use central difference at left boundary+1
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-1,  cn_theta_*nu_d_c_b*0.5) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+1, -cn_theta_*nu_d_c_b*0.5) );
                }
                else
                {
                    // use backward difference away from left boundary
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-2, -cn_theta_*nu_d_c_b*0.5) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-1,  cn_theta_*nu_d_c_b*2.0) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i,   -cn_theta_*nu_d_c_b*1.5) );
                }
            }
            else
            {
                // if nu_d_c_b is negative
                if(i == (n_node_-2))
                {
                    // use central difference at right boundary-1
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i-1,  cn_theta_*nu_d_c_b*0.5) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+1, -cn_theta_*nu_d_c_b*0.5) );
                }
                else
                {
                    // use forward difference away from right boundary
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i  ,  cn_theta_*nu_d_c_b*1.5) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+1, -cn_theta_*nu_d_c_b*2.0) );
                    triplet_list.push_back( T(phi_offset_+i, phi_offset_+i+2,  cn_theta_*nu_d_c_b*0.5) );
                }
            }

            // second term of chemotaxis
            triplet_list.push_back( T(phi_offset_+i, phi_offset_+i, -cn_theta_*p.nu*(c_b(0,i-1) - 2.0*c_b(0,i) + c_b(0,i+1))) );
        }

        // RHS boundary conditions/equations
        // --------------------------------------------------------------------

        // C_u / C_u
        triplet_list.push_back( T(c_u_offset_+n_node_-1, c_u_offset_+n_node_-1, 1.0*dx_*dx_) );

        // C_b / C_u
        triplet_list.push_back( T(c_b_offset_+n_node_-1, c_u_offset_+n_node_-1, cn_theta_*p.alpha*dx_*dx_) );

        // C_b / C_b
        triplet_list.push_back( T(c_b_offset_+n_node_-1, c_b_offset_+n_node_-1, -time_factor_*dx_*dx_/dt_ - cn_theta_*(p.beta + p.gamma_b*phi(0,n_node_-1))*dx_*dx_) );

        // C_b / C_s
        // (no entries)

        // C_b / phi
        triplet_list.push_back( T(c_b_offset_+n_node_-1, phi_offset_+n_node_-1, -cn_theta_*(p.gamma_b*c_b(0,n_node_-1))*dx_*dx_) );

        // C_s / C_s
        triplet_list.push_back( T(c_s_offset_+n_node_-1, c_s_offset_+n_node_-1, 1.0*dx_*dx_) );

        // phi / phi
        triplet_list.push_back( T(phi_offset_+n_node_-1, phi_offset_+n_node_-1, dx_*dx_) );

        // construct matrix
        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }
};

int main(int argc, char **argv)
{
    if(argc != 5 && argc != 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " config_file n_node dt t_max [ output_interval ]\n";
        std::exit(1);
    }

    const unsigned n_node = std::atoi(argv[2]);
    const double dt = std::atof(argv[3]);
    const double t_max = std::atof(argv[4]);

    unsigned output_interval = 1;

    if(argc == 6)
    {
        output_interval = std::atoi(argv[5]);
    }

    ChemokinesProblem1D problem(n_node, dt);
    ChemokinesProblem1D::Max_residual = 1e-10;

    std::ifstream config_file(argv[1]);
    ConfigFile cf(config_file);

    problem.p.p_u     = cf.get<double>("p_u");
    problem.p.alpha   = cf.get<double>("alpha");
    problem.p.beta    = cf.get<double>("beta");
    problem.p.gamma_u = cf.get<double>("gamma_u");
    problem.p.gamma_b = cf.get<double>("gamma_b");
    problem.p.D_su    = cf.get<double>("D_su");
    problem.p.D_ju    = cf.get<double>("D_ju");
    problem.p.nu      = cf.get<double>("nu");
    problem.p.lambda  = cf.get<double>("lambda");

    std::cout << "p_u     = " << problem.p.p_u     << '\n';
    std::cout << "alpha   = " << problem.p.alpha   << '\n';
    std::cout << "beta    = " << problem.p.beta    << '\n';
    std::cout << "gamma_u = " << problem.p.gamma_u << '\n';
    std::cout << "gamma_b = " << problem.p.gamma_b << '\n';
    std::cout << "D_su    = " << problem.p.D_su    << '\n';
    std::cout << "D_ju    = " << problem.p.D_ju    << '\n';
    std::cout << "nu      = " << problem.p.nu      << '\n';
    std::cout << "lambda  = " << problem.p.lambda  << '\n';

    problem.enable_terse_logging();

    char filename[200];
    std::ofstream outfile;

    bool dump = false;
    //bool dump = true;

    problem.Max_newton_iterations = 100;

    // perform a steady solve and output it
    problem.steady_solve(dump);
    std::sprintf(filename, "output_steady.csv");
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    // set initial conditions
    problem.set_initial_conditions();

    // output initial conditions
    std::sprintf(filename, "output_%05i.csv", 0);
    outfile.open(filename);
    problem.output(outfile);
    outfile.close();

    unsigned i = 1;

    // timestepping loop
    while(problem.time() <= t_max)
    {
        // solve for current timestep
        problem.unsteady_solve(dump);

        if(i % output_interval == 0)
        {
            // output current solution
            //std::cout << "Outputting solution at time = " << problem.time() << '\n';
            std::cout << ";\tOutputting";
            std::sprintf(filename, "output_%05i.csv", i/output_interval);
            outfile.open(filename);
            problem.output(outfile);
            outfile.close();
        }

        ++i;
    }
}
