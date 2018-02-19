#pragma once

#include <problem.h>
#include <stencil.h>

#include <iostream>

namespace mjrfd
{
    //class AdvectionDiffusionReactionProblem : public FiniteDifferenceProblem
    class AdvectionDiffusionReactionProblem : public Problem
    {
        typedef Eigen::Triplet<double> T;

    public:
        enum class Boundary
        {
            Left,
            Right
        };

        /// Constructor
        AdvectionDiffusionReactionProblem(const unsigned n_var,
                                          const unsigned n_node,
                                          const double dt = 0.0);

        /// Destructor
        virtual ~AdvectionDiffusionReactionProblem();

        void make_steady() override;

        void make_unsteady() override;

        /// Output the solution
        void output(std::ostream &out) const override;

        /// Calculate the residual vector
        void calculate_residual() override;

        /// Calculate the jacobian matrix
        void calculate_jacobian() override;

        void enable_bc(Boundary b, const std::vector<unsigned> &vars);
        void disable_bc(Boundary b, const std::vector<unsigned> &vars);

    private:
        /// Get the coefficients a1, a2, a3 for the (Robin) boundary conditions
        virtual void get_bc(Boundary b,
                            std::vector<double> &a1,
                            std::vector<double> &a2,
                            std::vector<double> &a3) const = 0;

        /// Get the vector of diffusivities
        virtual void get_d(std::vector<double> &d) const = 0;

        /// Get the vector of advection (or chemotaxis etc) velocities
        virtual void get_v(const unsigned i,
                           const unsigned var,
                           double &v) const = 0;

        /// Get the derivatives of then advection (or chemotaxis etc) velocities wrt dofs
        virtual void get_dv_du(const unsigned i,
                               const std::vector<double> &u,
                               std::vector<std::vector<double>> &dv_du) const = 0;
        
        /// Get the vector of reactions
        virtual void get_r(const std::vector<double> &u,
                           std::vector<double> &r) const = 0;
        
        /// Get the derivatives of reactions wrt dofs
        virtual void get_dr_du(const std::vector<double> &u,
                               std::vector<std::vector<double>> &dr_du) const = 0;

        const std::unordered_map<int, double>&
        upwind_stencil_weights(const unsigned i, const double v) const;

    protected:
        /// Calculate the x coordinate from the node number i
        const double x(const unsigned i) const;

        /// The number of nodes
        const unsigned n_node_;

    private:
        /// Crank-Nicolson theta
        double cn_theta_;

        /// Backup for Crank-Nicolson theta
        double cn_theta_backup_;

        /// Spatial step size
        const double dx_;

        /// Time factor for switching between steady/unsteady solutions
        double time_factor_;

        std::vector<bool> left_bc_;
        std::vector<bool> right_bc_;
    };

    /// Contructor
    AdvectionDiffusionReactionProblem::AdvectionDiffusionReactionProblem(
        const unsigned n_var,
        const unsigned n_node,
        const double dt) :
        Problem(n_var, n_node, dt),
        n_node_(n_node),
        cn_theta_(1.0),
        dx_(1.0/(n_node-1)),
        time_factor_(1.0),
        left_bc_(n_var, false),
        right_bc_(n_var, false)
    {
    }

    /// Destructor
    AdvectionDiffusionReactionProblem::~AdvectionDiffusionReactionProblem()
    {
    }

    void AdvectionDiffusionReactionProblem::make_steady()
    {
        Problem::make_steady();
        time_factor_ = 0.0;
        cn_theta_backup_ = cn_theta_;
        cn_theta_ = 1.0;
    }

    void AdvectionDiffusionReactionProblem::make_unsteady()
    {
        Problem::make_unsteady();
        time_factor_ = 1.0;
        cn_theta_ = cn_theta_backup_;
    }

    /// Output the solution
    void AdvectionDiffusionReactionProblem::output(std::ostream &out) const
    {
        out << "x ";

        for(unsigned var = 0; var < n_var_; ++var)
        {
            out << "c" << var;
            out << (var < (n_var_-1) ? ' ' : '\n');
        }

        for(unsigned i = 0; i < n_node_; ++i)
        {
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' ';

            for(unsigned var = 0; var < n_var_; ++var)
            {
                out << u(0, var, i);
                out << (var < (n_var_-1) ? ' ' : '\n');
            }
        }
    }

    /// Calculate the residual vector
    void AdvectionDiffusionReactionProblem::calculate_residual()
    {
        // Set the residuals to zero
        residual_.setZero();

        // Storage for the diffusion coefficients
        std::vector<double> d(n_var_);

        // Get the diffusion coefficients. Do it here because they are
        // constants for each variable
        get_d(d);

        // Get the left boundary condition coefficients
        std::vector<double> a1_left(n_var_);
        std::vector<double> a2_left(n_var_);
        std::vector<double> a3_left(n_var_);
        get_bc(Boundary::Left, a1_left, a2_left, a3_left);

        // Get the right boundary condition coefficients
        std::vector<double> a1_right(n_var_);
        std::vector<double> a2_right(n_var_);
        std::vector<double> a3_right(n_var_);
        get_bc(Boundary::Right, a1_right, a2_right, a3_right);

        // Storage for the advection velocity at a node
        double v = 0;

        std::vector<double> r(n_var_);
        std::vector<double> r_old(n_var_);

        // Storage for u (all vars) at a node
        std::vector<double> u_at_node(n_var_);
        std::vector<double> u_old_at_node(n_var_);

        // Loop over the nodes
        for(unsigned i = 0; i < n_node_; ++i)
        {
            for(unsigned var = 0; var < n_var_; ++var)
            {
                u_at_node[var]     = u(0, var, i);
                u_old_at_node[var] = u(1, var, i);
            }

            // Get the reaction term at the current node
            get_r(u_at_node, r);

            // Loop over the variables
            for(unsigned var = 0; var < n_var_; ++var)
            {
                // Calculate the index of the current dof
                const unsigned index = var*n_node_ + i;

                // Boundary conditions
                if(i == 0 && left_bc_[var] == true)
                {
                    residual_(index) += (a1_left[var]*u(0, var, i) - a3_left[var])*dx_*dx_;

                    for(const auto& [j, w] : stencil::forward_1::weights)
                    {
                        residual_(index) += w*a2_left[var]*u(0, var, i+j)*dx_;
                    }
                }
                else if(i == n_node_-1 && right_bc_[var] == true)
                {
                    residual_(index) += (a1_right[var]*u(0, var, i) - a3_right[var])*dx_*dx_;

                    for(const auto& [j, w] : stencil::backward_1::weights)
                    {
                        residual_(index) += w*a2_right[var]*u(0, var, i+j)*dx_;
                    }
                }
                else
                {
                    // Time derivatives
                    residual_(index) += time_factor_*dx_*dx_/dt_*(u(0, var, i) - u(1, var, i));

                    // Diffusion

                    // Loop over the stencil points and get the offset j (relative to i)
                    // and the weight w
                    for(const auto& [j, w] : stencil::central_2::weights)
                    {
                        // Sum the contributions from the stencil points
                        if(d[var] > 0.0)
                        {
                            residual_(index) += -d[var]*w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
                        }
                    }

                    // Advection (or chemotaxis etc)

                    // Get the velocity at the location of the target node
                    double v_at_node = 0.0;
                    get_v(i, var, v_at_node);

                    // Loop over the stencil points, using the upwind stencil helper,
                    // and get the offset j (relative to i) and the weight w
                    for(const auto& [j, w] : upwind_stencil_weights(i, v_at_node))
                    {
                        // Get the velocity at the current stencil point
                        get_v(i+j, var, v);

                        // Sum the contributions from the stencil points
                        if(v > 0.0)
                        {
                            residual_(index) += dx_*v*w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
                        }
                    }

                    // Reaction
                    residual_(index) += -dx_*dx_*(cn_theta_*r[var] + (1.0-cn_theta_)*r_old[var]);
                }
            }
        }
    }

    /// Calculate the jacobian matrix
    void AdvectionDiffusionReactionProblem::calculate_jacobian()
    {
        jacobian_.setZero();

        // Vector of triplets for constructing the sparse jacobian matrix
        // TODO maybe make triplet_list a member variable so that it doesn't
        // get reallocated every iteration
        std::vector<T> triplet_list;

        // TODO reserve the correct number of entries
        // triplet_list.reserve(*);

        std::vector<double> d(n_var_);

        //std::vector<std::vector<double>> dv_du(n_var_,);
        //std::vector<std::vector<double>> dr_du(n_var_,);

        // Get the diffusion coefficients
        get_d(d);

        // Get the left boundary condition coefficients
        std::vector<double> a1_left(n_var_);
        std::vector<double> a2_left(n_var_);
        std::vector<double> a3_left(n_var_);
        get_bc(Boundary::Left, a1_left, a2_left, a3_left);

        // Get the right boundary condition coefficients
        std::vector<double> a1_right(n_var_);
        std::vector<double> a2_right(n_var_);
        std::vector<double> a3_right(n_var_);
        get_bc(Boundary::Right, a1_right, a2_right, a3_right);

        // Storage for the advection velocity at a node
        double v = 0;

        std::vector<double> r(n_var_);
        std::vector<std::vector<double>> dr_du(n_var_, std::vector<double>(n_var_));

        // Storage for u (all vars) at a node
        std::vector<double> u_at_node(n_var_);

        // Loop over the nodes
        for(unsigned i = 0; i < n_node_; ++i)
        {
            for(unsigned var = 0; var < n_var_; ++var)
            {
                u_at_node[var] = u(0, var, i);
            }

            // Get the reaction term and derivatives at the current node
            get_r(u_at_node, r);
            get_dr_du(u_at_node, dr_du);

            // Loop over the variables
            for(unsigned var = 0; var < n_var_; ++var)
            {
                // Calculate the index of the current dof
                const unsigned index = var*n_node_ + i;

                if(i == 0 && left_bc_[var] == true)
                {
                    triplet_list.push_back( T(index, index, a1_left[var]*dx_*dx_) );

                    for(const auto& [j, w] : stencil::forward_1::weights)
                    {
                        triplet_list.push_back( T(index, index+j, w*a2_left[var]*dx_) );
                    }
                }
                else if(i == n_node_-1 && right_bc_[var] == true)
                {
                    triplet_list.push_back( T(index, index, a1_right[var]*dx_*dx_) );

                    for(const auto& [j, w] : stencil::backward_1::weights)
                    {
                        triplet_list.push_back( T(index, index+j, w*a2_right[var]*dx_) );
                    }
                }
                else
                {
                    // Time derivatives
                    triplet_list.push_back( T(index, index, time_factor_*dx_*dx_/dt_) );

                    // Diffusion
                    for(const auto& [j, w] : stencil::central_2::weights)
                    {
                        if(d[var] > 0.0)
                        {
                            triplet_list.push_back( T(index, index+j, -d[var]*w*cn_theta_) );
                        }
                    }

                    // Advection (or chemotaxis etc)
                    // TODO

                    // Reaction

                    // Loop over the variables again
                    for(unsigned var2 = 0; var2 < n_var_; ++var2)
                    {
                        unsigned index2 = var2*n_node_ + i;
                        triplet_list.push_back( T(index, index2, -dx_*dx_*cn_theta_*dr_du[var][var2]) );
                    }
                }
            }
        }

        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }

    void AdvectionDiffusionReactionProblem::enable_bc(Boundary b, const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            if(b == Boundary::Left)
            {
                left_bc_[var] = true;
            }
            else if(b == Boundary::Right)
            {
                right_bc_[var] = true;
            }
        }
    }

    void AdvectionDiffusionReactionProblem::disable_bc(Boundary b, const std::vector<unsigned> &vars)
    {
        for(auto &var : vars)
        {
            if(b == Boundary::Left)
            {
                left_bc_[var] = false;
            }
            else if(b == Boundary::Right)
            {
                right_bc_[var] = false;
            }
        }
    }

    const std::unordered_map<int, double>&
    AdvectionDiffusionReactionProblem::upwind_stencil_weights(const unsigned i,
                                                              const double v) const
    {
        assert(i >= 0 && i <= n_node_-1);

        if(v > 0.0)
        {
            if(i >= 2)
                return stencil::backward_1::weights;
            else if(i == 1)
                return stencil::central_1::weights;
            else if(i == 0)
                return stencil::forward_1::weights;
        }
        else
        {
            if(i <= n_node_-3)
                return stencil::forward_1::weights;
            else if(i == n_node_-2)
                return stencil::central_1::weights;
            else if(i == n_node_-1)
                return stencil::backward_1::weights;
        }
    }

    const double AdvectionDiffusionReactionProblem::x(const unsigned i) const
    {
        return static_cast<double>(i)/static_cast<double>(n_node_-1);
    }
}
