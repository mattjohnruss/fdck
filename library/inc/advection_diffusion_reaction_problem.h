#pragma once

//#include <finite_difference_problem.h>
#include <problem.h>
#include <stencil.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>

namespace mjrfd
{
    //class AdvectionDiffusionReactionProblem : public FiniteDifferenceProblem
    class AdvectionDiffusionReactionProblem : public Problem
    {
        typedef Eigen::Triplet<double> T;

    public:
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

    private:
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

        double diffusion_res_helper(const unsigned i,
                                    const unsigned var);

        double advection_res_helper(const unsigned i,
                                    const unsigned var);

        void diffusion_jac_helper(std::vector<T> &triplet_list,
                                  const unsigned i,
                                  const unsigned var,
                                  const double d);

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
        time_factor_(1.0)
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
        residual_.setZero();
        // Storage for the diffusion coefficients
        std::vector<double> d(n_var_);

        // Storage for the advection velocity at a node
        double v = 0;

        std::vector<double> r(n_var_);


        // Storage for u (all vars) at a node
        std::vector<double> u_at_node(n_var_);

        //std::vector<std::vector<double>> dv_du(n_var_,);
        //std::vector<std::vector<double>> dr_du(n_var_,);

        // Get the diffusion coefficients. Do it here because they are
        // constants for each variable
        get_d(d);

        // Loop over the nodes
        for(unsigned i = 0; i < n_node_; ++i)
        {
            for(unsigned var = 0; var < n_var_; ++var)
            {
                u_at_node[var] = u(0, var, i);
            }

            // Get the reaction term at the current node
            get_r(u_at_node, r);

            // Loop over the variables
            for(unsigned var = 0; var < n_var_; ++var)
            {
                // Calculate the index of the current dof
                const unsigned index = var*n_node_ + i;

                if(i == 0)
                {
                    // FIXME this is a temporary hack
                    residual_(index) += (1.0 - u(0, var, i))*dx_*dx_;
                }
                else if(i == n_node_-1)
                {
                    // FIXME this is a temporary hack
                    residual_(index) += (0.0 - u(0, var, i))*dx_*dx_;
                }
                else
                {
                    // Time derivatives
                    residual_(index) += time_factor_*dx_*dx_/dt_*(u(0, var, i) - u(1, var, i));

                    // Diffusion
                    //residual_(index) += -d[var]*diffusion_res_helper(i, var);

                    // Loop over the stencil points and get the offset j (relative to i)
                    // and the weight w
                    for(const auto& [j, w] : stencil::central_2::weights)
                    {
                        // Sum the contributions from the stencil points
                        residual_(index) += -d[var]*w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
                    }

                    // Advection (or chemotaxis etc)
                    //residual_(index) += dx_*advection_res_helper(i, var);
                    for(const auto& [j, w] : stencil::central_1::weights)
                    {
                        // Get the velocity at the current stencil point
                        get_v(i+j, var, v);

                        // Sum the contributions from the stencil points
                        residual_(index) += dx_*v*w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
                    }

                    // Reaction
                    residual_(index) += -dx_*dx_*r[var];
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
        //std::vector<double> v(n_var_);
        //std::vector<double> r(n_var_);
        //std::vector<std::vector<double>> dv_du(n_var_,);
        //std::vector<std::vector<double>> dr_du(n_var_,);

        // Get the diffusion coefficients
        get_d(d);

        // Loop over the variables
        for(unsigned var = 0; var < n_var_; ++var)
        {
            // Loop over the nodes
            for(unsigned i = 0; i < n_node_; ++i)
            {
                // Calculate the index of the current dof
                unsigned index = var*n_node_ + i;

                if(i == 0)
                {
                    // FIXME this is a temporary hack
                    triplet_list.push_back( T(index, index, -1.0*dx_*dx_) );
                }
                else if(i == n_node_-1)
                {
                    // FIXME this is a temporary hack
                    triplet_list.push_back( T(index, index, -1.0*dx_*dx_) );
                }
                else
                {
                    // Diffusion
                    diffusion_jac_helper(triplet_list, i, var, d[var]);
                }
            }
        }

        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();
    }

    double AdvectionDiffusionReactionProblem::diffusion_res_helper(const unsigned i,
                                                                   const unsigned var)
    {
        double result = 0.0;

        //// Loop over the stencil points and get the offset j (relative to i)
        //// and the weight w
        //for(const auto& [j, w] : stencil::central_2::weights)
        //{
            //// Sum the contributions from the stencil points
            //result += w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
        //}

        return result;
    }

    double AdvectionDiffusionReactionProblem::advection_res_helper(const unsigned i,
                                                                   const unsigned var)
    {
        double result = 0.0;

        // Loop over the stencil points and get the offset j (relative to i)
        // and the weight w
        // TODO use upwinding for the advection terms!
        //for(const auto& [j, w] : stencil::central_1::weights)
        //{
            //double v = 0;
            //get_v(i+j, var, v);

            //result += w*v*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
        //}

        return result;
    }

    void AdvectionDiffusionReactionProblem::diffusion_jac_helper(
        std::vector<T> &triplet_list,
        const unsigned i,
        const unsigned var,
        const double d)
    {
        unsigned index = var*n_node_ + i;

        // Loop over the stencil points
        for(const auto& [j, w] : stencil::central_2::weights)
        {
            unsigned index2 = var*n_node_ + i + j;
            triplet_list.push_back( T(index, index2, -w*cn_theta_) );
        }
    }

    const double AdvectionDiffusionReactionProblem::x(const unsigned i) const
    {
        return static_cast<double>(i)/static_cast<double>(n_node_-1);
    }
}
