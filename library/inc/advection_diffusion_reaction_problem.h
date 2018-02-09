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
        virtual void get_v(const unsigned node,
                           const std::vector<double> &u,
                           std::vector<double> &v) const = 0;

        /// Get the derivatives of then advection (or chemotaxis etc) velocities wrt dofs
        virtual void get_dv_du(const unsigned node,
                               const std::vector<double> &u,
                               std::vector<std::vector<double>> &dv_du) const = 0;
        
        /// Get the vector of reactions
        virtual void get_r(const std::vector<double> &u,
                           std::vector<double> &r) const = 0;
        
        /// Get the derivatives of reactions wrt dofs
        virtual void get_dr_du(const std::vector<double> &u,
                               std::vector<std::vector<double>> &dr_du) const = 0;

        double diff_res_helper(const unsigned var,
                               const unsigned i);

        void diff_jac_helper(std::vector<T> &triplet_list,
                             const unsigned var,
                             const unsigned i,
                             const double d);

        /// The number of nodes
        const unsigned n_node_;

        /// Crank-Nicolson theta
        double cn_theta_;

        /// Spatial step size
        const double dx_;
    };

    /// Contructor
    AdvectionDiffusionReactionProblem::AdvectionDiffusionReactionProblem(
        const unsigned n_var,
        const unsigned n_node,
        const double dt) :
        Problem(n_var, n_node),
        n_node_(n_node),
        cn_theta_(1.0),
        dx_(1.0/(n_node-1))
    {
    }

    /// Destructor
    AdvectionDiffusionReactionProblem::~AdvectionDiffusionReactionProblem()
    {
    }

    /// Output the solution
    void AdvectionDiffusionReactionProblem::output(std::ostream &out) const
    {
        for(unsigned i = 0; i < n_node_; ++i)
        {
            out << static_cast<double>(i)/static_cast<double>(n_node_-1) << ' ';

            for(unsigned var = 0; var < n_var_; ++var)
            {
                out << u(0, var, i) << ' ';
            }

            out << '\n';
        }
    }

    /// Calculate the residual vector
    void AdvectionDiffusionReactionProblem::calculate_residual()
    {
        residual_.setZero();
        // Storage for the diffusion coefficients
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
            //for(unsigned i = 0; i < n_node_; ++i)
            for(unsigned i = 0; i < n_node_; ++i)
            {
                // Calculate the index of the current dof
                unsigned index = var*n_node_ + i;

                if(i == 0)
                {
                    // FIXME this is a temporary hack
                    residual_(index) = (1.0 - u(0, var, i))*dx_*dx_;
                }
                else if(i == n_node_-1)
                {
                    // FIXME this is a temporary hack
                    residual_(index) = (0.0 - u(0, var, i))*dx_*dx_;
                }
                else
                {
                    // Time derivatives
                    // TODO

                    // Diffusion
                    residual_(index) -= d[var]*diff_res_helper(var, i);

                    // Advection (or chemotaxis etc)
                    // TODO

                    // Reaction
                    // TODO
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
            //for(unsigned i = 0; i < n_node_; ++i)
            for(unsigned i = 0; i < n_node_; ++i)
            {
                // Calculate the index of the current dof
                unsigned index = var*n_node_ + i;

                // Diffusion
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
                    diff_jac_helper(triplet_list, var, i, d[var]);
                }
            }
        }

        jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        jacobian_.makeCompressed();

        //std::cout << Eigen::MatrixXd(jacobian_).determinant() << '\n';
        std::cout << std::setprecision(12) << Eigen::MatrixXd(jacobian_) << '\n';
    }

    double AdvectionDiffusionReactionProblem::diff_res_helper(const unsigned var,
                                                              const unsigned i)
    {
        double result = 0.0;

        // Loop over the stencil points and get the offset j (relative to i)
        // and the weight w
        for(const auto& [j, w] : stencil::central_2::weights)
        {
            // Sum the contributions from the stencil points
            result += w*(cn_theta_*u(0, var, i+j) + (1.0-cn_theta_)*u(1, var, i+j));
        }

        return result;
    }

    void AdvectionDiffusionReactionProblem::diff_jac_helper(std::vector<T> &triplet_list,
                                                            const unsigned var,
                                                            const unsigned i,
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
}
