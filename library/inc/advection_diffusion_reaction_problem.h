#pragma once

#include <problem.h>
#include <stencil.h>

#include <iostream>

namespace mjrfd
{
    class AdvectionDiffusionReactionProblem : public Problem
    {
        typedef Eigen::Triplet<double> T;

    public:
        /// Enum for the two boundaries in the 1D problem
        enum class Boundary
        {
            Left,
            Right
        };

        /// Constructor
        AdvectionDiffusionReactionProblem(const unsigned n_var,
                                          const unsigned n_node);

        /// Destructor
        virtual ~AdvectionDiffusionReactionProblem();

        /// Make the problem steady by disabling time derivatives
        void make_steady() override;

        /// Make the problem unsteady by enabling time derivatives
        void make_unsteady() override;

        /// Output the solution
        void output(std::ostream &out) const override;

        /// Output the exact solution
        void output_exact(std::ostream &out) const override;

        /// Calculate the residual vector
        void calculate_residual() override;

        /// Calculate the jacobian matrix
        void calculate_jacobian() override;

        /// Enable boundary condition on boundary b for all variables in vars
        void enable_bc(Boundary b, const std::vector<unsigned> &vars);

        /// Disable boundary condition on boundary b for all variables in vars
        void disable_bc(Boundary b, const std::vector<unsigned> &vars);

        /// Enable spatial terms for all variables in vars
        void enable_spatial_terms(const std::vector<unsigned> &vars);

        /// Disable spatial terms for all variables in vars
        void disable_spatial_terms(const std::vector<unsigned> &vars);

        /// Enable outputting the time in the first column of each file
        void enable_output_time_column();

        /// Enable outputting the time in the first column of each file
        void disable_output_time_column();

        /// Set the variable names
        void set_variable_names(const std::vector<std::string> &var_names);

        /// Calculate the spatial integral of a variable using trapizium rule
        /// (for now)
        double integrate_solution(const unsigned var) const;

    private:
        /// Get the coefficients a1, a2, a3 for the (Robin) boundary conditions
        virtual void get_bc(Boundary b,
                            std::vector<double> &a1,
                            std::vector<double> &a2,
                            std::vector<double> &a3) const = 0;

        /// Get the derivatives of the boundary condition coefficients wrt dofs
        /// at node i2
        virtual void get_dbc_du(Boundary b,
                                const unsigned i2,
                                std::vector<std::vector<double>> &da1_du,
                                std::vector<std::vector<double>> &da2_du,
                                std::vector<std::vector<double>> &da3_du) const;

        /// Get the vector of diffusivities
        virtual void get_d(const unsigned t,
                           const unsigned i,
                           std::vector<double> &d) const = 0;

        /// Get the derivatives of the diffusivities wrt dofs
        /// We assume that d at node i only depends on other quantities at node
        /// i, which means it cannot depend on derivatives of the variables etc
        /// This means we only need to pass node i to the derivatives, unlike
        /// dv_du which needs i and j
        /// Defaults to zero
        virtual void get_dd_du(const unsigned t,
                               const unsigned i,
                               std::vector<std::vector<double>> &dd_du) const;

        /// Get the vector of advection (or chemotaxis etc) velocities
        virtual void get_v(const unsigned i,
                           const unsigned var,
                           double &v) const = 0;

        /// Get the derivatives of the advection (or chemotaxis etc) velocities wrt dofs
        /// Defaults to zero
        virtual void get_dv_du(const unsigned i,
                               const unsigned var,
                               const unsigned i2,
                               const unsigned var2,
                               double &dv_du) const;

        /// Get the vector of reactions
        virtual void get_r(const std::vector<double> &u,
                           std::vector<double> &r) const = 0;
        
        /// Get the derivatives of reactions wrt dofs
        /// Defaults to zero
        virtual void get_dr_du(const std::vector<double> &u,
                               std::vector<std::vector<double>> &dr_du) const;

        virtual void exact_solution(const double time,
                                    const double x,
                                    std::vector<double> &sol) const;

    protected:
        /// Returns an upwinded stencil in the appropriate direction depending
        /// on the sign of v; if node is near or on a boundary, return a central
        /// or opposite direction stencil as appropriate
        const std::unordered_map<int, double>&
        upwind_stencil_weights(const unsigned i, const double v) const;

        /// Returns a central difference stencil if node i is in the bulk, and
        /// forward/backward stencils if i is on a boundary
        const std::unordered_map<int, double>&
        central_1_stencil_weights(unsigned i) const;

        /// Calculate the x coordinate from the node number i
        const double x(const unsigned i) const;

        /// The number of nodes
        const unsigned n_node_;

        /// Spatial step size
        const double dx_;

    private:
        /// Crank-Nicolson theta
        double cn_theta_;

        /// Backup for Crank-Nicolson theta
        double cn_theta_backup_;

        /// Time factor for switching between steady/unsteady solutions
        double time_factor_;

        /// Vector of boolean flags for whether the left boundary condition is
        /// enabled for each variable
        std::vector<bool> left_bc_;

        /// Vector of boolean flags for whether the right boundary condition is
        /// enabled for each variable
        std::vector<bool> right_bc_;

        /// Vector of boolean flags for whether spatial terms are enabled for
        /// each variable
        std::vector<bool> spatial_terms_;

        /// Bool for whether to put a column in each output file with the
        /// current time
        bool output_time_column_;

        /// Vector of "human-readable" variable names for each variable which
        /// are used as column headers in the output files
        std::vector<std::string> var_names_;
    };
}
