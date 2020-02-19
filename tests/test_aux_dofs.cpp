#include "advection_diffusion_reaction_problem.h"
#include "utilities.h"
#include "log.h"

#include <catch2/catch.hpp>

using namespace fdck;

enum AuxVariable
{
    y = 0,
};

class AuxDofProblem : public Problem
{
public:
    AuxDofProblem() :
        Problem(0, 0, 1),
        trace_count_(0)
    {
        Max_residual = 1.0e-14;
    }

    ~AuxDofProblem()
    {
    }

    void set_initial_conditions()
    {
        // Set everything to zero
        clear_solution();

        u_aux(y) = 1.0;
    }

    void trace(Eigen::VectorXd &solution, Eigen::VectorXd &exact_solution)
    {
        solution(trace_count_) = u_aux(y);
        exact_solution(trace_count_) = this->exact_solution(time());

        trace_count_++;
    }

    // dummy implementation of output
    void output(std::ostream &) const override
    {
    }

    // dummy implementation of output_exact
    void output_exact(std::ostream &) const override
    {
    }

    double exact_solution(double time) const
    {
        return std::exp(time);
    }

private:
    unsigned trace_count_;

    unsigned aux_dof_index(unsigned i) const
    {
        return n_var_*n_dof_per_var_ + i;
    }

    void calculate_residual(Eigen::VectorXd &residual) const override
    {
        // Add the ode
        residual(aux_dof_index(y)) += (u_aux(0, y) - u_aux(1, y)) - dt_*u_aux(y);
    }

    void calculate_jacobian(std::vector<Triplet> &triplet_list) const override
    {
        triplet_list.emplace_back(aux_dof_index(y), aux_dof_index(y), 1.0 - dt_);
    }
};

namespace fdck
{
    TEST_CASE( "Solve u' = u with backward Euler using auxiliary dofs", "[aux_dofs]" )
    {
        unsigned n_timestep = 1000;
        const double dt = 1.0/static_cast<double>(n_timestep);
        const double t_max = 1.0;

        Eigen::VectorXd solution(n_timestep + 1);
        Eigen::VectorXd exact_solution(n_timestep + 1);

        AuxDofProblem problem;

        // reduce log level to remove some noise during tests
        log::set_level("warn");

        problem.enable_terse_logging();

        // set initial conditions
        problem.set_initial_conditions();

        // output initial conditions
        problem.trace(solution, exact_solution);

        unsigned i = 1;

        // timestepping loop
        while(problem.time() < t_max)
        {
            // solve for current timestep
            problem.unsteady_solve(dt);
            problem.trace(solution, exact_solution);

            ++i;
        }

        FDCK_INFO("Reached t > t_max ({}) after performing {} timesteps", t_max, i-1);

        REQUIRE( utilities::l2_norm(solution - exact_solution, dt) <= 1.0e-3 );
    }
}
