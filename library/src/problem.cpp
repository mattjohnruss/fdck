#include <problem.h>

#include <iostream>

namespace mjrfd
{
    Problem::Problem(const unsigned n_dof, const double dt) :
        u_(n_dof),
        u_old_(n_dof),
        du_(n_dof),
        residual_(n_dof),
        jacobian_(n_dof, n_dof),
        time_(0.0),
        dt_(dt)
    {
        u_.setZero();
        du_.setZero();
    }

    Problem::~Problem()
    {
    }

    void Problem::solve()
    {
        // backup the current solution before solving
        u_old_ = u_;

        // the Newton iteration will continue until stop == true
        bool stop = false;

        // counter for the number of Newton iterations
        unsigned count = 0;

        // Newton method
        while(!stop)
        {
            // calculate the residuals and jacobian for the current solution
            calculate_residual();
            calculate_jacobian();

            // find the l^\infty norm of the residual vector
            double max_residual = residual_.lpNorm<Eigen::Infinity>();

            std::cout << "Maximum residual: " << max_residual << '\n';

            // test if the current max residual is greater than the threshold
            if(max_residual > Max_residual)
            {
                ++count; 

                std::cout << "Exceeds " << Max_residual
                          << ": performing Newton iteration\n";

                // solve the linear system jacobian_*du_ = -residual_ for du_
                linear_solve();

                // update the solution
                u_ += du_;
            }
            else
            {
                // the max residual is less than the threshold so stop
                stop = true;
                std::cout << "Newton method converged in " << count << " "
                          << ((count == 1) ? "iteration" : "iterations") << '\n';
            }
        }

    }

    void Problem::unsteady_solve()
    {
        time_ += dt_;
        Problem::solve();
    }

    const double Problem::time() const
    {
        return time_;
    }

    void Problem::linear_solve()
    {
        linear_solver_.compute(jacobian_);

        if(linear_solver_.info() != Eigen::Success)
        {
            std::cerr << "Compute failed!\n";
        }

        du_ = linear_solver_.solve(-residual_);
    }

    double Problem::Max_residual = 1.0e-8;
}
