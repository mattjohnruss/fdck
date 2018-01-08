#include <problem.h>

#include <iostream>
#include <iomanip>
#include <fstream>

namespace mjrfd
{
    Problem::Problem(const unsigned n_dof, const double dt) :
        u_(2, Eigen::VectorXd(n_dof)),
        //u_old_(n_dof),
        du_(n_dof),
        residual_(n_dof),
        jacobian_(n_dof, n_dof),
        time_(0.0),
        dt_(dt),
        n_dof_(n_dof)
    {
        for(unsigned i = 0; i < 2; ++i)
        {
            u_[i].setZero();
        }

        du_.setZero();
    }

    Problem::~Problem()
    {
    }

    void Problem::solve(bool dump)
    {
        // backup the current solution before solving
        //u_old_ = u_;
        u_[1] = u_[0];

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

            if(dump)
            {
                char filename[100];

                std::snprintf(filename, 100, "res_t=%f_%u.dat", time_, count);
                std::ofstream res_file(filename);

                std::snprintf(filename, 100, "jac_t=%f_%u.dat", time_, count);
                std::ofstream jac_file(filename);

                dump_res_and_jac(res_file, jac_file);

                res_file.close();
                jac_file.close();
            }

            // find the l^\infty norm of the residual vector
            //double max_residual = residual_.lpNorm<Eigen::Infinity>();
            double max_residual = 0.0;
            unsigned max_residual_index = 0;

            for(unsigned i = 0; i < n_dof_; ++i)
            {
                if(std::abs(residual_(i)) > max_residual)
                {
                    max_residual = std::abs(residual_(i));
                    max_residual_index = i;
                }
            }

            std::cout << "Maximum residual (" << max_residual_index << "): "
                      << max_residual << '\n';

            // test if the current max residual is greater than the threshold
            if(max_residual > Max_residual)
            {
                ++count; 

                std::cout << "Exceeds " << Max_residual
                          << ": performing Newton iteration\n";

                // solve the linear system jacobian_*du_ = -residual_ for du_
                linear_solve();

                // update the solution
                u_[0] += du_;
            }
            else
            {
                // the max residual is less than the threshold so stop
                stop = true;
                std::cout << "Newton method converged (with maximum residual "
                          << max_residual << ")"
                          << " in " << count << " "
                          << ((count == 1) ? "iteration" : "iterations") << '\n';
            }

            if(count >= Max_newton_iterations)
            {
                // the maximum number of newton steps has been exceeded so stop
                stop = true;
                std::cout << "The maximum number of Newton iterations ("
                          << Max_newton_iterations << ") "
                          << "has been exceeded\nExiting\n";

                std::exit(1);
            }
        }
    }

    void Problem::unsteady_solve(bool dump)
    {
        time_ += dt_;
        std::cout << "\nSolving at time = " << time_ << "\n\n";
        Problem::solve(dump);
    }

    void Problem::dump_res_and_jac(std::ostream &res_stream, std::ostream &jac_stream) const
    {
        // sensible output format for Eigen matrices
        const Eigen::IOFormat plain_fmt(16, Eigen::DontAlignCols);

        // output residuals
        res_stream << residual_.format(plain_fmt);

        // output jacobian
        jac_stream << std::setprecision(16);

        for(int i = 0; i < jacobian_.outerSize(); ++i)
        {
            for(Eigen::SparseMatrix<double>::InnerIterator it(jacobian_,i); it; ++it)
            {
                jac_stream << it.row() << " " << it.col() << " " << it.value() << "\n";
            }
        }

        //jac_stream << Eigen::MatrixXd(jacobian_).format(plain_fmt);
    }

    const double Problem::time() const
    {
        return time_;
    }

    const double Problem::u(unsigned i) const
    {
        return u_[0](i);
    }

    const double Problem::u(unsigned t, unsigned i) const
    {
        assert(t < 2);

        return u_[t](i);
    }

    double& Problem::u(unsigned i)
    {
        return u_[0](i);
    }

    double& Problem::u(unsigned t, unsigned i)
    {
        assert(t < 2);

        return u_[t](i);
    }

    void Problem::linear_solve()
    {
        linear_solver_.compute(jacobian_);

        if(linear_solver_.info() != Eigen::Success)
        {
            std::cerr << "Problem::linear_solve() - Compute failed!\n";
        }

        du_ = linear_solver_.solve(-residual_);
    }

    double Problem::Max_residual = 1.0e-8;
    unsigned Problem::Max_newton_iterations = 20;
}
