#include <problem.h>

#include <iostream>
#include <iomanip>
#include <fstream>

namespace mjrfd
{
    Problem::Problem(const unsigned n_var,
                     const unsigned n_dof_per_var,
                     const double dt) :
        n_dof_per_var_(n_dof_per_var),
        n_dof_(n_var*n_dof_per_var),
        n_var_(n_var),
        n_time_history_(2),
        u_(n_time_history_, Eigen::VectorXd(n_dof_)),
        du_(Eigen::VectorXd(n_dof_)),
        residual_(n_dof_),
        jacobian_(n_dof_, n_dof_),
        time_(0.0),
        dt_(dt),
        steady_(false),
        terse_logging_(false),
        use_fd_jacobian_(false),
        dump_jacobian_(false)
    {
        clear_solution();
    }

    Problem::~Problem()
    {
    }

    void Problem::solve()
    {
        // backup the current solution before solving, unless the problem is steady
        if(steady_ == false)
        {
            // TODO generalise this for n_time_history_ > 1
            u_[1] = u_[0];
        }

        // the Newton iteration will continue until stop == true
        bool stop = false;

        // counter for the number of Newton iterations
        unsigned count = 0;

        // Newton method
        while(!stop)
        {
            // calculate the residuals for the current state
            calculate_residual();

            //if(dump_jacobian_ == true)
            //{
                //// if we're dumping the jacobian, calculate it here
                //// this is inefficient, but less so than dumping the jacobian
                //if(use_fd_jacobian_ == true)
                //{
                    //Problem::calculate_jacobian();
                //}
                //else
                //{
                    //calculate_jacobian();
                //}

                //char filename[100];

                //std::snprintf(filename, 100, "res_t=%f_%u.dat", time_, count);
                //std::ofstream res_file(filename);

                //std::snprintf(filename, 100, "jac_t=%f_%u.dat", time_, count);
                //std::ofstream jac_file(filename);

                //dump_res_and_jac(res_file, jac_file);

                //res_file.close();
                //jac_file.close();
            //}

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

            if(terse_logging_ == true && count == 0)
            {
                std::cout << "Max res = " << max_residual << ";\t[";
            }

            if(terse_logging_ == false)
            {
                std::cout << "Maximum residual (" << max_residual_index << "): "
                          << max_residual << '\n';
            }

            // test if the current max residual is greater than the threshold
            // TODO always do one Newton iteration regardless of the residuals
            if(max_residual > Max_residual)
            {
                ++count; 

                // calculate the jacobian for the current state
                if(use_fd_jacobian_ == true)
                {
                    Problem::calculate_jacobian();
                }
                else
                {
                    calculate_jacobian();
                }

                if(terse_logging_ == true)
                {
                    std::cout << "x";
                }
                else
                {
                    std::cout << "Exceeds " << Max_residual
                              << ": performing Newton iteration\n";
                }

                if(dump_jacobian_ == true)
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

                // solve the linear system jacobian_*du_ = -residual_ for du_
                linear_solve();

                // update the solution
                u_[0] += du_;
            }
            else
            {
                // the max residual is less than the threshold so stop
                stop = true;

                if(terse_logging_ == true)
                {
                    std::cout << "];\tMax res = " << max_residual;
                }
                else
                {
                    std::cout << "Newton method converged (with maximum residual "
                              << max_residual << ")"
                              << " in " << count << " "
                              << ((count == 1) ? "iteration" : "iterations") << '\n';
                }
            }

            if(count >= Max_newton_iterations)
            {
                // the maximum number of newton steps has been exceeded so stop
                stop = true;
                std::cout << "\nThe maximum number of Newton iterations ("
                          << Max_newton_iterations << ") "
                          << "has been exceeded\n";

                // If we're doing a steady solve, keep going in case we want to
                // timestep after it fails. Otherwise, exit
                if(steady_ == false)
                {
                    std::cout << "Exiting\n";
                    std::exit(1);
                }
            }
        }
    }

    void Problem::steady_solve()
    {
        // check if time derivatives enabled
        bool steady = is_steady();

        // disable time derivatives if they are enabled
        if(steady == false)
        {
            make_steady();
        }

        // solve the steady problem
        solve();

        // re-enable time derivatives if they were enabled before
        if(steady == false)
        {
            make_unsteady();
        }
    }

    void Problem::unsteady_solve()
    {
        time_ += dt_;

        if(terse_logging_ == true)
        {
            std::cout << "\nTime = " << time_ << ";\t";
        }
        else
        {
            std::cout << "\nSolving at time = " << time_ << "\n\n";
        }

        Problem::solve();
    }

    const bool Problem::is_steady() const
    {
        return steady_;
    }

    void Problem::make_steady()
    {
        steady_ = true;
    }

    void Problem::make_unsteady()
    {
        steady_ = false;
    }

    const double Problem::time() const
    {
        return time_;
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

        // if there isn't an entry in the bottom right-hand corner, output a zero
        if(jacobian_.rows() != n_dof_ || jacobian_.cols() != n_dof_)
        {
            std::cout << "Problem::dump_res_and_jac - outputting zero\n";
            jac_stream << n_dof_-1 << " " << n_dof_-1 << " 0.0\n";
        }

        //jac_stream << Eigen::MatrixXd(jacobian_).format(plain_fmt);
    }

    void Problem::enable_terse_logging()
    {
        terse_logging_ = true;
    }

    void Problem::disable_terse_logging()
    {
        terse_logging_ = false;
    }

    void Problem::enable_fd_jacobian()
    {
        use_fd_jacobian_ = true;
    }

    void Problem::disable_fd_jacobian()
    {
        use_fd_jacobian_ = false;
    }

    void Problem::enable_dump_jacobian()
    {
        dump_jacobian_ = true;
    }

    void Problem::disable_dump_jacobian()
    {
        dump_jacobian_ = false;
    }

    void Problem::clear_solution()
    {
        for(unsigned t = 0; t < n_time_history_; ++t)
        {
            u_[t].setZero();
        }
    }

    // Default implementation which calculates the jacobian by
    // finite-differencing the residuals
    void Problem::calculate_jacobian()
    {
        // Storage for the dense jacobian matrix
        Eigen::MatrixXd dense_jacobian(n_dof_, n_dof_);

        // Backup the residuals before doing anything
        Eigen::VectorXd residual_backup = residual_;

        // Storage for the residuals after incrementing/decrementing a dof
        Eigen::VectorXd residual_plus(residual_.size());
        Eigen::VectorXd residual_minus(residual_.size());

        // Derivative wrt dof i
        for(unsigned i = 0; i < n_dof_; ++i)
        {
            // Backup the i-th dof
            double u_i_backup = u_[0](i);

            // Calculate plus residual
            u_[0](i) += Jacobian_fd_step;
            calculate_residual();
            residual_plus = residual_;
            u_[0](i) = u_i_backup;

            // Calculate minus residual
            u_[0](i) -= Jacobian_fd_step;
            calculate_residual();
            residual_minus = residual_;
            u_[0](i) = u_i_backup;

            // Set the column of the jacobian matrix
            dense_jacobian.col(i) =
                (residual_plus - residual_minus)/(2.0*Jacobian_fd_step);
        }

        //std::cout << std::setprecision(12) << dense_jacobian << '\n';

        // Set the sparse jacobian from the dense one
        // TODO check the tolerance for throwing away terms close to zero
        jacobian_ = dense_jacobian.sparseView();
        jacobian_.makeCompressed();

        // Restore the residuals to their previous state
        residual_ = residual_backup;
    }

    void Problem::linear_solve()
    {
        linear_solver_.compute(jacobian_);

        if(linear_solver_.info() != Eigen::Success)
        {
            std::cerr << "Problem::linear_solve() - Compute failed!\n";
        }

        // get the update vector as a single array
        du_ = linear_solver_.solve(-residual_);
    }

    double Problem::Max_residual = 1.0e-8;
    unsigned Problem::Max_newton_iterations = 20;

    double Problem::Jacobian_fd_step = 1.0e-8;
}
