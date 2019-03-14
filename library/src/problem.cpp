#include "problem.h"
#include "log.h"

#include <iostream>
#include <iomanip>
#include <fstream>

namespace mjrfd
{
    Problem::Problem(const unsigned n_var,
                     const unsigned n_dof_per_var,
                     const unsigned n_aux_dof,
                     const unsigned n_previous_values) :
        n_dof_per_var_(n_dof_per_var),
        n_dof_(n_var*n_dof_per_var + n_aux_dof),
        n_var_(n_var),
        n_aux_dof_(n_aux_dof),
        n_previous_values_(n_previous_values),
        u_(n_previous_values_+1, Eigen::VectorXd(n_dof_)),
        du_(Eigen::VectorXd(n_dof_)),
        residual_(n_dof_),
        jacobian_(n_dof_, n_dof_),
        time_(0.0),
        dt_(0.01),
        steady_(false),
        Max_residual(1.0e-8),
        Max_newton_iterations(20),
        Jacobian_fd_step(1.0e-8),
        Fd_jacobian_reserved_entries(10*n_dof_),
        Fd_jacobian_threshold(1.0e-14),
        terse_logging_(true),
        use_fd_jacobian_(false),
        dump_jacobian_(false),
        exit_on_solve_fail_(false),
        jac_filename_prefix_("")
    {
        clear_solution();
    }

    Problem::~Problem()
    {
    }

    void Problem::solve()
    {
        // shift the time history before solving, unless the problem is steady
        if(steady_ == false)
        {
            for(unsigned t = 0; t < n_previous_values_; ++t)
            {
                u_[t+1] = u_[t];
            }
        }

        // the Newton iteration will continue until stop == true
        bool stop = false;

        // counter for the number of Newton iterations
        unsigned count = 0;

        actions_before_solve();

        // keep the max residual before solving as we output after solving when
        // terse logging in enabled
        double max_residual_before_solve = 0;

        // Newton method
        while(!stop)
        {
            // clear then calculate the residuals for the current state
            residual_.setZero();
            calculate_residual(residual_);

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

            if(count == 0)
                max_residual_before_solve = max_residual;

            if(terse_logging_ == false)
            {
                MJRFD_LIB_INFO("Maximum residual ({}): {}",
                               max_residual_index,
                               max_residual);
            }

            // test if the current max residual is greater than the threshold
            // always do one Newton iteration regardless of the residuals
            if(max_residual > Max_residual || count == 0)
            {
                std::vector<Triplet> triplet_list;

                // calculate the jacobian for the current state
                if(use_fd_jacobian_ == true)
                {
                    // Reserve enough memory for an estimated 10*n_dof entries
                    // by default. Should be enough for any discretisation that
                    // creates (block) banded jacobians but can be customised
                    // in derived classes.
                    jacobian_.reserve(Fd_jacobian_reserved_entries);
                    calculate_jacobian_fd(triplet_list);
                }
                else
                {
                    calculate_jacobian(triplet_list);
                }

                jacobian_.setFromTriplets(triplet_list.begin(), triplet_list.end());
                jacobian_.makeCompressed();

                if(terse_logging_ == false)
                {
                    MJRFD_LIB_INFO("Exceeds {}: performing Newton iteration",
                                   Max_residual);
                }

                if(dump_jacobian_ == true)
                {
                    char filename[100];

                    std::snprintf(filename, 100, "%sres_t=%f_%u.dat",
                                  jac_filename_prefix_.c_str(), time_, count);
                    std::ofstream res_file(filename);

                    std::snprintf(filename, 100, "%sjac_t=%f_%u.dat",
                                  jac_filename_prefix_.c_str(), time_, count);
                    std::ofstream jac_file(filename);

                    dump_res_and_jac(res_file, jac_file);

                    res_file.close();
                    jac_file.close();
                }

                // solve the linear system jacobian_*du_ = -residual_ for du_
                linear_solve();

                // update the solution
                u_[0] += du_;

                // Increment the counter
                ++count;
            }
            else
            {
                // the max residual is less than the threshold so stop
                stop = true;

                if(terse_logging_ == true)
                {
                    if(steady_)
                    {
                        MJRFD_LIB_INFO("Max res = {}\t[{:x<{}}]\tMax res = {}",
                                       max_residual_before_solve,
                                       "",
                                       count,
                                       max_residual);
                    }
                    else
                    {
                        MJRFD_LIB_INFO("Time = {}\tMax res = {}\t[{:x<{}}]\tMax res = {}",
                                       time_,
                                       max_residual_before_solve,
                                       "",
                                       count,
                                       max_residual);
                    }
                }
                else
                {
                    if(count == 1)
                    {
                        MJRFD_LIB_INFO("Newton method converged "
                                       "(with maximum residual {}) "
                                       "in {} iteration", max_residual, count);
                    }
                    else
                    {
                        MJRFD_LIB_INFO("Newton method converged "
                                       "(with maximum residual {}) "
                                       "in {} iterations", max_residual, count);
                    }

                }
            }

            if(count >= Max_newton_iterations)
            {
                // the maximum number of newton steps has been exceeded so stop
                stop = true;
                MJRFD_LIB_ERROR("The maximum number of Newton iterations "
                                "({}) has been exceeded", Max_newton_iterations);

                // Exit if the flag is set
                if(exit_on_solve_fail_ == true)
                {
                    MJRFD_LIB_FATAL("Exiting");
                    std::exit(1);
                }
            }
        }

        actions_after_solve();
    }

    void Problem::steady_solve()
    {
        // check if time derivatives are enabled
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

    void Problem::unsteady_solve(const double dt)
    {
        dt_ = dt;
        time_ += dt;

        // only output time here if using verbose logging
        // for terse logging, it's output on the same line as the solve info,
        // so it's done in solve()
        if(terse_logging_ == false)
        {
            MJRFD_LIB_INFO("Solving at time = {}", time_);
        }

        actions_before_timestep();
        Problem::solve();
        actions_after_timestep();
    }

    bool Problem::is_steady() const
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

    double& Problem::time()
    {
        return time_;
    }

    double Problem::time() const
    {
        return time_;
    }

    void Problem::output(std::ostream &) const
    {
    }

    void Problem::output_exact(std::ostream &) const
    {
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
            MJRFD_LIB_INFO("Problem::dump_res_and_jac - outputting zero");
            jac_stream << n_dof_-1 << " " << n_dof_-1 << " 0.0\n";
        }

        //jac_stream << Eigen::MatrixXd(jacobian_).format(plain_fmt);
    }

    double Problem::u_aux(const unsigned t, const unsigned i) const
    {
        assert(t < n_previous_values_+1);
        assert(i < n_aux_dof_);

        unsigned n_nodal_dof = n_var_*n_dof_per_var_;
        return u_[t](n_nodal_dof + i);
    }

    double Problem::u_aux(const unsigned i) const
    {
        assert(i < n_aux_dof_);

        unsigned n_nodal_dof = n_var_*n_dof_per_var_;
        return u_[0](n_nodal_dof + i);
    }

    double& Problem::u_aux(const unsigned t, const unsigned i)
    {
        assert(t < n_previous_values_+1);
        assert(i < n_aux_dof_);

        unsigned n_nodal_dof = n_var_*n_dof_per_var_;
        return u_[t](n_nodal_dof + i);
    }

    double& Problem::u_aux(const unsigned i)
    {
        assert(i < n_aux_dof_);

        unsigned n_nodal_dof = n_var_*n_dof_per_var_;
        return u_[0](n_nodal_dof + i);
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

    void Problem::enable_dump_jacobian(const std::string &filename_prefix)
    {
        jac_filename_prefix_ = filename_prefix;
        dump_jacobian_ = true;
    }

    void Problem::disable_dump_jacobian()
    {
        dump_jacobian_ = false;
    }

    void Problem::enable_exit_on_solve_fail()
    {
        exit_on_solve_fail_ = true;
    }

    void Problem::disable_exit_on_solve_fail()
    {
        exit_on_solve_fail_ = false;
    }

    void Problem::clear_solution()
    {
        for(unsigned t = 0; t < n_previous_values_+1; ++t)
        {
            u_[t].setZero();
        }
    }

    void Problem::calculate_jacobian(std::vector<Triplet> &) const
    {
    }

    // Default implementation which calculates the jacobian by
    // finite-differencing the residuals
    void Problem::calculate_jacobian_fd(std::vector<Triplet> &triplet_list)
    {
        triplet_list.reserve(Fd_jacobian_reserved_entries);

        // Storage for the residuals after incrementing/decrementing a dof
        Eigen::VectorXd residual_plus(n_dof_);
        Eigen::VectorXd residual_minus(n_dof_);

        // Derivative wrt dof i
        for(unsigned i = 0; i < n_dof_; ++i)
        {
            // Backup the i-th dof
            double u_i_backup = u_[0](i);

            // Calculate plus residual, setting to zero first
            u_[0](i) += Jacobian_fd_step;
            residual_plus.setZero();
            calculate_residual(residual_plus);
            u_[0](i) = u_i_backup;

            // Calculate minus residual, setting to zero first
            u_[0](i) -= Jacobian_fd_step;
            residual_minus.setZero();
            calculate_residual(residual_minus);
            u_[0](i) = u_i_backup;

            // Calculate the column of the jacobian matrix
            Eigen::VectorXd residual_column =
                (residual_plus - residual_minus)/(2.0*Jacobian_fd_step);

            // Loop over the rows of the current column
            for(unsigned j = 0; j < n_dof_; ++j)
            {
                // If the current entry is above the threshold, add it to the
                // triplet list
                if(std::abs(residual_column(j)) > Fd_jacobian_threshold)
                {
                    triplet_list.emplace_back(j, i, residual_column(j));
                }
            }
        }
    }

    void Problem::actions_before_timestep()
    {
    }

    void Problem::actions_after_timestep()
    {
    }

    void Problem::actions_before_solve()
    {
    }

    void Problem::actions_after_solve()
    {
    }

    void Problem::linear_solve()
    {
        linear_solver_.compute(jacobian_);

        if(linear_solver_.info() != Eigen::Success)
        {
            MJRFD_LIB_ERROR("Problem::linear_solve() - Compute failed!");

            if(exit_on_solve_fail_ == true)
            {
                MJRFD_LIB_FATAL("Exiting");
                std::exit(1);
            }
        }

        // get the update vector as a single array
        du_ = linear_solver_.solve(-residual_);
    }
}
