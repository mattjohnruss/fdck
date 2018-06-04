#pragma once

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <vector>
#include <ostream>

#include <residual.h>

namespace mjrfd
{
    class Problem
    {
    public:
        Problem(const unsigned n_var,
                const unsigned n_dof_per_var,
                const double dt = 0.0);

        virtual ~Problem();

        void solve();
        void steady_solve();
        void unsteady_solve();

        const bool is_steady() const;
        virtual void make_steady();
        virtual void make_unsteady();

        double& time();
        const double time() const;

        virtual void output(std::ostream &out) const = 0;
        virtual void output_exact(std::ostream &out) const = 0;

        void dump_res_and_jac(std::ostream &res_stream, std::ostream &jac_stream) const;

        template<class T>
        const double u(const T variable, const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            assert(v >= 0 && v < n_var_);
            assert(i >= 0 && i < n_dof_per_var_);

            return u_[0](v*n_dof_per_var_ + i);
        }

        template<class T>
        const double u(const unsigned t, const T variable, const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            assert(t >= 0 && t < n_time_history_);
            assert(v >= 0 && v < n_var_);
            assert(i >= 0 && i < n_dof_per_var_);

            return u_[t](v*n_dof_per_var_ + i);
        }

        template<class T>
        double& u(const T variable, const unsigned i)
        {
            unsigned v = static_cast<unsigned>(variable);
            assert(v >= 0 && v < n_var_);
            assert(i >= 0 && i < n_dof_per_var_);

            return u_[0](v*n_dof_per_var_ + i);
        }

        template<class T>
        double& u(const unsigned t, const T variable, const unsigned i)
        {
            unsigned v = static_cast<unsigned>(variable);
            assert(t >= 0 && t < n_time_history_);
            assert(v >= 0 && v < n_var_);
            assert(i >= 0 && i < n_dof_per_var_);

            return u_[t](v*n_dof_per_var_ + i);
        }

        void enable_terse_logging();
        void disable_terse_logging();

        void enable_fd_jacobian();
        void disable_fd_jacobian();

        void enable_dump_jacobian(const std::string &filename_prefix = "");
        void disable_dump_jacobian();

        void enable_exit_on_solve_fail();
        void disable_exit_on_solve_fail();

        void clear_solution();

    protected:
        const unsigned n_dof_per_var_;
        const unsigned n_dof_;
        const unsigned n_var_;

        const unsigned n_time_history_;

        std::vector<Eigen::VectorXd> u_;

    private:
        Eigen::VectorXd du_;

    protected:
        Eigen::VectorXd residual_;
        Eigen::SparseMatrix<double> jacobian_;

        double time_;
        double dt_;

        bool steady_;

        // TODO rewrite so that these accept a vector/matrix/vector of triplets etc by reference
        virtual void calculate_residual() = 0;
        virtual void calculate_jacobian() = 0;

        double Max_residual;
        unsigned Max_newton_iterations;
        double Jacobian_fd_step;

    private:
        void linear_solve();

        Eigen::SparseLU<Eigen::SparseMatrix<double> > linear_solver_;

        bool terse_logging_;
        bool use_fd_jacobian_;
        bool dump_jacobian_;
        bool exit_on_solve_fail_;

        std::string jac_filename_prefix_;
    };
}
