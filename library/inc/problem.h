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
        typedef Eigen::Triplet<double> Triplet;

        Problem(const unsigned n_var,
                const unsigned n_dof_per_var,
                const unsigned n_aux_dof = 0);

        virtual ~Problem();

        void solve();
        void steady_solve();
        void unsteady_solve(const double dt);

        bool is_steady() const;
        virtual void make_steady();
        virtual void make_unsteady();

        double& time();
        double time() const;

        virtual void output(std::ostream &out) const = 0;
        virtual void output_exact(std::ostream &out) const = 0;

        void dump_res_and_jac(std::ostream &res_stream, std::ostream &jac_stream) const;

        template<class T>
        double u(const T variable, const unsigned i) const
        {
            unsigned v = static_cast<unsigned>(variable);
            assert(v >= 0 && v < n_var_);
            assert(i >= 0 && i < n_dof_per_var_);

            return u_[0](v*n_dof_per_var_ + i);
        }

        template<class T>
        double u(const unsigned t, const T variable, const unsigned i) const
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

        double u_aux(const unsigned t, const unsigned i) const;
        double u_aux(const unsigned i) const;

        double& u_aux(const unsigned t, const unsigned i);
        double& u_aux(const unsigned i);

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
        const unsigned n_aux_dof_;

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

        virtual void calculate_residual(Eigen::VectorXd &residual) const = 0;
        virtual void calculate_jacobian(std::vector<Triplet> &triplet_list) const;

        void calculate_jacobian_fd(Eigen::SparseMatrix<double> &jacobian);

        virtual void actions_before_timestep();
        virtual void actions_after_timestep();

        virtual void actions_before_solve();
        virtual void actions_after_solve();

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
