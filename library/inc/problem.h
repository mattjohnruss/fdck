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
        Problem(const unsigned n_dof, const double dt = 0.0);
        virtual ~Problem();

        static double Max_residual;
        static unsigned Max_newton_iterations;

        void solve(const bool dump = false);
        void steady_solve(const bool dump = false);
        void unsteady_solve(const bool dump = false);

        const bool is_steady() const;
        virtual void make_steady();
        virtual void make_unsteady();

        const double time() const;

        virtual void output(std::ostream &out) const = 0;

        void dump_res_and_jac(std::ostream &res_stream, std::ostream &jac_stream) const;

        const double u(const unsigned i) const;
        const double u(const unsigned t, const unsigned i) const;

        double& u(const unsigned i);
        double& u(const unsigned t, const unsigned i);

        void enable_terse_logging();
        void disable_terse_logging();

    protected:
        std::vector<Eigen::VectorXd> u_;

    private:
        Eigen::VectorXd du_;

    protected:
        Eigen::VectorXd residual_;
        Eigen::SparseMatrix<double> jacobian_;

        double time_;
        double dt_;

        unsigned n_dof_;
        bool steady_;

        virtual void calculate_residual() = 0;
        virtual void calculate_jacobian() = 0;

    private:
        void linear_solve();

        Eigen::SparseLU<Eigen::SparseMatrix<double> > linear_solver_;

        bool terse_logging_;

        //std::vector<Residual*> m_residuals;
    };
}
