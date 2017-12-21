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

        void solve();
        void unsteady_solve();

        const double time() const;

        virtual void output(std::ostream &out) const = 0;

        static double Max_residual;

    protected:
        Eigen::VectorXd u_;

        // FIXME temporary hack for timestepping
        Eigen::VectorXd u_old_;

    private:
        Eigen::VectorXd du_;

    protected:
        Eigen::VectorXd residual_;
        Eigen::SparseMatrix<double> jacobian_;

        double time_;
        double dt_;

        virtual void calculate_residual() = 0;
        virtual void calculate_jacobian() = 0;

    private:
        void linear_solve();

        Eigen::SparseLU<Eigen::SparseMatrix<double> > linear_solver_;

        //unsigned n_dof_;

        //std::vector<Residual*> m_residuals;
    };
}
