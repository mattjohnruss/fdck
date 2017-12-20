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
        Problem(unsigned n_dof);
        virtual ~Problem();

        void solve();
        virtual void output(std::ostream &out) const = 0;

        static double Max_residual;
    protected:
        Eigen::VectorXd u_;
        Eigen::VectorXd residual_;

        //Eigen::MatrixXd jacobian_;
        Eigen::SparseMatrix<double> jacobian_;

        virtual void calculate_residual() = 0;
        virtual void calculate_jacobian() = 0;
        
    private:
        void linear_solve();

        Eigen::VectorXd du_;
        Eigen::SparseLU<Eigen::SparseMatrix<double> > linear_solver_;

        //unsigned n_dof_;

        //std::vector<Residual*> m_residuals;
    };
}
