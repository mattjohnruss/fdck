#pragma once

#include <Eigen/Core>

#include <iostream>
#include <vector>

namespace mjrfd
{
    namespace utilities
    {
        /// Integrate an array of function values sampled at intervals dx using
        /// the trapezium rule
        double integrate_trapezium(const Eigen::VectorXd &v, const double dx);

        /// Integrate a function f on the inteval [a,b] using the trapezium
        /// rule with n points
        template<class F>
        double integrate_trapezium(const F &f,
                                   const double a,
                                   const double b,
                                   const unsigned n)
        {
            assert(n >= 3);
            assert(b > a);

            Eigen::VectorXd v(n);

            const double dx = (b - a)/static_cast<double>(n - 1);

            for(unsigned i = 0; i < n; ++i)
            {
                const double x = a + i*dx;
                v(i) = f(x);
            }

            return integrate_trapezium(v, dx);
        }

        double l2_norm(const Eigen::VectorXd &v, const double dx);

        template<class F>
        double l2_norm(const F &f,
                       const double a,
                       const double b,
                       const unsigned n)
        {
            auto f_sq = [&](double x) { return std::pow(f(x), 2); };
            return std::sqrt(integrate_trapezium(f_sq, a, b, n));
        }

        // Evaluate a polynomial with coeffs in ascending power order at x using
        // Horner's method
        double evaluate_polynomial(const double x, const std::vector<double> &coeffs);

        double lerp(double s, double v0, double v1);

        void lerp_mesh(const Eigen::VectorXd &x,
                       const Eigen::VectorXd &v,
                       const Eigen::VectorXd &x_interp,
                       Eigen::VectorXd &v_interp);
    }
}
