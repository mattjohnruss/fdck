#pragma once

#include <Eigen/Core>

#include <iostream>

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
    }
}