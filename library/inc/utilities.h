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

        /// Calculate the L2 norm of a function sampled at intervals dx using
        /// the trapezium rule
        double l2_norm(const Eigen::VectorXd &v, const double dx);

        /// Calculate the L2 norm of a function on the interval [a, b] using
        /// the trapezium rule with n points
        template<class F>
        double l2_norm(const F &f,
                       const double a,
                       const double b,
                       const unsigned n)
        {
            auto f_sq = [&](double x) { return std::pow(f(x), 2); };
            return std::sqrt(integrate_trapezium(f_sq, a, b, n));
        }

        /// Evaluate a polynomial with coeffs in ascending power order at x
        /// using Horner's method
        double evaluate_polynomial(const double x, const std::vector<double> &coeffs);

        /// Simple 1d linear interpolation
        double lerp(double s, double v0, double v1);

        /// Resample data v defined at points x at the new set of points
        /// x_interp, using linear interpolation. Templated to avoid
        /// difficulties with binding Ref<VectorXd>& to "col"s from a
        /// Map<...> that can be either row- or column-major. We only need to
        /// call e.g. x(i), v(i) etc. in the function body.
        template<class V1, class V2, class V3>
        void lerp_mesh(const V1 &x,
                       const V2 &v,
                       const V3 &x_interp,
                       Eigen::Ref<Eigen::VectorXd> v_interp)
        {
            // get the number of points in the original grid
            const unsigned n = x.size();

            // check it matches the number of data points
            assert(n == v.size());

            // get the number of points in the new grid
            const unsigned n_interp = x_interp.size();

            // check it matches the number of data points
            assert(n_interp == v_interp.size());

            // loop over the points in x_interp
            for(unsigned i = 0; i < n_interp; ++i)
            {
                // Binary search:
                // Start the search at the (approx) middle of the original grid.
                // There are n-1 cells, so the index of the last cell is n-2
                unsigned cell = 0;

                {
                    unsigned l_cell = 0;
                    unsigned r_cell = n - 2;

                    // we assume this can't fail
                    while(l_cell <= r_cell)
                    {
                        cell = (l_cell + r_cell)/2;

                        if(x_interp(i) > x(cell+1))
                        {
                            l_cell = cell + 1;
                        }
                        else if(x_interp(i) < x(cell))
                        {
                            r_cell = cell - 1;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                // After the loop, cell is the index of the interval that
                // contains x_interp(i)
                // Now do the actual lerp:
                const double s = (x_interp(i) - x(cell))/(x(cell+1) - x(cell));
                v_interp(i) = lerp(s, v(cell), v(cell+1));
            }
        }

        /// Read CSV-style data from the stream is into a row-major flat vector
        /// Returns a tuple of (data, n_rows, n_cols)
        std::tuple<std::vector<double>, unsigned, unsigned>
        read_csv_to_flat_vector(std::istream &is,
                                char delimiter = ' ',
                                unsigned skip_rows = 0);
    }
}
