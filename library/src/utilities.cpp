#include <utilities.h>

namespace mjrfd
{
    namespace utilities
    {
        double integrate_trapezium(const Eigen::VectorXd &v, const double dx)
        {
            double integral = 0.0;

            unsigned n = v.size();

            assert(n >= 3);

            integral += 0.5*dx*v(0);
            integral += (dx*v).segment(1, n - 2).sum();
            integral += 0.5*dx*v(n - 1);

            return integral;
        }

        double l2_norm(const Eigen::VectorXd &v, const double dx)
        {
            return std::sqrt(integrate_trapezium(v.array().square(), dx));
        }

        double evaluate_polynomial(const double x, const std::vector<double> &coeffs)
        {
            double result = 0.0;
            unsigned n = coeffs.size();

            for(int i = n-1; i >= 0; --i)
            {
                assert(static_cast<unsigned>(i) <= (n-1) && i >= 0);
                result = result*x + coeffs[i];
            }

            return result;
        }

        /// Simple lerp function
        double lerp(double s, double v0, double v1)
        {
            return (1.0 - s)*v0 + s*v1;
        }

        /// Linearly interpolate data v from grid x onto grid x_interp, returning
        /// v_interp
        void lerp_mesh(const Eigen::Ref<Eigen::VectorXd> &x,
                       const Eigen::Ref<Eigen::VectorXd> &v,
                       const Eigen::Ref<Eigen::VectorXd> &x_interp,
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
    }
}
