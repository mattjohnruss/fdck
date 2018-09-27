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
    }
}
