#include "utilities.h"

#include <sstream>

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

        double lerp(double s, double v0, double v1)
        {
            return (1.0 - s)*v0 + s*v1;
        }

        std::tuple<std::vector<double>, unsigned, unsigned>
        read_csv_to_flat_vector(std::istream &is,
                                char delimiter,
                                unsigned skip_rows)
        {
            std::vector<double> data;
            std::string line;

            unsigned n_rows = 0;
            unsigned n_cols = 0;

            // loop over the rows
            while(std::getline(is, line))
            {
                std::istringstream line_stream(line);
                std::string x;

                if(n_rows < skip_rows)
                {
                    n_rows += 1;
                    continue;
                }

                // loop over the columns
                while(std::getline(line_stream, x, delimiter))
                {
                    data.emplace_back(std::stod(x));
                    if(n_rows == skip_rows)
                    {
                        n_cols += 1;
                    }
                }

                n_rows += 1;
            }

            // We rely on RVO to avoid copying into the return value. It's
            // actually eliding the copy of the tuple, but we also std::move
            // the vector into the tuple, so that is never copied either
            return { std::move(data), n_rows, n_cols };
        }
    }
}
