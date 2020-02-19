#include "functions.h"
#include "utilities.h"

#include <catch2/catch.hpp>

namespace fdck
{
    TEST_CASE( "Numerical integration using trapezium rule", "[integration]" )
    {
        using utilities::integrate_trapezium;

        SECTION( "Array of function values" )
        {
            double dx = 0.15;
            Eigen::VectorXd v(8);
            v << 0.2, 0.5, 0.6, 0.9, 1.6, 3.1, 2.7, 2.4;

            REQUIRE( integrate_trapezium(v, dx) == Approx(1.605) );
        }

        SECTION( "Define functions with lambdas" )
        {
            SECTION( "Constant function" )
            {
                REQUIRE( integrate_trapezium([] (double) { return 1.0; }, 0.0, 1.0, 101) == Approx(1.0) );
            }

            SECTION( "Linear function" )
            {
                REQUIRE( integrate_trapezium([] (double x) { return x; }, 0.0, 1.0, 101) == Approx(0.5) );
            }
        }

        SECTION( "Define functions with DifferentiableFunction" )
        {
            SECTION( "Constant function" )
            {
                ConstantFunction f(1.0);
                REQUIRE( integrate_trapezium(f, 0.0, 1.0, 101) == Approx(1.0) );
            }

            SECTION( "Linear function" )
            {
                class LinearFunction : public DifferentiableFunction
                {
                public:
                    LinearFunction(const double a_0, const double a_1) : a_0_(a_0), a_1_(a_1) {}

                    double value(const double c) const override { return a_0_ + a_1_*c; }
                    double deriv(const double) const override { return a_1_; }

                private:
                    double a_0_;
                    double a_1_;
                };

                LinearFunction f(3.0, 2.0);
                REQUIRE( integrate_trapezium(f, -5.0, 1.0, 11) == Approx(-6.0) );
            }
        }
    }

    TEST_CASE( "Numerical compuatation of L2 norms" )
    {
        using utilities::l2_norm;

        SECTION( "Define functions with lambdas" )
        {
            auto f = [](double x) { return 1.0/3.0 - 2.0*x + 0.5*std::pow(x, 2); };
            REQUIRE( l2_norm(f, 0.0, 1.0, 1001) == Approx(std::sqrt(15.8)/6.0) );
        }

        SECTION( "Array of function values" )
        {
            double dx = 0.15;
            Eigen::VectorXd v(8);
            v << 0.2, 0.5, 0.6, 0.9, 1.6, 3.1, 2.7, 2.4;

            REQUIRE( l2_norm(v, dx) == Approx(1.888650311730576) );
        }
    }
}
