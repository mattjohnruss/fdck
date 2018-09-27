#include <functions.h>
#include <utilities.h>

#include <catch2/catch.hpp>

namespace mjrfd
{
    TEST_CASE( "Numerical integration using trapezium rule", "[integration]" )
    {
        using utilities::integrate_trapezium;

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
}
