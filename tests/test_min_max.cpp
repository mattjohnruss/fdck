#include "utilities.h"

#include <catch2/catch.hpp>
#include <limits>

namespace fdck
{
    TEST_CASE( "Minimum and all negative", "[min]" )
    {
        using utilities::min;
        using utilities::all_negative;

        // can't do this because catch2 seems unable to handle the abort signal
        // issued by assert
        //SECTION( "Empty vector" )
        //{
            //REQUIRE_THROWS( min<double>({}) );
        //}

        SECTION( "Regular values via initialiser list" )
        {
            const auto result = min<double>({ -0.1, -100.3, 0.0, 1.0, 235.1 });
            REQUIRE( result == Approx(-100.3) );
        }

        SECTION( "Mixture of regular and limiting values" )
        {
            auto result = min<double>({ -0.1, -100.3, 0.0, 1.0, 235.1, std::numeric_limits<double>::max() });
            REQUIRE( result == Approx(-100.3) );

            result = min<double>({ -0.1, -100.3, 0.0, 1.0, 235.1, std::numeric_limits<double>::min() });
            REQUIRE( result == Approx(-100.3) );

            result = min<double>({ -0.1, -100.3, 0.0, 1.0, 235.1, std::numeric_limits<double>::lowest() });
            REQUIRE( result == Approx(std::numeric_limits<double>::lowest()) );
        }

        SECTION( "all_negative(...)" )
        {
            bool result = all_negative<double>({ -10.3, -21.1 -1004.3111, -0.001, -1.0e-15 });
            REQUIRE( result == true );

            result = all_negative<double>({ -10.3, -21.1 -1004.3111, -0.001, -1.0e-15, 0.0 });
            REQUIRE( result == false );

            result = all_negative<double>({ -10.3, -21.1 -1004.3111, -0.001, -1.0e-15, 2.5 });
            REQUIRE( result == false );

            result = all_negative<double>({ -10.3, -21.1 -1004.3111, -0.001, -1.0e-15, std::numeric_limits<double>::lowest() });
            REQUIRE( result == true );

            result = all_negative<double>({ -10.3, -21.1 -1004.3111, -0.001, -1.0e-15, std::numeric_limits<double>::min() });
            REQUIRE( result == false );
        }
    }

    TEST_CASE( "Maximum and all positive", "[max]" )
    {
        using utilities::max;
        using utilities::all_positive;

        SECTION( "Regular values via initialiser list" )
        {
            const auto result = max<double>({ -0.1, -100.3, 0.0, 1.0, 235.1 });
            REQUIRE( result == Approx(235.1) );
        }

        SECTION( "Mixture of regular and limiting values" )
        {
            auto result = max<double>({ -0.1, -100.3, 0.0, 1.0, 235.1, std::numeric_limits<double>::max() });
            REQUIRE( result == Approx(std::numeric_limits<double>::max()) );

            result = max<double>({ -0.1, -100.3, 0.0, 1.0, 235.1, std::numeric_limits<double>::lowest() });
            REQUIRE( result == Approx(235.1) );
        }

        SECTION( "all_positive(...)" )
        {
            bool result = all_positive<double>({ 1.0, 3.0, 1.0e15, 1.0e-15, 0.12351231 });
            REQUIRE( result == true );

            result = all_positive<double>({ 1.0, 3.0, 1.0e15, 1.0e-15, 0.12351231, 0.0 });
            REQUIRE( result == false );

            result = all_positive<double>({ 1.0, 3.0, 1.0e15, 1.0e-15, 0.12351231, -2.5 });
            REQUIRE( result == false );

            result = all_positive<double>({ 1.0, 3.0, 1.0e15, 1.0e-15, 0.12351231, std::numeric_limits<double>::max() });
            REQUIRE( result == true );

            result = all_positive<double>({ 1.0, 3.0, 1.0e15, 1.0e-15, 0.12351231, std::numeric_limits<double>::min() });
            REQUIRE( result == true );
        }
    }
}
