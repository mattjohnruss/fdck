#include <catch2/catch.hpp>

//#include <advection_diffusion_reaction_problem.h>

unsigned factorial(unsigned x)
{
    if(x == 0 || x == 1)
        return 1;
    else
        return factorial(x-1)*x;
}

TEST_CASE( "Factorials are computed", "[factorial]" )
{
    REQUIRE( factorial(0) == 1 );
    REQUIRE( factorial(1) == 1 );
    REQUIRE( factorial(2) == 2 );
    REQUIRE( factorial(3) == 6 );
    REQUIRE( factorial(10) == 3628800 );
}

