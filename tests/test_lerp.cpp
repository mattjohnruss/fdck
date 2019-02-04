#include "utilities.h"

#include <catch2/catch.hpp>
#include <Eigen/Dense>

namespace mjrfd
{
    TEST_CASE( "Linear interpolation" )
    {
        SECTION( "Simple lerp" )
        {
            double v0 = 3.5;
            double v1 = 4.0;
            double s = 0.5;
            REQUIRE( utilities::lerp(s, v0, v1) == Approx(3.75) );
        }

        SECTION( "Lerp between meshes" )
        {
            double x_max = 7.5;
            unsigned n = 123;
            unsigned n_interp = 11;

            // original grid
            Eigen::VectorXd x(n);
            x.setLinSpaced(n, 0.0, x_max);

            // original data sampled on x
            Eigen::VectorXd v(n);

            const double pi = 2.0*std::acos(0.0);

            // assign a sine wave
            for(unsigned i = 0; i < n; ++i)
            {
                v(i) = std::sin(2.0*pi*x(i)/x_max);
            }

            // grid onto which we interpolate
            Eigen::VectorXd x_interp(n_interp);
            x_interp.setLinSpaced(n_interp, 0.0, x_max);

            // vector the interpolated data will be stored in
            Eigen::VectorXd v_interp(n_interp);

            // perform the interpolation
            mjrfd::utilities::lerp_mesh(x, v, x_interp, v_interp);

            // hardcode the manually checked results
            Eigen::VectorXd x_interp_reference(n_interp);
            x_interp_reference << 0, 0.75, 1.5, 2.25, 3, 3.75, 4.5, 5.25, 6, 6.75, 7.5;

            Eigen::VectorXd v_interp_reference(n_interp);
            v_interp_reference << 0, 0.587659, 0.950753, 0.950753, 0.587659, 5.66554e-16, -0.587659, -0.950753, -0.950753, -0.587659, -1.13311e-15;

            // check the computed results match the hardcoded ones
            for(unsigned i = 0; i < n_interp; ++i)
            {
                REQUIRE( x_interp(i) == Approx(x_interp_reference(i)) );
                REQUIRE( v_interp(i) == Approx(v_interp_reference(i)) );
            }
        }

        SECTION( "Lerp between meshes (using Eigen::Map)" )
        {
            double x_max = 7.5;
            unsigned n = 123;
            unsigned n_interp = 11;

            // original grid
            Eigen::VectorXd x(n);
            x.setLinSpaced(n, 0.0, x_max);

            // original data sampled on x
            std::vector<double> v(n);

            const double pi = 2.0*std::acos(0.0);

            // assign a sine wave
            for(unsigned i = 0; i < n; ++i)
            {
                v[i] = std::sin(2.0*pi*x(i)/x_max);
            }

            // make an Eigen::Map mapping the data from the underlying vector
            Eigen::Map<Eigen::VectorXd> v_map(v.data(), n);

            // grid onto which we interpolate
            Eigen::VectorXd x_interp(n_interp);
            x_interp.setLinSpaced(n_interp, 0.0, x_max);

            // vector the interpolated data will be stored in
            std::vector<double> v_interp(n_interp);

            // make an Eigen::Map mapping the data from the underlying vector
            Eigen::Map<Eigen::VectorXd> v_interp_map(v_interp.data(), n_interp);

            // perform the interpolation
            mjrfd::utilities::lerp_mesh(x, v_map, x_interp, v_interp_map);

            // hardcode the manually checked results
            Eigen::VectorXd x_interp_reference(n_interp);
            x_interp_reference << 0, 0.75, 1.5, 2.25, 3, 3.75, 4.5, 5.25, 6, 6.75, 7.5;

            Eigen::VectorXd v_interp_reference(n_interp);
            v_interp_reference << 0, 0.587659, 0.950753, 0.950753, 0.587659, 5.66554e-16, -0.587659, -0.950753, -0.950753, -0.587659, -1.13311e-15;

            // check the computed results match the hardcoded ones
            for(unsigned i = 0; i < n_interp; ++i)
            {
                REQUIRE( x_interp(i) == Approx(x_interp_reference(i)) );
                REQUIRE( v_interp[i] == Approx(v_interp_reference(i)) );
            }
        }
    }
}
