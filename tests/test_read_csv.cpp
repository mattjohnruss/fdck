#include "utilities.h"
#include "log.h"

#include <catch2/catch.hpp>

#include <fstream>

namespace fdck
{
    TEST_CASE( "Read and parse CSV files", "[csv]" )
    {
        std::vector<double> t = { 0.1, 0.3, 0.4, 0.7, 0.9, 1.5, 2.1 };
        std::vector<double> c = { 2.4, 3.2, 3.1, 3.4, 4.7, 9.1, 23.4 };

        SECTION( "CSV file with space delimiter" )
        {
            std::ifstream csv_file("../../tests/test_read_csv.csv");

            REQUIRE( csv_file.is_open() == true );

            auto [m_vec, n_rows, n_cols] =
                fdck::utilities::read_csv_to_flat_vector(csv_file, ' ');

            REQUIRE( n_rows == 7 );
            REQUIRE( n_cols == 2 );
            REQUIRE( m_vec.size() == n_rows * n_cols );

            typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> S;

            Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, S> m(m_vec.data(), n_rows, n_cols, S(1, n_cols));

            for(unsigned row = 0; row < n_rows; ++row)
            {
                REQUIRE( m(row, 0) == Approx( t[row] ) );
                REQUIRE( m(row, 1) == Approx( c[row] ) );
            }
        }

        SECTION( "CSV file with space delimiter and header" )
        {
            std::ifstream csv_file("../../tests/test_read_csv_with_header.csv");

            REQUIRE( csv_file.is_open() == true );

            unsigned skip_rows = 1;

            auto [m_vec, n_rows, n_cols] =
                fdck::utilities::read_csv_to_flat_vector(csv_file, ' ', skip_rows);

            REQUIRE( n_rows == 7 );
            REQUIRE( n_cols == 2 );
            REQUIRE( m_vec.size() == n_rows * n_cols );

            typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> S;

            Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, S> m(m_vec.data(), n_rows, n_cols, S(1, n_cols));

            for(unsigned row = 0; row < n_rows; ++row)
            {
                REQUIRE( m(row, 0) == Approx( t[row] ) );
                REQUIRE( m(row, 1) == Approx( c[row] ) );
            }
        }

        SECTION( "CSV file with comma delimiter" )
        {
            std::ifstream csv_file("../../tests/test_read_csv_comma.csv");

            REQUIRE( csv_file.is_open() == true );

            auto [m_vec, n_rows, n_cols] =
                fdck::utilities::read_csv_to_flat_vector(csv_file, ',');

            REQUIRE( n_rows == 7 );
            REQUIRE( n_cols == 2 );
            REQUIRE( m_vec.size() == n_rows * n_cols );

            typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> S;

            Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, S> m(m_vec.data(), n_rows, n_cols, S(1, n_cols));

            for(unsigned row = 0; row < n_rows; ++row)
            {
                REQUIRE( m(row, 0) == Approx( t[row] ) );
                REQUIRE( m(row, 1) == Approx( c[row] ) );
            }
        }
    }
}
