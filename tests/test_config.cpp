#include <config.h>

#include <catch2/catch.hpp>

#include <fstream>
#include <sstream>
#include <string>

namespace mjrfd
{

    TEST_CASE( "Configuration file and command line arguments", "[config]" )
    {
        mjrfd::Config config;

        SECTION( "Configuration from file" )
        {
            std::ifstream file("../../tests/test_config.ini");
            config.parse_config_file(file);

            REQUIRE( config.get<double>("a") == 1.3 );
            REQUIRE( config.get<bool>("b") == false );
            REQUIRE( config.get<bool>("c") == true );
            REQUIRE( config.get<std::string>("d") == "AString" );

            SECTION( "Override configuration file with command line" )
            {
                int argc = 3;
                const char * const argv[] = { "test_config", "--a", "4.6" };

                config.parse_command_line(argc, const_cast<char **>(argv));

                REQUIRE( config.get<double>("a") == 4.6 );
                REQUIRE( config.get<bool>("b") == false );
                REQUIRE( config.get<bool>("c") == true );
                REQUIRE( config.get<std::string>("d") == "AString" );
            }
        }

        SECTION( "Configuration from stringstream" )
        {
            std::stringstream config_ss;
            config_ss << "a = 1.3\n"
                      << "b = 0\n"
                      << "c = 1\n"
                      << "d = AString";

            config.parse_config_file(config_ss);

            REQUIRE( config.get<double>("a") == 1.3 );
            REQUIRE( config.get<bool>("b") == false );
            REQUIRE( config.get<bool>("c") == true );
            REQUIRE( config.get<std::string>("d") == "AString" );
        }

        SECTION( "Configuration from command line arguments" )
        {
            SECTION( "Only \"--key value\" arguments" )
            {
                int argc = 9;
                const char * const argv[] = { "test_config", "--a", "1.3", "--b", "0", "--c", "1", "--d", "AString" };

                config.parse_command_line(argc, const_cast<char **>(argv));

                REQUIRE( config.get<double>("a") == 1.3 );
                REQUIRE( config.get<bool>("b") == false );
                REQUIRE( config.get<bool>("c") == true );
                REQUIRE( config.get<std::string>("d") == "AString" );

                SECTION( "Override command line with configuration file" )
                {
                    std::stringstream config_ss;
                    config_ss << "a = 4.6";

                    config.parse_config_file(config_ss);

                    REQUIRE( config.get<double>("a") == 4.6 );
                    REQUIRE( config.get<bool>("b") == false );
                    REQUIRE( config.get<bool>("c") == true );
                    REQUIRE( config.get<std::string>("d") == "AString" );
                }
            }

            SECTION( "Positional then \"--key value\" arguments" )
            {
                int argc = 12;
                const char * const argv[] = { "test_config", "should", "be", "ignored", "--a", "1.3", "--b", "0", "--c", "1", "--d", "AString" };

                config.parse_command_line(argc, const_cast<char **>(argv));

                REQUIRE( config.get<double>("a") == 1.3 );
                REQUIRE( config.get<bool>("b") == false );
                REQUIRE( config.get<bool>("c") == true );
                REQUIRE( config.get<std::string>("d") == "AString" );
            }
        }
    }
}
