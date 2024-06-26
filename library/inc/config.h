#pragma once

#include <cassert>
#include <iostream>
#include <istream>
#include <string>
#include <unordered_map>

namespace fdck
{
    /// Simple command-line arguments and configuration file loader and parser.
    /// Items are added to a map in the order they are found. This allows to
    /// override the values in the config file with the cmdline, say.
    /// Currently supports long-style cmdline args, e.g. "--key value", and
    /// ini-style key-value pairs, ignoring whitespace within each line, e.g.:
    /// a=2.4
    /// b = 4.134
    class Config
    {
    public:
        /// Default constructor
        Config() = default;

        /// Constructor to immediately parse the command line arguments
        Config(int argc, char **argv);

        /// Constructor to immediately parse an istream of the config file
        explicit Config(std::istream &is);

        /// Delete default copy constructor and assignment operator to prevent copying
        Config(const Config &) = delete;
        Config& operator=(const Config&) = delete;

        /// Get the value associated with key from the parameters, interpreted
        /// as templated type T. Must be specialised for each required type.
        template<class T>
        T get(const std::string &key)
        {
            // check the key is actually in the map
            assert(params_.count(key) == 1);

            // Call the actual conversion and return the result
            return get_impl<T>(key);
        }

        /// Get the value associated with key from the parameters as with
        /// get(), unless the key does not exist in which case the default
        /// value def is returned
        template<class T>
        T get_or(const std::string &key, const T &def)
        {
            unsigned count = params_.count(key);

            if(count == 1)
            {
                return get_impl<T>(key);
            }
            else // count must be 0
            {
                return def;
            }
        }

        /// Read the config file (or any istream) and parse its contents,
        /// adding the params to the map
        void parse_config_file(std::istream &is);

        /// Read and parse the command-line arguments adding the params to the
        /// map. Currently limited to simple "--key value" pairs of long-style
        /// options which always have values. Ignores items in argv until it
        /// finds an argument starting with "--", but not just "-" to allow for
        /// negative numbers, (allows drivers to have custom positional
        /// arguments without much hassle here), and then all arguments are
        /// assumed to be of the above format
        void parse_command_line(int argc, char **argv);

        /// Print all key-value pairs to the default config_logger
        void print_all() const;

    private:
        /// Map parameter names to values
        std::unordered_map<std::string, std::string> params_;

        /// Helper function for removing all whitespace characters from a string
        static void trim_whitespace(std::string &s);
        static void trim_left(std::string &s);
        static void trim_right(std::string &s);

        /// Function to perform actual conversion to templated type - must be
        /// specialised for each type required
        template<class T>
        T get_impl(const std::string &key);
    };
}
