#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace mjrfd
{
    // Simple configuration file loader and parser. Currently supports
    // ini-style key-value pairs, ignoring whitespace within each line, e.g.:
    // a=2.4
    // b = 4.134
    class Config
    {
    public:
        Config()
        {
        }

        // Constructor to immediately parse the command line arguments
        Config(int argc, char **argv)
        {
            parse_command_line(argc, argv);
        }

        // Constructor to immediately parse an istream of the config file
        Config(std::istream &is)
        {
            parse_config_file(is);
        }

        // Delete default copy constructor and assignment operator to prevent copying
        Config(const Config &) = delete;
        Config& operator=(const Config&) = delete;

        // Empty destructor
        ~Config()
        {
        }

        // Get the value associated with key from the parameters, interpreted
        // as templated type T. Must be specialised for each required type.
        template<class T>
        const T get(const std::string &key)
        {
            // check the key is actually in the map
            assert(params_.count(key) == 1);

            // Call the actual conversion and return the result
            return get_impl<T>(key);
        }

        // Read the config file (or any istream) and parse its contents,
        // storing the params in a map
        void parse_config_file(std::istream &is)
        {
            // storage for a line of the config file
            std::string line;

            // get each line of the config file
            while(std::getline(is, line))
            {
                // put the into a stream
                std::istringstream line_stream(line);

                // storage for a key
                std::string key;

                // if the line contains an '=', put the part before into key
                if(std::getline(line_stream, key, '='))
                {
                    // remove whitespace from key
                    remove_whitespace(key);

                    // storage for the value assoicated with key
                    std::string value;

                    // put the rest of the line into value
                    if(std::getline(line_stream, value))
                    {
                        // remove whitespace from value
                        remove_whitespace(value);

                        // add the (key,value) pair to the params_ map
                        params_[key] = value;
                    }
                }
            }
        }

        void parse_command_line(int argc, char **argv)
        {
            // Find the index of the first command line argument that starts with "--"
            // Arguments without "--" can be used directly by the driver
            // program and are ignored by the Config class
            unsigned first_param_index = 0;
            for(unsigned i = 1; i < argc; ++i)
            {
                std::string arg(argv[i]);

                if(arg.size() >= 2 && arg[0] == '-' && arg[1] == '-')
                {
                    first_param_index = i;
                    break;
                }
            }

            // return immediately if no command line arguments were found
            if(first_param_index == 0)
            {
                return;
            }

            // we assume that every flag has an argument, so that the number of
            // parameters to parse is argc/2 where argc is even
            assert((argc - first_param_index)%2 == 0 && "Not every parameter has an argument!");

            // loop over the supplied parameters
            for(unsigned i = first_param_index; i < argc; i += 2)
            {
                // get the key and value as strings
                std::string key = argv[i];
                std::string value = argv[i+1];

                // check if the parameter key starts with "--"
                assert(key[0] == '-' && key[1] == '-' && "Each parameter must begin with \"--\"");

                // add the (key,value) pair to the params_ map, ignoring the
                // initial "--" of each key
                params_[key.substr(2)] = value;
            }
        }

        // Print all key-value pairs to stream out
        void print_all(std::ostream &out = std::cout) const
        {
            out << "Recognised parameters:\n";

            // find the length of the longest key
            std::size_t max_key_length = 0;

            for(const auto& p : params_)
            {
                if(p.first.size() > max_key_length)
                {
                    max_key_length = p.first.size();
                }
            }

            // print all the keys and values
            for(const auto& [key, value] : params_)
            {
                out << std::setw(max_key_length) << std::left << key << " = " << value << '\n';
            }
        }

    private:
        // Map parameter names to values
        std::unordered_map<std::string, std::string> params_;

        // Helper function for removing all whitespace characters from a string
        static void remove_whitespace(std::string &s)
        {
            s.erase(std::remove_if(s.begin(),
                                   s.end(),
                                   [](unsigned char c) { return std::isspace(c); }),
                    s.end());
        }

        // Function to perform actual conversion to templated type - must be
        // specialised for each type required
        template<class T>
        const T get_impl(const std::string &key);
    };

    // Specialise get_impl for double
    template<>
    const double Config::get_impl(const std::string &key)
    {
        // convert value to a double and return it
        return std::atof(params_[key].c_str());
    }

    // Specialise get_impl for bool
    template<>
    const bool Config::get_impl(const std::string &key)
    {
        // convert value to a bool and return it
        return static_cast<bool>(std::atoi(params_[key].c_str()));
    }

    // Specialise get_impl for string
    template<>
    const std::string Config::get_impl(const std::string &key)
    {
        // convert value to a bool and return it
        return params_[key];
    }
}
