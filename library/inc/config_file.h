#pragma once

#include <algorithm>
#include <istream>
#include <map>
#include <sstream>
#include <string>
#include <iostream>

namespace mjrfd
{
    // Simple configuration file loader and parser. Currently supports
    // ini-style key-value pairs, ignoring whitespace within each line, e.g.:
    // a=2.4
    // b = 4.134
    class ConfigFile
    {
    public:
        // Default constructor takes an istream of the config file
        ConfigFile(std::istream &is)
        {
            read_and_parse(is);
        }

        // Delete default copy constructor and assignment operator to prevent copying
        ConfigFile(const ConfigFile &) = delete;
        ConfigFile& operator=(const ConfigFile&) = delete;

        // Empty destructor
        ~ConfigFile()
        {
        }

        // Get the value associated with key from the parameters, interpreted
        // as templated type T. Must be specialised for each required type.
        template<class T>
        const T get(const std::string &key);

        // Read the config file (or any istream) and parse its contents,
        // storing the params in a map
        void read_and_parse(std::istream &is)
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

        // Print all key-value pairs to stream out
        void print_all(std::ostream &out = std::cout) const
        {
            for(const auto& [key, value] : params_)
            {
                out << key << " = " << value << '\n';
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
    };

    // Specialise get for double
    template<>
    const double ConfigFile::get(const std::string &key)
    {
        // check the key is actually in the map
        assert(params_.count(key) == 1);

        // convert value to a double and return it
        return std::atof(params_[key].c_str());
    }
}
