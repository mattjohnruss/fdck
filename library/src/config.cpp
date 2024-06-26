#include "config.h"
#include "log.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace fdck
{
    Config::Config(int argc, char **argv)
    {
        parse_command_line(argc, argv);
    }

    Config::Config(std::istream &is)
    {
        parse_config_file(is);
    }

    void Config::parse_config_file(std::istream &is)
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
                trim_whitespace(key);

                // storage for the value assoicated with key
                std::string value;

                // put the rest of the line into value
                if(std::getline(line_stream, value))
                {
                    // remove whitespace from value
                    trim_whitespace(value);

                    // add the (key,value) pair to the params_ map
                    params_[key] = value;
                }
            }
        }
    }

    void Config::parse_command_line(int argc, char **argv)
    {
        // Find the index of the first command line argument that starts with "--"
        // Arguments without "--" can be used directly by the driver
        // program and are ignored by the Config class. Note that an argument
        // that starts with just "-" could be a negative number and is not
        // caught here.
        unsigned first_param_index = 0;
        for(int i = 1; i < argc; ++i)
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
        for(int i = first_param_index; i < argc; i += 2)
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

    void Config::print_all() const
    {
        FDCK_CONFIG_INFO("Recognised parameters:");

        // find the length of the longest key
        std::size_t max_key_length = 0;

        for(const auto& p : params_)
        {
            if(p.first.size() > max_key_length)
            {
                max_key_length = p.first.size();
            }
        }

        // print all the keys and values (nested placeholder {} for max_key_length)
        for(const auto& [key, value] : params_)
        {
            FDCK_CONFIG_INFO("{:>{}} = {}", key, max_key_length, value);
        }
    }

    void Config::trim_whitespace(std::string &s)
    {
        trim_left(s);
        trim_right(s);
    }

    void Config::trim_left(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(),
                                        s.end(),
                                        [](int c) { return !std::isspace(c); }));
    }

    void Config::trim_right(std::string &s)
    {
        s.erase(std::find_if(s.rbegin(),
                             s.rend(),
                             [](int c) { return !std::isspace(c); }).base(),
                s.end());
    }

    /// Specialise get_impl for double
    template<>
    double Config::get_impl(const std::string &key)
    {
        // convert value to a double and return it
        return std::atof(params_[key].c_str());
    }

    /// Specialise get_impl for unsigned
    template<>
    unsigned Config::get_impl(const std::string &key)
    {
        return std::atoi(params_[key].c_str());
    }

    /// Specialise get_impl for bool
    template<>
    bool Config::get_impl(const std::string &key)
    {
        // convert value to a bool and return it
        return static_cast<bool>(std::atoi(params_[key].c_str()));
    }

    /// Specialise get_impl for string
    template<>
    std::string Config::get_impl(const std::string &key)
    {
        // convert value to a bool and return it
        return params_[key];
    }

    /// Specialise get_impl for char. We require char values to be wrapped in
    /// single quotes to avoid disappearing spaces etc. when stripping whitespace
    template<>
    char Config::get_impl(const std::string &key)
    {
        if(params_[key].size() != 3 && params_[key][0] == '\'' && params_[key][2] == '\'')
        {
            FDCK_LIB_WARN("Expected single character wrapped in single quotes"
                          "for parameter \"{}\"; found \"{}\"",
                          key, params_[key]);
        }

        // return the second character, inside the single quotes
        return params_[key][1];
    }
}
