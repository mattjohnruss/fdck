#include "log.h"

#include <spdlog/sinks/stdout_color_sinks.h>

namespace fdck
{
    namespace log
    {
        // the single instantiations of the loggers
        std::shared_ptr<spdlog::logger> library_logger = spdlog::stdout_color_st("library");
        std::shared_ptr<spdlog::logger> config_logger = spdlog::stdout_color_st("config");
        std::shared_ptr<spdlog::logger> driver_logger = spdlog::stdout_color_st("driver");

        void enable_labels()
        {
            spdlog::set_pattern("[%n] [%^%l%$] %v");
        }

        void disable_labels()
        {
            spdlog::set_pattern("%v");
        }

        void set_level(const std::string &level)
        {
            if(level == "trace")
                spdlog::set_level(spdlog::level::trace);
            else if(level == "debug")
                spdlog::set_level(spdlog::level::debug);
            else if(level == "info")
                spdlog::set_level(spdlog::level::info);
            else if(level == "warn")
                spdlog::set_level(spdlog::level::warn);
            else if(level == "error")
                spdlog::set_level(spdlog::level::err);
            else if(level == "fatal")
                spdlog::set_level(spdlog::level::critical);
            else
            {
                FDCK_LIB_ERROR("Log level \"{}\" not recognised!", level);
            }
        }
    }
}
