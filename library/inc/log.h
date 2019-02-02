#pragma once

#include <memory>

#include <spdlog/spdlog.h>

namespace mjrfd
{
    namespace log
    {
        // declare loggers as extern as the (unique) definitions will appear in
        // the .cpp file
        extern std::shared_ptr<spdlog::logger> library_logger;
        extern std::shared_ptr<spdlog::logger> config_logger;
        extern std::shared_ptr<spdlog::logger> driver_logger;

        void enable_labels();
        void disable_labels();

        void set_level(const std::string &level);
    }
}

#define MJRFD_LIB_TRACE(...) ::mjrfd::log::library_logger->trace(__VA_ARGS__)
#define MJRFD_LIB_DEBUG(...) ::mjrfd::log::library_logger->debug(__VA_ARGS__)
#define MJRFD_LIB_INFO(...)  ::mjrfd::log::library_logger->info(__VA_ARGS__)
#define MJRFD_LIB_WARN(...)  ::mjrfd::log::library_logger->warn(__VA_ARGS__)
#define MJRFD_LIB_ERROR(...) ::mjrfd::log::library_logger->error(__VA_ARGS__)
#define MJRFD_LIB_FATAL(...) ::mjrfd::log::library_logger->critical(__VA_ARGS__)

#define MJRFD_CONFIG_TRACE(...) ::mjrfd::log::config_logger->trace(__VA_ARGS__)
#define MJRFD_CONFIG_DEBUG(...) ::mjrfd::log::config_logger->debug(__VA_ARGS__)
#define MJRFD_CONFIG_INFO(...)  ::mjrfd::log::config_logger->info(__VA_ARGS__)
#define MJRFD_CONFIG_WARN(...)  ::mjrfd::log::config_logger->warn(__VA_ARGS__)
#define MJRFD_CONFIG_ERROR(...) ::mjrfd::log::config_logger->error(__VA_ARGS__)
#define MJRFD_CONFIG_FATAL(...) ::mjrfd::log::config_logger->critical(__VA_ARGS__)

#define MJRFD_TRACE(...) ::mjrfd::log::driver_logger->trace(__VA_ARGS__)
#define MJRFD_DEBUG(...) ::mjrfd::log::driver_logger->debug(__VA_ARGS__)
#define MJRFD_INFO(...)  ::mjrfd::log::driver_logger->info(__VA_ARGS__)
#define MJRFD_WARN(...)  ::mjrfd::log::driver_logger->warn(__VA_ARGS__)
#define MJRFD_ERROR(...) ::mjrfd::log::driver_logger->error(__VA_ARGS__)
#define MJRFD_FATAL(...) ::mjrfd::log::driver_logger->critical(__VA_ARGS__)
