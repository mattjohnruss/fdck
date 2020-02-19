#pragma once

#include <memory>

#include <spdlog/spdlog.h>

namespace fdck
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

#define FDCK_LIB_TRACE(...) ::fdck::log::library_logger->trace(__VA_ARGS__)
#define FDCK_LIB_DEBUG(...) ::fdck::log::library_logger->debug(__VA_ARGS__)
#define FDCK_LIB_INFO(...)  ::fdck::log::library_logger->info(__VA_ARGS__)
#define FDCK_LIB_WARN(...)  ::fdck::log::library_logger->warn(__VA_ARGS__)
#define FDCK_LIB_ERROR(...) ::fdck::log::library_logger->error(__VA_ARGS__)
#define FDCK_LIB_FATAL(...) ::fdck::log::library_logger->critical(__VA_ARGS__)

#define FDCK_CONFIG_TRACE(...) ::fdck::log::config_logger->trace(__VA_ARGS__)
#define FDCK_CONFIG_DEBUG(...) ::fdck::log::config_logger->debug(__VA_ARGS__)
#define FDCK_CONFIG_INFO(...)  ::fdck::log::config_logger->info(__VA_ARGS__)
#define FDCK_CONFIG_WARN(...)  ::fdck::log::config_logger->warn(__VA_ARGS__)
#define FDCK_CONFIG_ERROR(...) ::fdck::log::config_logger->error(__VA_ARGS__)
#define FDCK_CONFIG_FATAL(...) ::fdck::log::config_logger->critical(__VA_ARGS__)

#define FDCK_TRACE(...) ::fdck::log::driver_logger->trace(__VA_ARGS__)
#define FDCK_DEBUG(...) ::fdck::log::driver_logger->debug(__VA_ARGS__)
#define FDCK_INFO(...)  ::fdck::log::driver_logger->info(__VA_ARGS__)
#define FDCK_WARN(...)  ::fdck::log::driver_logger->warn(__VA_ARGS__)
#define FDCK_ERROR(...) ::fdck::log::driver_logger->error(__VA_ARGS__)
#define FDCK_FATAL(...) ::fdck::log::driver_logger->critical(__VA_ARGS__)
