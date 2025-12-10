#include "io_log.h"

#include "io_environment.h"

#include "spdlog/sinks/rotating_file_sink.h"

namespace SPH
{
namespace Log
{
//=================================================================================================//
std::shared_ptr<spdlog::logger> logger; // Global logger pointer
//=================================================================================================//
std::shared_ptr<spdlog::logger> init()
{
    if (!logger) // Create a logger with a file sink
    {
        // Create a file rotating logger with 5 MB size max and 3 rotated files
        auto max_size = 1048576 * 5;
        auto max_files = 3;
        logger = spdlog::rotating_logger_mt("sphinxsys_logger", "sphinxsys.log", max_size, max_files);
        logger->flush_on(spdlog::level::info);                          // Flush logs on info level
        spdlog::set_default_logger(logger);                             // Set the default logger to the created logger
        spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%l] [thread %t] %v"); // Set the log format
        spdlog::flush_every(std::chrono::seconds(1));                   // Flush logs every second
    }
    return logger;
}
//=================================================================================================//
std::shared_ptr<spdlog::logger> get()
{
    if (!logger)
    {
        throw std::runtime_error("Logger not initialized. Call Log::init() first.");
    }
    return logger;
}
//=================================================================================================//
} // namespace Log
} // namespace SPH
