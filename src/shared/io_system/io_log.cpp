#include "io_log.h"

#include "io_environment.h"

#include <spdlog/sinks/basic_file_sink.h>

namespace SPH
{
namespace Log
{
//=================================================================================================//
std::shared_ptr<spdlog::logger> logger; // Global logger pointer
//=================================================================================================//
void init(IOEnvironment &io_environment)
{
    if (!logger) // Create a logger with a file sink
    {
        logger = spdlog::basic_logger_mt("sphinxsys_logger", io_environment.outputFolder() + "/sphinxsys.log");
        logger->flush_on(spdlog::level::info);                          // Flush logs on info level
        spdlog::set_default_logger(logger);                             // Set the default logger to the created logger
        spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%l] [thread %t] %v"); // Set the log format
        spdlog::flush_every(std::chrono::seconds(1));                   // Flush logs every second
    }
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
