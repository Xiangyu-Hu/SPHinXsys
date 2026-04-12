
#include "io_environment.h"

#include "parameterization.h"
#include "spdlog/sinks/rotating_file_sink.h"

namespace SPH
{
//=============================================================================================//
IOEnvironment::~IOEnvironment() = default;
//=============================================================================================//
IOEnvironment::IOEnvironment()
    : input_folder_("./input"), output_folder_("./output"),
      restart_folder_("./restart"), reload_folder_("./reload")
{
    setOutputFolder(output_folder_);

    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }

    if (!fs::exists(reload_folder_))
    {
        fs::create_directory(reload_folder_);
    }
}
//=================================================================================================//
void IOEnvironment::setOutputFolder(const std::string &folder_name)
{
    if (!fs::exists(folder_name))
    {
        fs::create_directory(folder_name);
    }
    else
    {
        std::string backup_folder = folder_name + "_backup";
        fs::path output_dir = output_folder_;
        fs::path backup_dir = backup_folder;
        if (fs::exists(backup_dir))
        {
            fs::remove_all(backup_dir);
        }
        fs::rename(output_dir, backup_dir);
        std::cout << "Moved existing output to: " << backup_dir << std::endl;
        fs::create_directory(folder_name);
    }
}
//=================================================================================================//
void IOEnvironment::resetForRestart()
{
    fs::remove_all(restart_folder_);
    fs::create_directory(restart_folder_);
}
//=================================================================================================//
void IOEnvironment::reinitializeReloadFolder()
{
    fs::remove_all(reload_folder_);
    fs::create_directory(reload_folder_);
}
//=================================================================================================//
void IOEnvironment::appendOutputFolder(const std::string &append_name)
{
    output_folder_ += "_" + append_name;
    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }
    else
    {
        fs::remove_all(output_folder_);
        fs::create_directory(output_folder_);
    }
}
//=================================================================================================//
void IOEnvironment::resetOutputFolder(const std::string &new_name, bool keep_existing)
{
    if (fs::exists(output_folder_) && !keep_existing)
    {
        fs::remove_all(output_folder_);
    }

    output_folder_ = new_name;
    setOutputFolder(output_folder_);
}
//=================================================================================================//
void IOEnvironment::resetRestartFolder(const std::string &new_name, bool keep_existing)
{
    if (fs::exists(restart_folder_) && !keep_existing)
    {
        fs::remove_all(restart_folder_);
    }

    restart_folder_ = new_name;
    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }
    else
    {
        fs::remove_all(restart_folder_);
        fs::create_directory(restart_folder_);
    }
}
//=================================================================================================//
void IOEnvironment::resetReloadFolder(const std::string &new_name, bool keep_existing)
{
    if (fs::exists(reload_folder_) && !keep_existing)
    {
        fs::remove_all(reload_folder_);
    }

    reload_folder_ = new_name;
    if (!fs::exists(reload_folder_))
    {
        fs::create_directory(reload_folder_);
    }
    else
    {
        fs::remove_all(reload_folder_);
        fs::create_directory(reload_folder_);
    }
}
//=============================================================================================//
ParameterizationIO *IOEnvironment::defineParameterizationIO()
{
    return parameterization_io_keeper_.createPtr<ParameterizationIO>(input_folder_);
}
//=============================================================================================//
namespace IO
{
//=============================================================================================//
SharedPtr<IOEnvironment> io_environment; // Global pointer to the IO environment
//=============================================================================================//
IOEnvironment &initEnvironment()
{
    if (!io_environment)
    {
        io_environment = makeShared<IOEnvironment>();
    }
    return *io_environment.get();
}
//=============================================================================================//
IOEnvironment &getEnvironment()
{
    if (!io_environment)
    {
        throw std::runtime_error("IOEnvironment not initialized. Call IO::initEnvironment() first.");
    }
    return *io_environment.get();
}
//=================================================================================================//
std::shared_ptr<spdlog::logger> logger; // Global logger pointer
//=================================================================================================//
std::shared_ptr<spdlog::logger> initLogger()
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
std::shared_ptr<spdlog::logger> getLogger()
{
    if (!logger)
    {
        throw std::runtime_error("Logger not initialized. Call IO::initLogger() first.");
    }
    return logger;
}
//=================================================================================================//
} // namespace IO
//=================================================================================================//
} // namespace SPH
