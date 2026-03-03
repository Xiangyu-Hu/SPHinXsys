
#include "io_environment.h"

#include "sph_system.h"

namespace SPH
{
//=============================================================================================//
IOEnvironment::IOEnvironment(SPHSystem &sph_system)
    : sph_system_(sph_system),
      input_folder_("./input"), output_folder_("./output"),
      restart_folder_("./restart"), reload_folder_("./reload")
{
    if (!fs::exists(input_folder_))
    {
        fs::create_directory(input_folder_);
    }

    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }
    else
    {
        std::string backup_name = output_folder_ + "_backup";
        fs::path output_dir = output_folder_;
        fs::path backup_path = output_dir.parent_path() / backup_name;
        if (fs::exists(backup_path))
        {
            fs::remove_all(backup_path);
        }
        fs::rename(output_dir, backup_path);
        std::cout << "Moved existing output to: " << backup_path << std::endl;
        fs::create_directory(output_folder_);
    }

    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }

    if (!fs::exists(reload_folder_))
    {
        fs::create_directory(reload_folder_);
    }

    sph_system.io_environment_ = this;
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
void IOEnvironment::resetOutputFolder(const std::string &new_name)
{
    output_folder_ = new_name;
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
//=============================================================================================//
ParameterizationIO *IOEnvironment::defineParameterizationIO()
{
    return parameterization_io_keeper_.createPtr<ParameterizationIO>(input_folder_);
}
//=================================================================================================//
} // namespace SPH
