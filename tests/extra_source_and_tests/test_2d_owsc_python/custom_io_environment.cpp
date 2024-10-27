#include "custom_io_environment.h"
#include "sph_system.h"
namespace fs = std::filesystem;

namespace SPH
{
//=================================================================================================//
CustomIOEnvironment::CustomIOEnvironment(SPHSystem &sph_system, bool delete_output, int parallel_env_number, int episode_number)
    : IOEnvironment(sph_system, delete_output) 
{
    // Append environment_number to the output_folder_
    output_folder_ += "_env_" + std::to_string(parallel_env_number) + "_episode_" + std::to_string(episode_number);

    // Check and create the output folder with the modified path
    if (!fs::exists(output_folder_)) {
        fs::create_directory(output_folder_);
    }

    // Handle deletion of contents in the output folder if required
    if (delete_output) {
        fs::remove_all(output_folder_);
        fs::create_directory(output_folder_);
    }
}
//=================================================================================================//
} // namespace SPH