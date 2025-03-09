#ifndef CUSTOM_IO_ENVIRONMENT_H
#define CUSTOM_IO_ENVIRONMENT_H

#include "io_environment.h"
#include "sph_system.h"

namespace SPH
{
class CustomIOEnvironment : public IOEnvironment
{
  public:
    // Constructor with an additional environment_number parameter
    CustomIOEnvironment(SPHSystem &sph_system, bool delete_output, int parallel_env_number, int episode_number);
};
} // namespace SPH
#endif // CUSTOM_IO_ENVIRONMENT_H
