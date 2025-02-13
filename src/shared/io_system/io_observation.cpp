
#include "io_observation.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
BaseQuantityRecording::BaseQuantityRecording(SPHSystem &sph_system,
                                             const std::string &dynamics_identifier_name,
                                             const std::string &quantity_name)
    : BaseIO(sph_system), plt_engine_(),
      dynamics_identifier_name_(dynamics_identifier_name),
      quantity_name_(quantity_name),
      filefullpath_output_(io_environment_.output_folder_ + "/" +
                           dynamics_identifier_name_ + "_" + quantity_name + ".dat") {}
//=================================================================================================//
} // namespace SPH
