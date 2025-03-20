
#include "io_observation.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
BaseQuantityRecording::BaseQuantityRecording(SPHSystem &sph_system,
                                             const std::string &dynamics_identifier_name)
    : BaseIO(sph_system), plt_engine_(), 
      quantity_name_("NeedAQuantityName"),
      dynamics_identifier_name_(dynamics_identifier_name),
      filefullpath_output_(io_environment_.output_folder_ + "/" +
                           dynamics_identifier_name_ + "_" + "NeedAQuantityName" + ".dat") {}
//=================================================================================================//
} // namespace SPH
