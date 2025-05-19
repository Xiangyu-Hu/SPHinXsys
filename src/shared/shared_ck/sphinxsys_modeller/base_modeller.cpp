#include "base_modeller.h"

#include "sph_system.hpp"

namespace SPH
{
//=================================================================================================//
SPHModeller::SPHModeller(SPHSystem &sph_system)
    : physical_time_(sph_system.getSystemVariableByName<Real>("PhysicalTime")) {}
//=================================================================================================//
} // namespace SPH
