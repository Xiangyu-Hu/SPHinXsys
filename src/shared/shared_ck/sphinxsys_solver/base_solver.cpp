#include "base_solver.h"

#include "sph_system.hpp"

namespace SPH
{
//=================================================================================================//
SPHSolver::SPHSolver(SPHSystem &sph_system)
    : physical_time_(sph_system.getSystemVariableByName<Real>("PhysicalTime")) {}
//=================================================================================================//
} // namespace SPH
