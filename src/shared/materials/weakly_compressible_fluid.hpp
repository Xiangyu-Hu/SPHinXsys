#ifndef WEAKLY_COMPRESSIBLE_FLUID_HPP
#define WEAKLY_COMPRESSIBLE_FLUID_HPP

#include "weakly_compressible_fluid.h"

namespace SPH
{
//=================================================================================================//
template <class FirstFluidType, typename... Args>
WeaklyCompressibleMultiPhase::WeaklyCompressibleMultiPhase(Real c0, Args &&...args)
    : c0_(c0)
{
    addPhase<FirstFluidType>(std::forward<Args>(args)...);
}
//=================================================================================================//
template <class FluidType, typename... Args>
void WeaklyCompressibleMultiPhase::addPhase(Args &&...args)
{
    if (is_phases_set_)
    {
        std::cout << "\n Error in WeaklyCompressibleMultiPhase::addPhase :"
                  << " Phases have been set, cannot add more phase ! \n ";
        exit(1);
    }
    auto &fluid = fluid_ptrs_.createPtr<FluidType>(std::forward<Args>(args)...);
    phase_list_.push_back(fluid);
    phase_name_list_.push_back(fluid->Name());
}
//=================================================================================================//
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_HPP