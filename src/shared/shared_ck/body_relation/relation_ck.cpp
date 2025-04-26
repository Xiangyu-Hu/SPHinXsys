#include "relation_ck.hpp"

namespace SPH
{
//=================================================================================================//
DiscreteVariable<Vecd> *Relation<Base>::assignConfigPosition(
    BaseParticles &particles, ConfigType config_type)
{
    if (config_type == ConfigType::Eulerian)
    {
        return particles.getVariableByName<Vecd>("Position");
    }
    else
    {
        return particles.registerStateVariableOnlyFrom<Vecd>(
            "InitialPosition", "Position");
    }
}
//=================================================================================================//
DiscreteVariable<Vecd> *Relation<Base>::getTargetPosition(UnsignedInt target_index)
{
    return dv_target_pos_[target_index];
}
//=================================================================================================//
DiscreteVariable<UnsignedInt> *Relation<Base>::getNeighborIndex(UnsignedInt target_index)
{
    return dv_target_neighbor_index_[target_index];
}
//=================================================================================================//
DiscreteVariable<UnsignedInt> *Relation<Base>::getParticleOffset(UnsignedInt target_index)
{
    return dv_target_particle_offset_[target_index];
}
//=================================================================================================//
void Relation<Base>::registerComputingKernel(
    execution::Implementation<Base> *implementation, UnsignedInt target_index)
{
    registered_computing_kernels_[target_index].push_back(implementation);
}
//=================================================================================================//
void Relation<Base>::resetComputingKernelUpdated(UnsignedInt target_index)
{
    auto &computing_kernels = registered_computing_kernels_[target_index];
    for (size_t k = 0; k != computing_kernels.size(); ++k)
    {
        computing_kernels[k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
