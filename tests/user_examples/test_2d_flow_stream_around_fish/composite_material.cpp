#include "composite_material.h"
#include "base_particles.hpp"

#include <numeric>

namespace SPH
{
//=================================================================================================//
void CompositeMaterial::initializeLocalParameters(BaseParticles *base_particles)
{
    ElasticSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(material_id_, "MaterialID");

    for (size_t i = 0; i < composite_materials_.size(); ++i)
    {
        composite_materials_[i]->initializeLocalParameters(base_particles);
    }
}
//=================================================================================================//
Matd ActiveModelSolid::StressPK2(Matd &F, size_t particle_index_i)
{
    return lambda0_ * F.trace() * Matd::Identity() + 2.0 * G0_ * F;
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//