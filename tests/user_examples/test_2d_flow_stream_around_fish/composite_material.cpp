#include "composite_material.h"
#include "base_local_dynamics.h"
#include "base_particles.hpp"

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
MaterialIdInitialization::MaterialIdInitialization(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      material_id_(*particles_->getVariableByName<int>("MaterialID")),
      pos0_(*particles_->getVariableByName<Vecd>("InitialPosition")){};
//=================================================================================================//
Matd ActiveModelSolid::StressPK2(Matd &F, size_t particle_index_i)
{
    return lambda0_ * F.trace() * Matd::Identity() + 2.0 * G0_ * F;
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//