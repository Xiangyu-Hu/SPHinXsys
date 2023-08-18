/**
 * @file 	complex_solid.cpp
 * @brief 	These are classes for define complex solid materials.
 * @author	Yaru Ren and Xiangyu Hu
 */

#include "complex_solid.h"

#include "base_local_dynamics.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
CompositeSolid::CompositeSolid(Real rho0)
    : ElasticSolid(rho0)
{
    material_type_name_ = "CompositeSolid";
}
//=================================================================================================//
Matd CompositeSolid::StressPK2(Matd &deformation, size_t index_i)
{
    return composite_materials_[material_id_[index_i]]->StressPK2(deformation, index_i);
}
//=================================================================================================//
Matd CompositeSolid::StressPK1(Matd &deformation, size_t index_i)
{
    return composite_materials_[material_id_[index_i]]->StressPK1(deformation, index_i);
}
//=================================================================================================//
Matd CompositeSolid::StressCauchy(Matd &almansi_strain, Matd &F, size_t index_i)
{
    return Matd::Identity();
}
//=================================================================================================//
Real CompositeSolid::CompositeDensity(size_t index_i)
{
    return composite_materials_[material_id_[index_i]]->ReferenceDensity();
}
//=================================================================================================//
void CompositeSolid::initializeLocalParameters(BaseParticles *base_particles)
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
} // namespace SPH
