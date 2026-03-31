#ifndef BASE_BODY_HPP
#define BASE_BODY_HPP

#include "base_body.h"

#include "adaptation.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
template <class AdaptationType, typename... Args>
void SPHBody::defineAdaptation(Args &&...args)
{
    sph_adaptation_ =
        sph_adaptation_keeper_.createPtr<AdaptationType>(
            sph_adaptation_->GlobalResolution(), std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename... Args>
LevelSetShape &SPHBody::defineComponentLevelSetShape(const std::string &shape_name, Args &&...args)
{
    ComplexShape *complex_shape = DynamicCast<ComplexShape>(this, initial_shape_);
    return complex_shape->defineLevelSetShape(*this, shape_name, std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename... Args>
LevelSetShape &SPHBody::defineBodyLevelSetShape(Args &&...args)
{
    initial_shape_ = shape_keeper_.resetPtr<LevelSetShape>(
        *this, *initial_shape_, std::forward<Args>(args)...);
    return *static_cast<LevelSetShape *>(initial_shape_);
}
//=================================================================================================//
template <typename ExecutionPolicy, typename... Args>
LevelSetShape &SPHBody::defineBodyLevelSetShape(const ExecutionPolicy &ex_policy, Args &&...args)
{
    initial_shape_ = shape_keeper_.resetPtr<LevelSetShape>(
        ex_policy, sph_system_, *sph_adaptation_, *initial_shape_, std::forward<Args>(args)...);
    return *static_cast<LevelSetShape *>(initial_shape_);
}
//=================================================================================================//
template <class MaterialType, typename... Args>
MaterialType &SPHBody::defineMaterial(Args &&...args)
{
    base_material_ = base_material_keeper_.createPtr<MaterialType>(
        std::forward<Args>(args)...);
    return *static_cast<MaterialType *>(base_material_);
}
//=================================================================================================//
template <class BaseModel, typename... AuxiliaryModels, typename... Args>
Closure<BaseModel, AuxiliaryModels...> &SPHBody::defineClosure(Args &&...args)
{
    base_material_ = base_material_keeper_.createPtr<Closure<BaseModel, AuxiliaryModels...>>(
        std::forward<Args>(args)...);
    return *static_cast<Closure<BaseModel, AuxiliaryModels...> *>(base_material_);
}
//=================================================================================================//
template <class ParticleType, class... Parameters, typename... Args>
ParticleType &SPHBody::generateParticles(Args &&...args)
{
    ParticleType *particles = base_particles_keeper_.createPtr<ParticleType>(*this, base_material_);
    ParticleGenerator<ParticleType, Parameters...> particle_generator(*this, *particles, std::forward<Args>(args)...);
    particle_generator.generateParticlesWithGeometricVariables();
    particles->initializeBasicDiscreteVariables();
    sph_adaptation_->initializeAdaptationVariables(*particles);
    base_material_->setLocalParameters(sph_system_, particles);
    return *particles;
}
//=================================================================================================//
template <class ParticleType, typename... Parameters, class ReserveType, typename... Args>
ParticleType &SPHBody::generateParticlesWithReserve(ReserveType &particle_reserve, Args &&...args)
{
    return generateParticles<ParticleType, ReserveType, Parameters...>(
        particle_reserve, std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename... Args>
RealBody::RealBody(Args &&...args)
    : SPHBody(std::forward<Args>(args)...), cell_linked_list_created_(false)
{
    addRealBodyToSPHSystem();
}
//=================================================================================================//
template <class BodyPartType, typename... Args>
BodyPartType &RealBody::addBodyPart(Args &&...args)
{
    return *body_parts_keeper_.createPtr<BodyPartType>(*this, std::forward<Args>(args)...);
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_BODY_HPP
