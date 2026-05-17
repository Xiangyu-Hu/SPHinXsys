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
MaterialType &SPHBody::defineMatterMaterial(Args &&...args)
{
    MaterialType *material = matter_material_keeper_.createPtr<MaterialType>(std::forward<Args>(args)...);
    all_material_properties_[0] = material; // set the first material as the matter material for the body
    return *material;
}
//=================================================================================================//
template <class MaterialType, typename... Args>
MaterialType &SPHBody::addMaterialProperty(Args &&...args)
{
    MaterialType *material = material_properties_keeper_.createPtr<MaterialType>(std::forward<Args>(args)...);
    all_material_properties_.push_back(material);
    return *material;
}
//=================================================================================================//
template <typename MaterialType>
StdVec<MaterialType *> SPHBody::collectMaterialProperties()
{
    StdVec<MaterialType *> materials;
    for (auto *material : all_material_properties_)
    {
        if (auto *cast_material = dynamic_cast<MaterialType *>(material))
        {
            materials.push_back(cast_material);
        }
    }
    return materials;
}
//=================================================================================================//
template <class MaterialType>
MaterialType &SPHBody::getMaterialProperty(const std::string &name)
{
    StdVec<MaterialType *> materials = collectMaterialProperties<MaterialType>();
    if (materials.empty())
    {
        throw std::runtime_error(
            std::string(type_name<MaterialType>()) + " not found in body: " + body_name_);
    }
    else if (materials.size() == 1)
    {
        return *materials[0];
    }
    else
    {
        for (auto *material : materials)
        {
            if (material->Name() == name)
            {
                return *material;
            }
        }
        throw std::runtime_error(std::string(type_name<MaterialType>()) + ": " + name +
                                 " not found in body: " + body_name_);
    }
}
//=================================================================================================//
template <class ParticleType, class... Parameters, typename... Args>
ParticleType &SPHBody::generateParticles(Args &&...args)
{
    ParticleType *particles = base_particles_keeper_.createPtr<ParticleType>(*this);
    ParticleGenerator<ParticleType, Parameters...> particle_generator(*this, *particles, std::forward<Args>(args)...);
    particle_generator.generateParticlesWithGeometricVariables();
    particles->initializeBasicDiscreteVariables();
    sph_adaptation_->initializeAdaptationVariables(*particles);
    for (auto *material : all_material_properties_)
    {
        material->setLocalParameters(sph_system_, particles);
    }
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
