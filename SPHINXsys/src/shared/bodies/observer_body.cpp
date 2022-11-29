#include "observer_body.h"

#include "sph_system.h"
#include "base_particles.hpp"
#include "base_material.h"
namespace SPH
{
    //=================================================================================================//
    ObserverBody::ObserverBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr)
        : SPHBody(sph_system, shape_ptr)
    {
        defineParticlesAndMaterial();
        sph_system.observation_bodies_.push_back(this);
    }
    //=================================================================================================//
    ObserverBody::ObserverBody(SPHSystem &sph_system, const std::string &name)
        : ObserverBody(sph_system, makeShared<DefaultShape>(name)) {}
    //=================================================================================================//
}
