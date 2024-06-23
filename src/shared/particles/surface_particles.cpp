#include "surface_particles.h"

#include "base_body.h"

namespace SPH
{
//=============================================================================================//
SurfaceParticles::SurfaceParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : BaseParticles(sph_body, base_material), n_(nullptr), thickness_(nullptr),
      transformation_matrix0_(nullptr)
{
    //----------------------------------------------------------------------
    //		modify kernel function for surface particles
    //----------------------------------------------------------------------
    sph_body.sph_adaptation_->getKernel()->reduceOnce();
    //----------------------------------------------------------------------
    // add geometries variable which will initialized particle generator
    //----------------------------------------------------------------------
    addSharedVariable<Vecd>("NormalDirection");
    addSharedVariable<Real>("Thickness");
    /**
     * add particle reload variables
     */
    addVariableToReload<Vecd>("NormalDirection");
    addVariableToReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::initializeBasicParticleVariables()
{
    BaseParticles::initializeBasicParticleVariables();
    n_ = getVariableDataByName<Vecd>("NormalDirection");
    thickness_ = getVariableDataByName<Real>("Thickness");
    registerTransformationMatrix();
}
//=================================================================================================//
void SurfaceParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerSharedVariable<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix((*n_)[index_i]); });
}
//=================================================================================================//
} // namespace SPH
