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
    //		register geometric data only
    //----------------------------------------------------------------------
    n_ = registerSharedVariable<Vecd>("NormalDirection");
    thickness_ = registerSharedVariable<Real>("Thickness");
    /**
     * add particle reload data
     */
    addVariableToReload<Vecd>("NormalDirection");
    addVariableToReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
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
