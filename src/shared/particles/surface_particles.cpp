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
}
//=================================================================================================//
void SurfaceParticles::initializeBasicParticleVariables()
{
    BaseParticles::initializeBasicParticleVariables();
    registerTransformationMatrix();
}
//=================================================================================================//
void SurfaceParticles::registerSurfaceProperties(StdLargeVec<Vecd> &n, StdLargeVec<Real> &thickness)
{
    n_ = registerStateVariableFrom<Vecd>("NormalDirection", n);
    thickness_ = registerStateVariableFrom<Real>("Thickness", thickness);
    addVariableToReload<Vecd>("NormalDirection");
    addVariableToReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::registerSurfacePropertiesFromReload()
{
    n_ = registerStateVariableFromReload<Vecd>("NormalDirection");
    thickness_ = registerStateVariableFromReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerStateVariable<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix(n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
