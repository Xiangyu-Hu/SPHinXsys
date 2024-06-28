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
    n_ = getVariableDataByName<Vecd>("NormalDirection");
    thickness_ = getVariableDataByName<Real>("Thickness");
    registerTransformationMatrix();
}
//=================================================================================================//
StdLargeVec<Real> *SurfaceParticles::registerParticleMass(StdLargeVec<Real> *rho)
{
    return registerSharedVariable<Real>(
        "Mass", [&](size_t i) -> Real { return (*rho)[i] * (*Vol_)[i] * (*thickness_)[i]; });
}
//=================================================================================================//
void SurfaceParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerSharedVariable<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd { return getTransformationMatrix((*n_)[index_i]); });
}
//=================================================================================================//
} // namespace SPH
