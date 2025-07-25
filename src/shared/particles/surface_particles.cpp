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
    sph_body.getSPHAdaptation().getKernel()->reduceOnce();
}
//=================================================================================================//
void SurfaceParticles::initializeBasicParticleVariables()
{
    BaseParticles::initializeBasicParticleVariables();
    registerTransformationMatrix();
}
//=================================================================================================//
void SurfaceParticles::registerSurfaceProperties(StdVec<Vecd> &n, StdVec<Real> &thickness)
{
    n_ = registerStateVariableDataFrom<Vecd>("NormalDirection", n);
    thickness_ = registerStateVariableDataFrom<Real>("Thickness", thickness);
    addEvolvingVariable<Vecd>("NormalDirection");
    addEvolvingVariable<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::registerSurfacePropertiesFromReload()
{
    n_ = registerStateVariableDataFromReload<Vecd>("NormalDirection");
    thickness_ = registerStateVariableDataFromReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerStateVariableData<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix(n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
