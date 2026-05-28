#include "surface_particles.h"

#include "adaptation.h"
#include "base_body.h"
#include "vector_functions.h"

namespace SPH
{
//=============================================================================================//
SurfaceParticles::SurfaceParticles(SPHBody &sph_body)
    : BaseParticles(sph_body), n_(nullptr), thickness_(nullptr),
      transformation_matrix0_(nullptr)
{
    //----------------------------------------------------------------------
    //		modify kernel function for surface particles
    //----------------------------------------------------------------------
    sph_body.getSPHAdaptation().getKernel()->reduceOnce();
}
//=================================================================================================//
void SurfaceParticles::initializeBasicDiscreteVariables()
{
    BaseParticles::initializeBasicDiscreteVariables();
    registerTransformationMatrix();
}
//=================================================================================================//
void SurfaceParticles::registerSurfaceProperties(StdVec<Vecd> &n, StdVec<Real> &thickness)
{
    n_ = registerStateVariableDataFrom<Vecd>("NormalDirection", n);
    thickness_ = registerStateVariableDataFrom<Real>("Thickness", thickness);
    addVariableToReload<Vecd>("NormalDirection");
    addVariableToReload<Real>("Thickness");
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
