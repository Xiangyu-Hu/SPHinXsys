#include "linear_particles.h"

#include "vector_functions.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body)
    : SurfaceParticles(sph_body), b_n_(nullptr), width_(nullptr) {}
//=================================================================================================//
void LinearParticles::registerLineProperties(StdVec<Vecd> &b_n, StdVec<Real> &width)
{
    b_n_ = registerStateVariableDataFrom<Vecd>("BinormalDirection", b_n);
    width_ = registerStateVariableDataFrom<Real>("Width", width);
    addEvolvingVariable<Vecd>("BinormalDirection");
    addEvolvingVariable<Real>("Width");
    addVariableToWrite<Vecd>("BinormalDirection");
}
//=================================================================================================//
void LinearParticles::registerLinePropertiesFromReload()
{
    b_n_ = registerStateVariableDataFromReload<Vecd>("BinormalDirection");
    width_ = registerStateVariableDataFromReload<Real>("Width");
}
//=================================================================================================//
void LinearParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerStateVariableData<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix(n_[index_i], b_n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
