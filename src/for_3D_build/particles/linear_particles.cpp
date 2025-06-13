#include "linear_particles.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SurfaceParticles(sph_body, base_material), b_n_(nullptr), width_(nullptr) {}
//=================================================================================================//
void LinearParticles::registerLineProperties(StdVec<Vecd> &b_n, StdVec<Real> &width)
{
    b_n_ = registerStateVariableFrom<Vecd>("BinormalDirection", b_n);
    width_ = registerStateVariableFrom<Real>("Width", width);
    addEvolvingVariable<Vecd>("BinormalDirection");
    addEvolvingVariable<Real>("Width");
    addVariableToWrite<Vecd>("BinormalDirection");
}
//=================================================================================================//
void LinearParticles::registerLinePropertiesFromReload()
{
    b_n_ = registerStateVariableFromReload<Vecd>("NormalDirection");
    width_ = registerStateVariableFromReload<Real>("Thickness");
}
//=================================================================================================//
void LinearParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerStateVariable<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix(n_[index_i], b_n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
