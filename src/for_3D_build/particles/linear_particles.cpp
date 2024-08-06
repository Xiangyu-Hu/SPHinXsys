#include "linear_particles.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SurfaceParticles(sph_body, base_material), b_n_(nullptr), width_(nullptr) {}
//=================================================================================================//
void LinearParticles::registerLineProperties(StdLargeVec<Vecd> &b_n, StdLargeVec<Real> &width)
{
    b_n_ = registerSharedVariableFrom<Vecd>("BinormalDirection", b_n);
    width_ = registerSharedVariableFrom<Real>("Width", width);
    addVariableToReload<Vecd>("BinormalDirection");
    addVariableToReload<Real>("Width");
    addVariableToWrite<Vecd>("BinormalDirection");
}
//=================================================================================================//
void LinearParticles::registerLinePropertiesFromReload()
{
    b_n_ = registerSharedVariableFromReload<Vecd>("NormalDirection");
    width_ = registerSharedVariableFromReload<Real>("Thickness");
}
//=================================================================================================//
void LinearParticles::registerTransformationMatrix()
{
    transformation_matrix0_ = registerSharedVariable<Matd>(
        "TransformationMatrix", [&](size_t index_i) -> Matd
        { return getTransformationMatrix((*n_)[index_i], (*b_n_)[index_i]); });
}
//=================================================================================================//
} // namespace SPH
