#include "linear_particles.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SurfaceParticles(sph_body, base_material), b_n_(nullptr), width_(nullptr)
{
    //----------------------------------------------------------------------
    // Indicate particle reload variable on geometries
    // without allocated memory for the variable.
    //----------------------------------------------------------------------
    addVariableToReload<Vecd>("BinormalDirection");
    addVariableToReload<Real>("Width");
}
//=================================================================================================//
void LinearParticles::initializeBasicParticleVariables()
{
    SurfaceParticles::initializeBasicParticleVariables();
    b_n_ = getVariableDataByName<Vecd>("BinormalDirection");
    width_ = getVariableDataByName<Real>("Width");
    addVariableToWrite<Vecd>("BinormalDirection");
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
