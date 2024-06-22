#include "linear_particles.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SurfaceParticles(sph_body, base_material), b_n_(nullptr), width_(nullptr)
{
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    b_n_ = registerSharedVariable<Vecd>("BinormalDirection");
    width_ = registerSharedVariable<Real>("Width");
    /**
     * add particle reload data
     */
    addVariableToReload<Vecd>("BinormalDirection");
    addVariableToReload<Real>("Width");
}
//=================================================================================================//
void LinearParticles::initializeBasicParticleVariables()
{
    SurfaceParticles::initializeBasicParticleVariables();

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
