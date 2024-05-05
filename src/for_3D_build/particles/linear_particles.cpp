#include "linear_particles.h"

namespace SPH
{
//=============================================================================================//
LinearParticles::LinearParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SurfaceParticles(sph_body, base_material)
{
    //----------------------------------------------------------------------
    //		modify kernel function for surface particles
    //----------------------------------------------------------------------
    // sph_body.sph_adaptation_->getKernel()->reduceTwice();
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    registerVariable(b_n_, "BinormalDirection");
    registerVariable(width_, "Width");
    /**
     * add particle reload data
     */
    // addVariableNameToList<Vecd>(variables_to_reload_, "BinormalDirection");
    // addVariableNameToList<Real>(variables_to_reload_, "Width");
}
//=================================================================================================//
void LinearParticles::initializeOtherVariables()
{
    SurfaceParticles::initializeOtherVariables();

    addVariableToWrite<Vecd>("BinormalDirection");
}
//=================================================================================================//
void LinearParticles::registerTransformationMatrix()
{
    registerVariable(transformation_matrix0_, "TransformationMatrix", [&](size_t index_i) -> Matd
                     { return getTransformationMatrix(n_[index_i], b_n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
