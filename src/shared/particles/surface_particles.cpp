#include "surface_particles.h"

#include "base_body.h"

namespace SPH
{
//=============================================================================================//
SurfaceParticles::SurfaceParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : BaseParticles(sph_body, base_material)
{
    //----------------------------------------------------------------------
    //		modify kernel function for surface particles
    //----------------------------------------------------------------------
    sph_body.sph_adaptation_->getKernel()->reduceOnce();
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    registerVariable(n_, "NormalDirection");
    registerVariable(thickness_, "Thickness");
    /**
     * add particle reload data
     */
    addVariableToReload<Vecd>("NormalDirection");
    addVariableToReload<Real>("Thickness");
}
//=================================================================================================//
void SurfaceParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
    registerTransformationMatrix();
}
//=================================================================================================//
void SurfaceParticles::registerTransformationMatrix()
{
    registerVariable(transformation_matrix0_, "TransformationMatrix", [&](size_t index_i) -> Matd
                     { return getTransformationMatrix(n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
