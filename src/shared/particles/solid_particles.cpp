#include "solid_particles.h"
#include "solid_particles_variable.h"

#include "base_body.h"
#include "elastic_solid.h"
#include "inelastic_solid.h"
#include "xml_engine.h"

namespace SPH
{
//=============================================================================================//
SolidParticles::SolidParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : BaseParticles(sph_body, base_material) {}
//=================================================================================================//
void SolidParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
    registerVariable(pos0_, "InitialPosition", [&](size_t i) -> Vecd
                     { return pos_[i]; });
    registerVariable(n_, "NormalDirection");
    registerVariable(n0_, "InitialNormalDirection", [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(B_, "LinearGradientCorrectionMatrix", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
}
//=============================================================================================//
ElasticSolidParticles::ElasticSolidParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : SolidParticles(sph_body, base_material) {}
//=================================================================================================//
void ElasticSolidParticles::initializeOtherVariables()
{
    SolidParticles::initializeOtherVariables();
    /**
     *	register particle data
     */
    registerVariable(F_, "DeformationGradient", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(dF_dt_, "DeformationRate");
    /**
     *	add restart output particle data
     */
    addVariableToRestart<Matd>("DeformationGradient");
    /**
     *	add restart output particle data
     */
    addVariableToWrite<Vecd>("NormalDirection");
    addDerivedVariableToWrite<Displacement>();
    addDerivedVariableToWrite<VonMisesStress>();
    addDerivedVariableToWrite<VonMisesStrain>();
    addVariableToRestart<Matd>("DeformationGradient");
}
//=============================================================================================//
ShellParticles::ShellParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : ElasticSolidParticles(sph_body, base_material), thickness_ref_(1.0)
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
    addVariableToList<Vecd>(variables_to_reload_, "NormalDirection");
    addVariableToList<Real>(variables_to_reload_, "Thickness");
}
//=================================================================================================//
void ShellParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
    /**
     * register particle data
     */
    registerTransformationMatrix();
    registerVariable(pos0_, "InitialPosition",
                     [&](size_t i) -> Vecd
                     { return pos_[i]; });
    registerVariable(n0_, "InitialNormalDirection",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(B_, "LinearGradientCorrectionMatrix", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(F_, "DeformationGradient", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(dF_dt_, "DeformationRate");
    registerVariable(pseudo_n_, "PseudoNormal",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(dpseudo_n_dt_, "PseudoNormalChangeRate");
    registerVariable(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
    registerVariable(rotation_, "Rotation");
    registerVariable(angular_vel_, "AngularVelocity");
    registerVariable(dangular_vel_dt_, "AngularAcceleration");
    registerVariable(F_bending_, "BendingDeformationGradient");
    registerVariable(dF_bending_dt_, "BendingDeformationGradientChangeRate");
    registerVariable(global_shear_stress_, "GlobalShearStress");
    registerVariable(global_stress_, "GlobalStress");
    registerVariable(global_moment_, "GlobalMoment");
    registerVariable(mid_surface_cauchy_stress_, "MidSurfaceCauchyStress");
    /**
     * for rotation.
     */
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToRestart<Vecd>("PseudoNormal");
    addVariableToRestart<Vecd>("Rotation");
    addVariableToRestart<Vecd>("AngularVelocity");
    /**
     * add basic output particle data
     */
    addVariableToWrite<Vecd>("NormalDirection");
    addDerivedVariableToWrite<Displacement>();
    addDerivedVariableToWrite<VonMisesStress>();
    addDerivedVariableToWrite<VonMisesStrain>();
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToWrite<Vecd>("Rotation");
    addDerivedVariableToWrite<MidSurfaceVonMisesStress>();
}
//=================================================================================================//
void ShellParticles::registerTransformationMatrix()
{
    registerVariable(transformation_matrix0_, "TransformationMatrix", [&](size_t index_i) -> Matd
                     { return getTransformationMatrix(n_[index_i]); });
}
//=================================================================================================//
} // namespace SPH
