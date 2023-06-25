#include "continuum_particles.h"
#include "general_continuum.h"
namespace SPH
{
    ContinuumParticles::
        ContinuumParticles(SPHBody &sph_body, GeneralContinuum *granular_material)
        : BaseParticles(sph_body, granular_material), granular_material_(*granular_material) {}
    //=================================================================================================//
    void ContinuumParticles::initializeOtherVariables()
    {
        BaseParticles::initializeOtherVariables();

        registerVariable(acc_shear_, "AccelerationByShear");
        registerVariable(stress_tensor_, "StressTensor");
        registerVariable(stress_tensor_rate_, "StressTensorRate");
        registerVariable(shear_stress_, "ShearStress");
        registerVariable(shear_stress_rate_, "ShearStressRate");
        registerVariable(von_mises_stress_, "VonMisesStress");
        registerVariable(velocity_gradient_, "VelocityGradient");
        registerVariable(strain_tensor_, "StrainTensor");
        registerVariable(strain_tensor_rate_, "StrainTensorRate");
        //----------------------------------------------------------------------
        //		register sortable particle data
        //----------------------------------------------------------------------

        registerSortableVariable<Vecd>("AccelerationByShear");
        registerSortableVariable<Matd>("StressTensor");
        registerSortableVariable<Matd>("StressTensorRate");
        registerSortableVariable<Matd>("ShearStress");
        registerSortableVariable<Matd>("ShearStressRate");
        registerSortableVariable<Real>("VonMisesStress");
        registerSortableVariable<Matd>("VelocityGradient");
        registerSortableVariable<Matd>("StrainTensor");
        registerSortableVariable<Matd>("StrainTensorRate");
        //----------------------------------------------------------------------
        registerVariable(pos0_, "InitialPosition", [&](size_t i) -> Vecd
                         { return pos_[i]; });
        registerVariable(n_, "NormalDirection");
        registerVariable(n0_, "InitialNormalDirection",
                         [&](size_t i) -> Vecd
                         { return n_[i]; });
    }

} // namespace SPH