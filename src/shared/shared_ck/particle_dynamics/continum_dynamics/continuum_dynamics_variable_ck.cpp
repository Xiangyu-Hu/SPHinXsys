#include "continuum_dynamics_variable_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=============================================================================================//
VonMisesStressCK::VonMisesStressCK(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStress"),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_derived_variable_(particles_->registerStateVariable<Real>("VonMisesStress")),
      dv_shear_stress_(particles_->getVariableByName<Matd>("ShearStress")) {}
//=============================================================================================//
VerticalStressCK::VerticalStressCK(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VerticalStress"),
      dv_stress_tensor_3D_(particles_-> registerStateVariable<Mat3d>("StressTensor3D")),
      dv_derived_variable_(particles_-> registerStateVariable<Real>("VerticalStress")) {}
//=============================================================================================//
AccDeviatoricPlasticStrainCK::AccDeviatoricPlasticStrainCK(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "AccDeviatoricPlasticStrain"),
      plastic_continuum_(DynamicCast<PlasticContinuum>(this, sph_body_->getBaseMaterial())),
      dv_stress_tensor_3D_(particles_-> registerStateVariable<Mat3d>("StressTensor3D")),
      dv_strain_tensor_3D_(particles_-> registerStateVariable<Mat3d>("StrainTensor3D")),
      dv_derived_variable_(particles_-> registerStateVariable<Real>("AccDeviatoricPlasticStrain")),
      E_(plastic_continuum_.getYoungsModulus()), nu_(plastic_continuum_.getPoissonRatio()) {}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
