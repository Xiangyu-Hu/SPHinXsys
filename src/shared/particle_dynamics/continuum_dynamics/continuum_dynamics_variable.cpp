#include "continuum_dynamics_variable.h"

namespace SPH
{
namespace continuum_dynamics
{
//=============================================================================================//
VonMisesStress::VonMisesStress(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStress"),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      shear_stress_(particles_->getVariableDataByName<Matd>("ShearStress")) {}
//=============================================================================================//
void VonMisesStress::update(size_t index_i, Real dt)
{
    Matd stress_tensor = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
    derived_variable_[index_i] = getVonMisesStressFromMatrix(stress_tensor);
}
//=============================================================================================//
VonMisesStrain::VonMisesStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStrain"),
      strain_tensor_(particles_->getVariableDataByName<Matd>("StrainTensor")) {}
//=============================================================================================//
void VonMisesStrain::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = getVonMisesStressFromMatrix(strain_tensor_[index_i]);
}
//=============================================================================================//
VerticalStress::VerticalStress(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VerticalStress"),
      stress_tensor_3D_(particles_->getVariableDataByName<Mat3d>("StressTensor3D")) {}
//=============================================================================================//
void VerticalStress::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = stress_tensor_3D_[index_i](1, 1);
}
//=============================================================================================//
AccDeviatoricPlasticStrain::AccDeviatoricPlasticStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "AccDeviatoricPlasticStrain"),
      plastic_continuum_(DynamicCast<PlasticContinuum>(this, sph_body_.getBaseMaterial())),
      stress_tensor_3D_(particles_->getVariableDataByName<Mat3d>("StressTensor3D")),
      strain_tensor_3D_(particles_->getVariableDataByName<Mat3d>("StrainTensor3D")),
      E_(plastic_continuum_.getYoungsModulus()), nu_(plastic_continuum_.getPoissonRatio()) {}
//=============================================================================================//
void AccDeviatoricPlasticStrain::update(size_t index_i, Real dt)
{
    Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1.0 / 3.0) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
    Real hydrostatic_pressure = (1.0 / 3.0) * stress_tensor_3D_[index_i].trace();
    Mat3d elastic_strain_tensor_3D = deviatoric_stress / (2.0 * plastic_continuum_.getShearModulus(E_, nu_)) +
                                     hydrostatic_pressure * Mat3d::Identity() / (9.0 * plastic_continuum_.getBulkModulus(E_, nu_));
    Mat3d plastic_strain_tensor_3D = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D;
    Mat3d deviatoric_strain_tensor = plastic_strain_tensor_3D - (1.0 / (Real)Dimensions) * plastic_strain_tensor_3D.trace() * Mat3d::Identity();
    Real sum = (deviatoric_strain_tensor.cwiseProduct(deviatoric_strain_tensor)).sum();
    derived_variable_[index_i] = sqrt(sum * 2.0 / 3.0);
}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
