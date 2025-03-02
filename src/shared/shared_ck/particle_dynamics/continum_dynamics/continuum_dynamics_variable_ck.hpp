#ifndef CONTINUUM_DYNAMICS_VARIABLE_CK_HPP
#define CONTINUUM_DYNAMICS_VARIABLE_CK_HPP


#include "continuum_dynamics_variable_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
    //=============================================================================================//
    VerticalStressCK::VerticalStressCK(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VerticalStress"),
    dv_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("StressTensor3D")),
    dv_derived_variable_(this->particles_->template registerStateVariableOnly<Real>("VerticalStress")){}
    //=================================================================================================//
    template <class ExecutionPolicy, class EncloserType>
    VerticalStressCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser):
    stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
    derived_variable_(encloser.dv_derived_variable_->DelegatedData(ex_policy)){}
    //=================================================================================================//   
    void VerticalStressCK::UpdateKernel::update(size_t index_i, Real dt)
    {
        derived_variable_[index_i] = stress_tensor_3D_[index_i](1, 1);
    }
    //=============================================================================================//
    AccDeviatoricPlasticStrainCK::AccDeviatoricPlasticStrainCK(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "AccDeviatoricPlasticStrain"),
    plastic_continuum_(DynamicCast<PlasticContinuum>(this, this->sph_body_.getBaseMaterial())),
    dv_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("StressTensor3D")),
    dv_strain_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("StrainTensor3D")),
    dv_derived_variable_(this->particles_->template registerStateVariableOnly<Real>("AccDeviatoricPlasticStrain")),
    E_(plastic_continuum_.getYoungsModulus()),nu_(plastic_continuum_.getPoissonRatio()){}
    //=================================================================================================//
    template <class ExecutionPolicy, class EncloserType>
    AccDeviatoricPlasticStrainCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser):
    plastic_kernel_(encloser.plastic_continuum_),
    stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
    strain_tensor_3D_(encloser.dv_strain_tensor_3D_->DelegatedData(ex_policy)),
    derived_variable_(encloser.dv_derived_variable_->DelegatedData(ex_policy)),
    E_(encloser.E_),nu_(encloser.nu_){}
    //=================================================================================================//   
    void AccDeviatoricPlasticStrainCK::UpdateKernel::update(size_t index_i, Real dt)
    {
        Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1.0 / 3.0) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
        Real hydrostatic_pressure = (1.0 / 3.0) * stress_tensor_3D_[index_i].trace();
        Mat3d elastic_strain_tensor_3D = deviatoric_stress / (2.0 * plastic_kernel_.getShearModulus(E_, nu_)) +
                                        hydrostatic_pressure * Mat3d::Identity() / (9.0 * plastic_kernel_.getBulkModulus(E_, nu_));
        Mat3d plastic_strain_tensor_3D = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D;
        Mat3d deviatoric_strain_tensor = plastic_strain_tensor_3D - (1.0 / (Real)Dimensions) * plastic_strain_tensor_3D.trace() * Mat3d::Identity();
        Real sum = (deviatoric_strain_tensor.cwiseProduct(deviatoric_strain_tensor)).sum();
        derived_variable_[index_i] = sqrt(sum * 2.0 / 3.0);
    }
} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_VARIABLE_CK_HPP