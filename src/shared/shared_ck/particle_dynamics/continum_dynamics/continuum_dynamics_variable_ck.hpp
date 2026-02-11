#ifndef CONTINUUM_DYNAMICS_VARIABLE_CK_HPP
#define CONTINUUM_DYNAMICS_VARIABLE_CK_HPP

#include "continuum_dynamics_variable_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=============================================================================================//
template <class ExecutionPolicy, class EncloserType>
VonMisesStressCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : p_(encloser.dv_p_->DelegatedData(ex_policy)),
      derived_variable_(encloser.dv_derived_variable_->DelegatedData(ex_policy)),
      shear_stress_(encloser.dv_shear_stress_->DelegatedData(ex_policy)) {}
//=============================================================================================//
inline void VonMisesStressCK::UpdateKernel::update(size_t index_i, Real dt)
{
    Matd stress_tensor = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
    derived_variable_[index_i] = getVonMisesStressFromMatrix(stress_tensor);
}
//=============================================================================================//
template <class ExecutionPolicy, class EncloserType>
VerticalStressCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      derived_variable_(encloser.dv_derived_variable_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void VerticalStressCK::UpdateKernel::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = stress_tensor_3D_[index_i](1, 1);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
AccDeviatoricPlasticStrainCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : constitute_(ex_policy, encloser.plastic_continuum_),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      strain_tensor_3D_(encloser.dv_strain_tensor_3D_->DelegatedData(ex_policy)),
      derived_variable_(encloser.dv_derived_variable_->DelegatedData(ex_policy)),
      E_(encloser.E_), nu_(encloser.nu_) {}
//=================================================================================================//
inline void AccDeviatoricPlasticStrainCK::UpdateKernel::update(size_t index_i, Real dt)
{
    Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1.0 / 3.0) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
    Real hydrostatic_pressure = (1.0 / 3.0) * stress_tensor_3D_[index_i].trace();
    Mat3d elastic_strain_tensor_3D = deviatoric_stress / (2.0 * constitute_.getShearModulus(E_, nu_)) +
                                     hydrostatic_pressure * Mat3d::Identity() / (9.0 * constitute_.getBulkModulus(E_, nu_));
    Mat3d plastic_strain_tensor_3D = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D;
    Mat3d deviatoric_strain_tensor = plastic_strain_tensor_3D - (1.0 / (Real)Dimensions) * plastic_strain_tensor_3D.trace() * Mat3d::Identity();
    Real sum = (deviatoric_strain_tensor.cwiseProduct(deviatoric_strain_tensor)).sum();
    derived_variable_[index_i] = sqrt(sum * 2.0 / 3.0);
}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_DYNAMICS_VARIABLE_CK_HPP