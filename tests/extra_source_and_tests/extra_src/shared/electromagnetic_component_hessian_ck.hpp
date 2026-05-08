#ifndef ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_HPP
#define ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_HPP

#include "electromagnetic_component_hessian_ck.h"

namespace SPH
{
namespace electromagnetics
{
//=================================================================================================//
CopyVectorFieldComponentsToScalarVariables::
    CopyVectorFieldComponentsToScalarVariables(
        SPHBody &sph_body,
        const std::string &vector_field_name,
        const ComponentVariableNames<> &component_variable_names)
    : LocalDynamics(sph_body),
      vector_field_(particles_->getVariableDataByName<Vecd>(vector_field_name))
{
    for (size_t axis = 0; axis != Components; ++axis)
    {
        component_variables_[axis] =
            particles_->getVariableDataByName<Real>(component_variable_names[axis]);
    }
}
//=================================================================================================//
void CopyVectorFieldComponentsToScalarVariables::update(size_t index_i, Real dt)
{
    (void)dt;
    for (size_t axis = 0; axis != Components; ++axis)
    {
        component_variables_[axis][index_i] = vector_field_[index_i][axis];
    }
}
//=================================================================================================//
ReconstructCurlNuBFromScalarComponentHessians::
    ReconstructCurlNuBFromScalarComponentHessians(
        SPHBody &sph_body,
        const ComponentVariableNames<> &component_variable_names,
        Real second_order_operator_scaling,
        const std::string &output_name,
        const std::string &magnetic_reluctivity_name)
    : LocalDynamics(sph_body),
      second_order_operator_scaling_(second_order_operator_scaling),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>(magnetic_reluctivity_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(output_name))
{
    for (size_t axis = 0; axis != Components; ++axis)
    {
        component_hessians_[axis] =
            particles_->getVariableDataByName<VecMatd>(component_variable_names[axis] + "Hessian");
    }
}
//=================================================================================================//
void ReconstructCurlNuBFromScalarComponentHessians::update(size_t index_i, Real dt)
{
    (void)dt;
    Vecd grad_div = ZeroData<Vecd>::value;
    Vecd laplacian = ZeroData<Vecd>::value;

    if constexpr (Components == 2)
    {
        const VecMatd &h_ax = component_hessians_[0][index_i];
        const VecMatd &h_ay = component_hessians_[1][index_i];
        grad_div[0] = h_ax[0] + h_ay[2];
        grad_div[1] = h_ax[2] + h_ay[1];
        laplacian[0] = h_ax[0] + h_ax[1];
        laplacian[1] = h_ay[0] + h_ay[1];
    }
    else if constexpr (Components == 3)
    {
        const VecMatd &h_ax = component_hessians_[0][index_i];
        const VecMatd &h_ay = component_hessians_[1][index_i];
        const VecMatd &h_az = component_hessians_[2][index_i];
        grad_div[0] = h_ax[0] + h_ay[3] + h_az[5];
        grad_div[1] = h_ax[3] + h_ay[1] + h_az[4];
        grad_div[2] = h_ax[5] + h_ay[4] + h_az[2];
        laplacian[0] = h_ax[0] + h_ax[1] + h_ax[2];
        laplacian[1] = h_ay[0] + h_ay[1] + h_ay[2];
        laplacian[2] = h_az[0] + h_az[1] + h_az[2];
    }

    curl_nu_b_[index_i] =
        magnetic_reluctivity_[index_i] * second_order_operator_scaling_ * (grad_div - laplacian);
}
//=================================================================================================//
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_HPP
