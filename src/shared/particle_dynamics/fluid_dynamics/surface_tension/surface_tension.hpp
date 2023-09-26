#pragma once

#include "surface_tension.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
ColorFunctionGradient<DataDelegationType>::ColorFunctionGradient(BaseRelationType &base_relation)
    : FreeSurfaceIndication<DataDelegationType>(base_relation)
{
    this->particles_->registerVariable(color_grad_, "ColorGradient");
    this->particles_->registerVariable(surface_norm_, "SurfaceNormal");
}
//=================================================================================================//
void SurfaceTensionAccelerationInner::
    interaction(size_t index_i, Real dt)
{
    Vecd n_i = surface_norm_[index_i];
    Real curvature(0.0);
    Real renormalized_curvature(0);
    Real pos_div(0);
    if (indicator_[index_i] == 1)
    {
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (indicator_[index_j] == 1)
            {
                Vecd n_j = surface_norm_[index_j];
                Vecd n_ij = n_i - n_j;
                curvature -= inner_neighborhood.dW_ijV_j_[n] * n_ij.dot(inner_neighborhood.e_ij_[n]);
                pos_div -= inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.r_ij_[n];
            }
        }
    }
    /**
     Adami et al. 2010 has a typo in equation.
     (dv / dt)_s = (1.0 / rho) (-sigma * k * n * delta)
                             = (1/rho) * curvature * color_grad
                             = (1/m) * curvature * color_grad * vol
     */
    renormalized_curvature = (Real)Dimensions * curvature / ABS(pos_div + TinyReal);
    Vecd acceleration = gamma_ * renormalized_curvature * color_grad_[index_i] * Vol_[index_i];
    acc_prior_[index_i] -= acceleration / mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
