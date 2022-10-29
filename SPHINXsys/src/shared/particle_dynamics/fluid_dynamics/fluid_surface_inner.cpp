/**
 * @file 	fluid_surface_inner.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "fluid_surface_inner.hpp"

namespace SPH
{
    //=================================================================================================//
    namespace fluid_dynamics
    {
        //=================================================================================================//
        FreeSurfaceIndicationInner::
            FreeSurfaceIndicationInner(BaseInnerRelation &inner_relation, Real threshold)
            : LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
              threshold_by_dimensions_(threshold * (Real)Dimensions),
              surface_indicator_(particles_->surface_indicator_),
              smoothing_length_(inner_relation.sph_body_.sph_adaptation_->ReferenceSmoothingLength())
        {
            particles_->registerVariable(pos_div_, "PositionDivergence");
        }
        //=================================================================================================//
        void FreeSurfaceIndicationInner::interaction(size_t index_i, Real dt)
        {
            Real pos_div = 0.0;
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                pos_div -= inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.r_ij_[n];
            }
            pos_div_[index_i] = pos_div;
        }
        //=================================================================================================//
        void FreeSurfaceIndicationInner::update(size_t index_i, Real dt)
        {
            bool is_free_surface = pos_div_[index_i] < threshold_by_dimensions_ ? true : false;

            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                /** Two layer particles.*/
                if (pos_div_[inner_neighborhood.j_[n]] < threshold_by_dimensions_ &&
                    inner_neighborhood.r_ij_[n] < smoothing_length_)
                {
                    is_free_surface = true;
                    break;
                }
            }
            surface_indicator_[index_i] = is_free_surface ? 1 : 0;
        }
        //=================================================================================================//
        FreeStreamBoundaryVelocityCorrection::FreeStreamBoundaryVelocityCorrection(SPHBody &sph_body)
            : LocalDynamics(sph_body), FluidDataSimple(sph_body),
              u_ref_(1.0), t_ref_(2.0), rho_ref_(particles_->fluid_.ReferenceDensity()),
              rho_sum(particles_->rho_sum_), vel_(particles_->vel_),
              surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")) {}
        //=================================================================================================//
        void FreeStreamBoundaryVelocityCorrection::update(size_t index_i, Real dt)
        {
            if (surface_indicator_[index_i] == 1)
            {
                // TODO: free stream condition should be summarized to a separated class.
                Real run_time_ = GlobalStaticVariables::physical_time_;
                Real u_freestream = run_time_ < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time_ / t_ref_)) : u_ref_;
                vel_[index_i][0] = u_freestream + SMIN(rho_sum[index_i], rho_ref_) * (vel_[index_i][0] - u_freestream) / rho_ref_;
            }
        }
        //=================================================================================================//
        ColorFunctionGradientInner::ColorFunctionGradientInner(BaseInnerRelation &inner_relation)
            : LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
              surface_indicator_(particles_->surface_indicator_),
              pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
              threshold_by_dimensions_((0.75 * (Real)Dimensions))
        {
            particles_->registerVariable(color_grad_, "ColorGradient");
            particles_->registerVariable(surface_norm_, "SurfaceNormal");
        }
        //=================================================================================================//
        void ColorFunctionGradientInner::interaction(size_t index_i, Real dt)
        {
            Vecd gradient(0);
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            if (pos_div_[index_i] < threshold_by_dimensions_)
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    gradient -= inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                }
            }
            color_grad_[index_i] = gradient;
            surface_norm_[index_i] = gradient / (gradient.norm() + TinyReal);
        }
        //=================================================================================================//
        ColorFunctionGradientInterpolationInner::ColorFunctionGradientInterpolationInner(BaseInnerRelation &inner_relation)
            : LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation), Vol_(particles_->Vol_),
              surface_indicator_(particles_->surface_indicator_),
              color_grad_(*particles_->getVariableByName<Vecd>("ColorGradient")),
              surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")),
              pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
              threshold_by_dimensions_((0.75 * (Real)Dimensions))

        {
            particles_->addVariableToWrite<Vecd>("SurfaceNormal");
            particles_->addVariableToWrite<Vecd>("ColorGradient");
        }
        //=================================================================================================//
        void ColorFunctionGradientInterpolationInner::interaction(size_t index_i, Real dt)
        {
            Vecd grad(0);
            Real weight(0);
            Real total_weight(0);
            if (surface_indicator_[index_i] == 1 && pos_div_[index_i] > threshold_by_dimensions_)
            {
                Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (surface_indicator_[index_j] == 1 && pos_div_[index_j] < threshold_by_dimensions_)
                    {
                        weight = inner_neighborhood.W_ij_[n] * Vol_[index_j];
                        grad += weight * color_grad_[index_j];
                        total_weight += weight;
                    }
                }
                Vecd grad_norm = grad / (total_weight + TinyReal);
                color_grad_[index_i] = grad_norm;
                surface_norm_[index_i] = grad_norm / (grad_norm.norm() + TinyReal);
            }
        }
        //=================================================================================================//
        SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation, Real gamma)
            : LocalDynamics(inner_relation.sph_body_), FluidDataInner(inner_relation),
              gamma_(gamma), Vol_(particles_->Vol_), mass_(particles_->mass_),
              acc_prior_(particles_->acc_prior_), surface_indicator_(particles_->surface_indicator_),
              color_grad_(*particles_->getVariableByName<Vecd>("ColorGradient")),
              surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")) {}
        //=================================================================================================//
        SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation)
            : SurfaceTensionAccelerationInner(inner_relation, 1.0) {}
        //=================================================================================================//
        void SurfaceTensionAccelerationInner::interaction(size_t index_i, Real dt)
        {
            Vecd n_i = surface_norm_[index_i];
            Real curvature(0.0);
            Real renormalized_curvature(0);
            Real pos_div(0);
            if (surface_indicator_[index_i] == 1)
            {
                Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (surface_indicator_[index_j] == 1)
                    {
                        Vecd n_j = surface_norm_[index_j];
                        Vecd n_ij = n_i - n_j;
                        curvature -= inner_neighborhood.dW_ijV_j_[n] * dot(n_ij, inner_neighborhood.e_ij_[n]);
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
    }
    //=================================================================================================//
}
//=================================================================================================//