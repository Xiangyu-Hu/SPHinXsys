/**
 * @file 	fluid_dynamics_complex.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"
#include "fluid_dynamics_inner.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(BaseBodyRelationInner &inner_relation,
										 BaseBodyRelationContact &contact_relation, Real thereshold)
			: FreeSurfaceIndicationInner(inner_relation, thereshold), FluidContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				Real rho0_k = contact_particles_[k]->rho0_;
				contact_inv_rho0_.push_back(1.0 / rho0_k);
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(ComplexBodyRelation &complex_relation, Real thereshold)
			: FreeSurfaceIndicationComplex(complex_relation.inner_relation_,
										   complex_relation.contact_relation_, thereshold) {}
		//=================================================================================================//
		void FreeSurfaceIndicationComplex::Interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIndicationInner::Interaction(index_i, dt);

			Real pos_div = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
				Real contact_inv_rho0_k = contact_inv_rho0_[k];
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					pos_div -= contact_neighborhood.dW_ij_[n] * contact_neighborhood.r_ij_[n] *
							   contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
				}
			}
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(BaseBodyRelationInner &inner_relation,
											   BaseBodyRelationContact &contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(inner_relation, contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation)
			: TransportVelocityCorrectionComplex(complex_relation.inner_relation_,
												 complex_relation.contact_relation_) {}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation,
											   BaseBodyRelationContact &extra_contact_relation)
			: ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(complex_relation, extra_contact_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::prepareContactData()
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::Interaction(size_t index_i, Real dt)
		{
			TransportVelocityCorrectionInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];

			Vecd acceleration_trans(0);
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					//acceleration for transport velocity
					acceleration_trans -= 2.0 * p_background_ * Vol_k[index_j] * nablaW_ij / rho_i;
				}
			}

			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
				pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationRiemannWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxation<PressureRelaxationRiemannInnerOldroyd_B>::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ij * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		void DensityRelaxationRiemannWithWallOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxation<DensityRelaxationRiemannInnerOldroyd_B>::Interaction(index_i, dt);

			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(wall_Vol_[k]);
				StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ij = wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];

					Matd velocity_gradient = -SimTK::outer((vel_i - vel_ave_k[index_j]), nablaW_ij) * Vol_k[index_j] * 2.0;
					stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
								   (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
				}
			}
			dtau_dt_[index_i] += stress_rate;
		}
		//=================================================================================================//
		ColorFunctionGradientComplex::ColorFunctionGradientComplex(BaseBodyRelationInner &inner_relation,
																   BaseBodyRelationContact &contact_relation)
			: ColorFunctionGradientInner(inner_relation), FluidContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		ColorFunctionGradientComplex::ColorFunctionGradientComplex(ComplexBodyRelation &complex_relation)
			: ColorFunctionGradientComplex(complex_relation.inner_relation_,
										   complex_relation.contact_relation_) {}
		//=================================================================================================//
		void ColorFunctionGradientComplex::Interaction(size_t index_i, Real dt)
		{
			ColorFunctionGradientInner::Interaction(index_i, dt);

			Vecd gradient(0.0);
			if (pos_div_[index_i] < thereshold_by_dimensions_)
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real> &contact_vol_k = *(contact_Vol_[k]);
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						gradient -= contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n] *
									contact_vol_k[contact_neighborhood.j_[n]];
					}
				}
			}
			color_grad_[index_i] += gradient;
			surface_norm_[index_i] = color_grad_[index_i] / (color_grad_[index_i].norm() + TinyReal);
		}
		//=================================================================================================//
		SurfaceNormWithWall::SurfaceNormWithWall(BaseBodyRelationContact &contact_relation, Real contact_angle)
			: InteractionDynamics(*contact_relation.sph_body_), FSIContactData(contact_relation),
			  contact_angle_(contact_angle),
			  surface_indicator_(particles_->surface_indicator_),
			  surface_norm_(*particles_->getVariableByName<indexVector, Vecd>("SurfaceNormal")),
			  pos_div_(*particles_->getVariableByName<indexScalar, Real>("PositionDivergence"))
		{
			particle_spacing_ = contact_relation.sph_body_->sph_adaptation_->ReferenceSpacing();
			smoothing_length_ = contact_relation.sph_body_->sph_adaptation_->ReferenceSmoothingLength();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				wall_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void SurfaceNormWithWall::Interaction(size_t index_i, Real dt)
		{
			Real large_dist(1.0e6);
			Vecd n_i = surface_norm_[index_i];
			Real smoothing_factor(1.0);
			Vecd smooth_norm(0);
			Vecd n_i_w(0);
			/** Contact interaction. */
			if (surface_indicator_[index_i] == 1)
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
					Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
					{
						size_t index_j = wall_neighborhood.j_[n];
						if (wall_neighborhood.r_ij_[n] < large_dist)
						{
							Vecd n_w_t = n_i - dot(n_i, n_k[index_j]) * n_k[index_j];
							Vecd n_t = n_w_t / (n_w_t.norm() + TinyReal);
							n_i_w = n_t * sin(contact_angle_) + cos(contact_angle_) * n_k[index_j];
							/** No change for multi-resolution. */
							Real r_ij = wall_neighborhood.r_ij_[n] * dot(n_k[index_j], wall_neighborhood.e_ij_[n]);
							if (r_ij <= smoothing_length_)
							{
								smoothing_factor = 0.0;
							}
							else
							{
								smoothing_factor = (r_ij - smoothing_length_) / smoothing_length_;
							}
							large_dist = wall_neighborhood.r_ij_[n];
							smooth_norm = smoothing_factor * n_i + (1.0 - smoothing_factor) * n_i_w;
							surface_norm_[index_i] = smooth_norm / (smooth_norm.norm() + TinyReal);
						}
					}
				}
			}
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//