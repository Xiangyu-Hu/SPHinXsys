#include "fluid_surface_complex.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(BaseInnerRelation &inner_relation,
										 BaseContactRelation &contact_relation, Real threshold)
			: FreeSurfaceIndicationInner(inner_relation, threshold), FluidContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
				contact_inv_rho0_.push_back(1.0 / rho0_k);
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
			}
		}
		//=================================================================================================//
		FreeSurfaceIndicationComplex::
			FreeSurfaceIndicationComplex(ComplexRelation &complex_relation, Real threshold)
			: FreeSurfaceIndicationComplex(complex_relation.getInnerRelation(),
										   complex_relation.getContactRelation(), threshold) {}
		//=================================================================================================//
		void FreeSurfaceIndicationComplex::interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIndicationInner::interaction(index_i, dt);

			Real pos_div = 0.0;
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					pos_div -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.r_ij_[n];
				}
			}
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
		ColorFunctionGradientComplex::ColorFunctionGradientComplex(BaseInnerRelation &inner_relation,
																   BaseContactRelation &contact_relation)
			: ColorFunctionGradientInner(inner_relation), FluidContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		ColorFunctionGradientComplex::ColorFunctionGradientComplex(ComplexRelation &complex_relation)
			: ColorFunctionGradientComplex(complex_relation.getInnerRelation(),
										   complex_relation.getContactRelation()) {}
		//=================================================================================================//
		void ColorFunctionGradientComplex::interaction(size_t index_i, Real dt)
		{
			ColorFunctionGradientInner::interaction(index_i, dt);

			Vecd gradient = Vecd::Zero();
			if (pos_div_[index_i] < threshold_by_dimensions_)
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						gradient -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
					}
				}
			}
			color_grad_[index_i] += gradient;
			surface_norm_[index_i] = color_grad_[index_i] / (color_grad_[index_i].norm() + TinyReal);
		}
		//=================================================================================================//
		SurfaceNormWithWall::SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle)
			: LocalDynamics(contact_relation.sph_body_), FSIContactData(contact_relation),
			  contact_angle_(contact_angle),
			  surface_indicator_(particles_->surface_indicator_),
			  surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")),
			  pos_div_(*particles_->getVariableByName<Real>("PositionDivergence"))
		{
			particle_spacing_ = contact_relation.sph_body_.sph_adaptation_->ReferenceSpacing();
			smoothing_length_ = contact_relation.sph_body_.sph_adaptation_->ReferenceSmoothingLength();
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				wall_n_.push_back(&(contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		void SurfaceNormWithWall::interaction(size_t index_i, Real dt)
		{
			Real large_dist(1.0e6);
			Vecd n_i = surface_norm_[index_i];
			Real smoothing_factor(1.0);
			Vecd smooth_norm = Vecd::Zero();
			Vecd n_i_w = Vecd::Zero();
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
							Vecd n_w_t = n_i - n_i.dot(n_k[index_j]) * n_k[index_j];
							Vecd n_t = n_w_t / (n_w_t.norm() + TinyReal);
							n_i_w = n_t * sin(contact_angle_) + cos(contact_angle_) * n_k[index_j];
							/** No change for multi-resolution. */
							Real r_ij = wall_neighborhood.r_ij_[n] * n_k[index_j].dot(wall_neighborhood.e_ij_[n]);
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