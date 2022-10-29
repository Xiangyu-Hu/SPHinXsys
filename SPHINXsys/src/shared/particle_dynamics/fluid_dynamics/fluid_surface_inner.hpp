/**
 * @file 	fluid_surface_inner.hpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "fluid_surface_inner.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		template <typename... ConstructorArgs>
		SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			SpatialTemporalFreeSurfaceIdentification(ConstructorArgs &&...args)
			: FreeSurfaceIdentification(std::forward<ConstructorArgs>(args)...)
		{
			this->particles_->registerVariable(previous_surface_indicator_, "PreviousSurfaceIndicator", 1);
			this->particles_->template registerSortableVariable<int>("PreviousSurfaceIndicator");
		}
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIdentification::interaction(index_i, dt);

			if (this->pos_div_[index_i] < this->threshold_by_dimensions_)
			{
				checkNearPreviousFreeSurface(index_i);
			}
		}
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			checkNearPreviousFreeSurface(size_t index_i)
		{
			if (previous_surface_indicator_[index_i] != 1)
			{
				bool is_near_previous_surface = false;
				const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					if (previous_surface_indicator_[inner_neighborhood.j_[n]] == 1)
					{
						is_near_previous_surface = true;
					}
				}
				if (!is_near_previous_surface)
				{
					this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
				}
			}
		}
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			update(size_t index_i, Real dt)
		{
			FreeSurfaceIdentification::update(index_i, dt);

			previous_surface_indicator_[index_i] = this->surface_indicator_[index_i];
		}
		//=================================================================================================//
		template <class DensitySummationType>
		void DensitySummationFreeSurface<DensitySummationType>::update(size_t index_i, Real dt)
		{
			this->rho_[index_i] = ReinitializedDensity(this->rho_sum_[index_i], this->rho0_, this->rho_[index_i]);
		}
		//=================================================================================================//
		template <class DensitySummationFreeSurfaceType>
		template <typename... ConstructorArgs>
		DensitySummationFreeStream<DensitySummationFreeSurfaceType>::
			DensitySummationFreeStream(ConstructorArgs &&...args)
			: DensitySummationFreeSurfaceType(std::forward<ConstructorArgs>(args)...),
			  surface_indicator_(*this->particles_->template getVariableByName<int>("SurfaceIndicator")){};
		//=================================================================================================//
		template <class DensitySummationFreeSurfaceType>
		void DensitySummationFreeStream<DensitySummationFreeSurfaceType>::update(size_t index_i, Real dt)
		{
			if (this->rho_sum_[index_i] < this->rho0_ && isNearSurface(index_i))
			{
				this->rho_[index_i] = this->ReinitializedDensity(this->rho_sum_[index_i], this->rho0_, this->rho_[index_i]);
			}
			else
			{
				this->rho_[index_i] = this->rho_sum_[index_i];
			}
		}
		//=================================================================================================//
		template <class DensitySummationFreeSurfaceType>
		bool DensitySummationFreeStream<DensitySummationFreeSurfaceType>::isNearSurface(size_t index_i)
		{
			bool is_near_surface = true;
			if (surface_indicator_[index_i] != 1)
			{
				is_near_surface = false;
				const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					if (surface_indicator_[inner_neighborhood.j_[n]] == 1)
					{
						is_near_surface = true;
						break;
					}
				}
			}
			return is_near_surface;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
