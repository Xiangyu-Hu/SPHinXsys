/**
 * @file 	fluid_surface_inner.hpp
 * @author	Chi ZHang and Xiangyu Hu
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
			this->particles_->registerAVariable(previous_surface_indicator_, "PreviousSurfaceIndicator", 1);
			this->particles_->template registerASortableVariable<int>("PreviousSurfaceIndicator");
		}
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			Interaction(size_t index_i, Real dt)
		{
			FreeSurfaceIdentification::Interaction(index_i, dt);

			if (this->pos_div_[index_i] < this->thereshold_by_dimensions_)
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
					this->pos_div_[index_i] = 2.0 * this->thereshold_by_dimensions_;
				}
			}
		}
		//=================================================================================================//
		template <class FreeSurfaceIdentification>
		void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentification>::
			Update(size_t index_i, Real dt)
		{
			FreeSurfaceIdentification::Update(index_i, dt);

			previous_surface_indicator_[index_i] = this->surface_indicator_[index_i];
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
