/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	diffusion_reaction_particles.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef	TWO_PHASE_HEAT_TRANSFER_PARTICLES_H
#define TWO_PHASE_HEAT_TRANSFER_PARTICLES_H

//#include "base_particles.h"
//#include "base_body.h"
//#include "base_material.h"
//#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"
namespace SPH
{

	template <class BaseParticlesType, class BaseMaterialType = BaseMaterial, int NUM_SPECIES = 1>
	class TwoPhaseDiffusionReactionParticles : public DiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	public:
		StdLargeVec<Real> external_diffusion_dt_; /**< array of the time derivative of diffusion species */
		StdLargeVec<Real> external_diffusion_dt_sum_; /**< array of the time derivative of diffusion species */
		StdLargeVec<Real> thermal_conductivity_; /**< array of the time derivative of diffusion species */

		TwoPhaseDiffusionReactionParticles(SPHBody& sph_body,
			DiffusionReaction<BaseMaterialType, NUM_SPECIES>* diffusion_reaction_material)
			: DiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(sph_body, diffusion_reaction_material)
		{
			thermal_conductivity_.resize(this->number_of_diffusion_species_);
			external_diffusion_dt_.resize(this->number_of_diffusion_species_);
			external_diffusion_dt_sum_.resize(this->number_of_diffusion_species_);
		};
		virtual ~TwoPhaseDiffusionReactionParticles() {};


		virtual void initializeOtherVariables() override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::initializeOtherVariables();

			this->registerVariable(external_diffusion_dt_sum_, "HeatFlux");
			this->registerVariable(thermal_conductivity_, "ThermalConductivity");
			this->template addVariableToWrite<Real>("HeatFlux");
			this->template addVariableToWrite<Real>("ThermalConductivity");
			for (size_t m = 0; m < this->number_of_diffusion_species_; ++m)
			{
				constexpr int type_index = DataTypeIndex<Real>::value;
				/**
				 * register reactive change rate terms without giving variable name
				 */
				std::get<type_index>(this->all_particle_data_).push_back(&external_diffusion_dt_);
				external_diffusion_dt_.resize(this->real_particles_bound_, Real(0));
			}
		};

		virtual TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>* ThisObjectPtr() override { return this; };
	};
}
#endif // TWO_PHASE_HEAT_TRANSFER_PARTICLES_H