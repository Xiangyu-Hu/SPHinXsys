/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	diffussion_reaction.h
 * @brief 	Describe the diffusive and reaction in which 
 *          the dynamics is characterized by diffusion equation and reactive source terms.
 *			Typical physical processes are diffusion, heat conduction 
 *			and chemical and biological reactions. 
 * @author  Xiangyu Hu, Chi Zhang
 */
#pragma once

#include "base_material.h"
#include "solid_particles.h"

#include <map>
#include <functional>
using namespace std::placeholders;

namespace SPH 
{
	template<class BaseParticlesType, class BaseMaterialType>
	class DiffusionReactionParticles;
	class ElectroPhysiologyParticles;
	/**
	 * @class BaseDiffusion
	 * @brief diffusion property abstract base class.
	 */
	class BaseDiffusion
	{
	public:
		BaseDiffusion(size_t diffusion_species_index, size_t gradient_species_index)
		: diffusion_species_index_(diffusion_species_index), gradient_species_index_(gradient_species_index) {};
		virtual ~BaseDiffusion() {};

		size_t diffusion_species_index_;
		size_t gradient_species_index_;

		virtual Real getReferenceDiffusivity() = 0;
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd& direction_from_j_to_i) = 0;
		virtual void initializeLocalDiffusionProperties(BaseParticles* base_particles) {};
		virtual void setupLocalDiffusionProperties(StdVec<Vecd>& material_fiber) {};
	};

	/**
	 * @class IsotropicDiffusion
	 * @brief isotropic diffusion property.
	 */
	class IsotropicDiffusion : public BaseDiffusion
	{
	protected:
		Real diff_cf_; /**< diffusion coefficient. */

	public:
		IsotropicDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf = 1.0) : BaseDiffusion(diffusion_species_index, gradient_species_index),
			diff_cf_(diff_cf) {};
		virtual ~IsotropicDiffusion() {};

		virtual Real getReferenceDiffusivity() override { return diff_cf_; };
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd& direction_from_j_to_i) override
		{
			return diff_cf_;
		};
	};

	/**
	 * @class DirectionalDiffusion
	 * @brief Diffussion is biased along a specific direction.
	 */
	class DirectionalDiffusion : public IsotropicDiffusion
	{
	protected:
		Vecd bias_direction_; 	/**< Reference bias direction. */
		Real bias_diff_cf_;		/**< The bias diffusion coefficient along the fiber direction. */
		Matd transformed_diffusivity_;	/**< The transformed diffusivity with inverse Cholesky decomposition. */

		void initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
	public:
		DirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf, Real bias_diff_cf, Vecd bias_direction) : 
			IsotropicDiffusion(diffusion_species_index, gradient_species_index, diff_cf), 
			bias_direction_(bias_direction), bias_diff_cf_(bias_diff_cf), 
			transformed_diffusivity_(1.0) 
		{
			initializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
		};
		virtual ~DirectionalDiffusion() {};

		virtual Real getReferenceDiffusivity() override
		{
			return SMAX(diff_cf_, diff_cf_ + bias_diff_cf_);
		};

		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, 
			size_t particle_index_j, Vecd& inter_particle_direction) override
		{
			Vecd grad_ij = transformed_diffusivity_ * inter_particle_direction;
			return 1.0 / grad_ij.scalarNormSqr();
		};
	};

	/**
	 * @class LocalDirectionalDiffusion
	 * @brief Diffusion is biased along a specific direction.
	 */
	class LocalDirectionalDiffusion : public DirectionalDiffusion
	{
	protected:
		StdVec<Vecd> local_bias_direction_;
		StdVec<Matd> local_transformed_diffusivity_;
	public:
		LocalDirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
			: DirectionalDiffusion(diffusion_species_index, gradient_species_index, diff_cf, bias_diff_cf, bias_direction) {};
		virtual ~LocalDirectionalDiffusion() {};
		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, size_t particle_index_j, Vecd& inter_particle_direction) override
		{
			Matd trans_diffusivity = getAverageValue(local_transformed_diffusivity_[particle_index_i], local_transformed_diffusivity_[particle_index_j]);
			Vecd grad_ij = trans_diffusivity * inter_particle_direction;
			return 1.0 / grad_ij.scalarNormSqr();
		};
		virtual void initializeLocalDiffusionProperties(BaseParticles* base_particles) override;
		virtual void setupLocalDiffusionProperties(StdVec<Vecd>& material_fiber) override;
	};

	/** Reaction functor . */
	typedef std::function<Real(StdVec<StdLargeVec<Real>>&, size_t particle_i)> ReactionFunctor;
	/**
	 * @class BaseReactionModel
	 * @brief Base class for all reaction models.
	 */
	class BaseReactionModel
	{
	protected:
		virtual void assignDerivedReactionParameters() = 0;
	public:
		BaseReactionModel() {};
		virtual ~BaseReactionModel() {};

		IndexVector reactive_species_;
		StdVec<ReactionFunctor> get_production_rates_;
		StdVec<ReactionFunctor> get_loss_rates_;
	};

	/**
	 * @class AlievPanfilowModel
 	 * @brief The simplest Electrophysiology Reaction model,
	 * which reduces the complex of array of ion currents to two variables that
	 * describe excitation and recovery.
	 */
	class ElectroPhysiologyReaction : public BaseReactionModel
	{
	protected:
		Real k_a_;
		size_t voltage_;
		size_t gate_variable_;
		size_t active_contraction_stress_;

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i) = 0;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i) = 0;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i) = 0;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i) = 0;
		virtual Real getProductionActiveContractionStress(StdVec<StdLargeVec<Real>>& species, size_t particle_i);
		virtual Real getLossRateActiveContractionStress(StdVec<StdLargeVec<Real>>& species, size_t particle_i);
		virtual void assignDerivedReactionParameters() override {};
	public:
		ElectroPhysiologyReaction() : BaseReactionModel(), k_a_(1.0),
			voltage_(0), gate_variable_(1), active_contraction_stress_(2) {};
		virtual ~ElectroPhysiologyReaction() {};
		void initializeElectroPhysiologyReaction(size_t voltage,
			size_t gate_variable, size_t active_contraction_stress);
	};

	class AlievPanfilowModel : public ElectroPhysiologyReaction
	{
	protected:
		/** Parameters for two variable cell model. */
		Real k_, a_, b_, mu_1_, mu_2_, epsilon_, c_m_;

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i) override;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i) override;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i) override;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i) override;
		virtual void assignDerivedReactionParameters() override 
		{
			ElectroPhysiologyReaction::assignDerivedReactionParameters();
		};
	public:
		AlievPanfilowModel() : ElectroPhysiologyReaction(),
			k_(0.0), a_(0.0), b_(0.0), mu_1_(0.0), mu_2_(0.0),
			epsilon_(0.0), c_m_(0.0) {};
		virtual ~AlievPanfilowModel() {};
	};
    
	/**
	 * @class DiffusionReactionMaterial
	 * @brief Complex material for diffusion or/and reactions.
	 */
	template<class BaseParticlesType = BaseParticles, class BaseMaterialType = BaseMaterial>
	class DiffusionReactionMaterial : public BaseMaterialType
	{
	protected:

		size_t number_of_species_;
		size_t number_of_diffusion_species_;
		map<string, size_t> species_indexes_map_;
		StdVec<BaseDiffusion*> species_diffusion_;
		BaseReactionModel* species_reaction_;
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* diffusion_reaction_particles_;

		virtual void assignDerivedMaterialParameters() override 
		{
			BaseMaterialType::assignDerivedMaterialParameters();
		};

		void insertASpecies(string species_name) 
		{
			species_indexes_map_.insert(make_pair(species_name, number_of_species_));
			number_of_species_++;
		};
	public:
		/** Constructor for material only with diffusion. */
		DiffusionReactionMaterial() 
			: BaseMaterialType(), number_of_species_(0), species_reaction_(NULL) 
		{
			BaseMaterialType::material_name_ = "DiffusionMaterial";
		};
		/** Constructor for material with diffusion and reaction. */
		DiffusionReactionMaterial(BaseReactionModel* species_reaction)
			: BaseMaterialType(), number_of_species_(0), species_reaction_(species_reaction) 
		{
			BaseMaterialType::material_name_ = "DiffusionReactionMaterial";
		};
		virtual ~DiffusionReactionMaterial() {};

		size_t NumberOfSpecies() { return number_of_species_; };
		size_t NumberOfSpeciesDiffusion() { return species_diffusion_.size(); };
		StdVec<BaseDiffusion*> SpeciesDiffusion() { return species_diffusion_; };
		BaseReactionModel* SpeciesReaction() { return species_reaction_; };
		map<string, size_t> SpeciesIndexMap() { return  species_indexes_map_; };
		void assignDiffusionReactionParticles(DiffusionReactionParticles<BaseParticlesType, 
			BaseMaterialType>* diffusion_reaction_particles) 
		{
			diffusion_reaction_particles_ = diffusion_reaction_particles;
			for (size_t k = 0; k < species_diffusion_.size(); ++k)
				species_diffusion_[k]->initializeLocalDiffusionProperties(diffusion_reaction_particles);
		};
		/**
		 * @brief Get diffusion time step size. Here, I follow the reference:
		 * https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf 
		 */
		Real getDiffusionTimeStepSize(Real smoothing_length) 
		{
			Real diff_coff_max = 0.0;
			for (size_t k = 0; k < species_diffusion_.size(); ++k)
				diff_coff_max = SMAX(diff_coff_max, species_diffusion_[k]->getReferenceDiffusivity());
			Real dimension = Real(Vecd(0).size());
			return 0.5 * smoothing_length * smoothing_length / diff_coff_max / dimension;
		};
		/** Initialize diffusion material. */
		virtual void initializeDiffusion() = 0;
		virtual DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>* 
			pointToThisObject() override { return this; };	
	};

	/**
	 * @class MonoFieldElectroPhysiology
	 * @brief material class for electro_physiology.
	 */
	class MonoFieldElectroPhysiology 
		: public DiffusionReactionMaterial<SolidParticles, Solid>
	{
	protected:
		Real diff_cf_;
		Real bias_diff_cf_; 
		Vecd bias_direction_;

		virtual void assignDerivedMaterialParameters() override
		{
			DiffusionReactionMaterial<SolidParticles, Solid>::assignDerivedMaterialParameters();
		};
	public:
		MonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction);
		virtual ~MonoFieldElectroPhysiology() {};

		virtual void initializeDiffusion() override;
	};

	/**
	 * @class LocalMonoFieldElectroPhysiology
	 * @brief material class for electro_physiology with locally oriented fibers.
	 */
	class LocalMonoFieldElectroPhysiology 
		: public MonoFieldElectroPhysiology
	{
	protected:
		virtual void assignDerivedMaterialParameters() override
		{
			MonoFieldElectroPhysiology::assignDerivedMaterialParameters();
		};
	public:
		LocalMonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction)
			:MonoFieldElectroPhysiology(electro_physiology_reaction) {
			MonoFieldElectroPhysiology::material_name_ = "LocalMonoFieldElectroPhysiology";
		};
		virtual ~LocalMonoFieldElectroPhysiology() {};

		virtual void initializeDiffusion() override;
		void assignFiberProperties(StdVec<Vecd> &material_fiber);
	};
}
