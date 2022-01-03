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

#ifndef DIFFUSION_REACTION_H
#define DIFFUSION_REACTION_H

#include "base_material.h"
#include "solid_particles.h"

#include <map>
#include <functional>
using namespace std::placeholders;

namespace SPH
{
	template <class BaseParticlesType, class BaseMaterialType>
	class DiffusionReactionParticles;
	class ElectroPhysiologyParticles;
	/**
	 * @class BaseDiffusion
	 * @brief diffusion property abstract base class.
	 */
	class BaseDiffusion : public BaseMaterial
	{
	public:
		BaseDiffusion(size_t diffusion_species_index, size_t gradient_species_index)
			: BaseMaterial(), diffusion_species_index_(diffusion_species_index),
			  gradient_species_index_(gradient_species_index)
		{
			material_type_ = "BaseDiffusion";
		};
		virtual ~BaseDiffusion(){};

		size_t diffusion_species_index_;
		size_t gradient_species_index_;

		virtual Real getReferenceDiffusivity() = 0;
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd &direction_from_j_to_i) = 0;
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
						   Real diff_cf = 1.0)
			: BaseDiffusion(diffusion_species_index, gradient_species_index),
			  diff_cf_(diff_cf)
		{
			material_type_ = "IsotropicDiffusion";
		};
		virtual ~IsotropicDiffusion(){};

		virtual Real getReferenceDiffusivity() override { return diff_cf_; };
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd &direction_from_j_to_i) override
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
		Vecd bias_direction_;		   /**< Reference bias direction. */
		Real bias_diff_cf_;			   /**< The bias diffusion coefficient along the fiber direction. */
		Matd transformed_diffusivity_; /**< The transformed diffusivity with inverse Cholesky decomposition. */

		void initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction);

	public:
		DirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
							 Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
			: IsotropicDiffusion(diffusion_species_index, gradient_species_index, diff_cf),
			  bias_direction_(bias_direction), bias_diff_cf_(bias_diff_cf),
			  transformed_diffusivity_(1.0)
		{
			material_type_ = "DirectionalDiffusion";
			initializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
		};
		virtual ~DirectionalDiffusion(){};

		virtual Real getReferenceDiffusivity() override
		{
			return SMAX(diff_cf_, diff_cf_ + bias_diff_cf_);
		};

		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i,
												   size_t particle_index_j, Vecd &inter_particle_direction) override
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
		StdLargeVec<Vecd> local_bias_direction_;
		StdLargeVec<Matd> local_transformed_diffusivity_;

		void initializeFiberDirection();

	public:
		LocalDirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
								  Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
			: DirectionalDiffusion(diffusion_species_index, gradient_species_index, diff_cf, bias_diff_cf, bias_direction)
		{
			material_type_ = "LocalDirectionalDiffusion";
		};
		virtual ~LocalDirectionalDiffusion(){};
		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, size_t particle_index_j, Vecd &inter_particle_direction) override
		{
			Matd trans_diffusivity = getAverageValue(local_transformed_diffusivity_[particle_index_i], local_transformed_diffusivity_[particle_index_j]);
			Vecd grad_ij = trans_diffusivity * inter_particle_direction;
			return 1.0 / grad_ij.scalarNormSqr();
		};
		virtual void assignBaseParticles(BaseParticles *base_particles) override;
		virtual void readFromXmlForLocalParameters(const std::string &filefullpath) override;
	};

	/** Reaction functor . */
	typedef std::function<Real(StdVec<StdLargeVec<Real>> &, size_t particle_i)> ReactionFunctor;
	/**
	 * @class BaseReactionModel
	 * @brief Base class for all reaction models.
	 */
	class BaseReactionModel
	{
	protected:
		StdVec<std::string> species_name_list_;
		std::map<std::string, size_t> species_indexes_map_;
		std::string reaction_model_;

	public:
		explicit BaseReactionModel(StdVec<std::string> species_name_list)
			: reaction_model_("BaseReactionModel"), species_name_list_(species_name_list)
		{
			for (size_t i = 0; i != species_name_list.size(); ++i)
			{
				species_indexes_map_.insert(make_pair(species_name_list[i], i));
			}
		};
		virtual ~BaseReactionModel(){};

		IndexVector reactive_species_;
		StdVec<ReactionFunctor> get_production_rates_;
		StdVec<ReactionFunctor> get_loss_rates_;

		StdVec<std::string> getSpeciesNameList() { return species_name_list_; };
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

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getProductionActiveContractionStress(StdVec<StdLargeVec<Real>> &species, size_t particle_i);
		virtual Real getLossRateActiveContractionStress(StdVec<StdLargeVec<Real>> &species, size_t particle_i);

	public:
		explicit ElectroPhysiologyReaction(Real k_a)
			: BaseReactionModel({"Voltage", "GateVariable", "ActiveContractionStress"}),
			  k_a_(k_a), voltage_(species_indexes_map_["Voltage"]),
			  gate_variable_(species_indexes_map_["GateVariable"]),
			  active_contraction_stress_(species_indexes_map_["ActiveContractionStress"])
		{
			reaction_model_ = "ElectroPhysiologyReaction";
			initializeElectroPhysiologyReaction();
		};
		virtual ~ElectroPhysiologyReaction(){};

		void initializeElectroPhysiologyReaction();
	};

	class AlievPanfilowModel : public ElectroPhysiologyReaction
	{
	protected:
		/** Parameters for two variable cell model. */
		Real k_, a_, b_, mu_1_, mu_2_, epsilon_, c_m_;

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;

	public:
		explicit AlievPanfilowModel(Real k_a, Real c_m, Real k, Real a, Real b, Real mu_1, Real mu_2, Real epsilon)
			: ElectroPhysiologyReaction(k_a), k_(k), a_(a), b_(b), mu_1_(mu_1), mu_2_(mu_2),
			  epsilon_(epsilon), c_m_(c_m)
		{
			reaction_model_ = "AlievPanfilowModel";
		};
		virtual ~AlievPanfilowModel(){};
	};

	/**
	 * @class DiffusionReaction
	 * @brief Complex material for diffusion or/and reactions.
	 */
	template <class BaseParticlesType = BaseParticles, class BaseMaterialType = BaseMaterial>
	class DiffusionReaction : public BaseMaterialType
	{
	private:
		UniquePtrVectorKeeper<BaseDiffusion> diffusion_ptr_keeper_;

	protected:
		StdVec<std::string> species_name_list_;
		size_t number_of_species_;
		std::map<std::string, size_t> species_indexes_map_;
		StdVec<BaseDiffusion *> species_diffusion_;
		BaseReactionModel *species_reaction_;
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType> *diffusion_reaction_particles_;

	public:
		/** Constructor for material with diffusion only. */
		template <typename... ConstructorArgs>
		DiffusionReaction(StdVec<std::string> species_name_list, ConstructorArgs &&...args)
			: BaseMaterialType(std::forward<ConstructorArgs>(args)...),
			  species_name_list_(species_name_list),
			  number_of_species_(species_name_list.size()),
			  species_reaction_(nullptr),
			  diffusion_reaction_particles_(nullptr)
		{
			BaseMaterialType::material_type_ = "Diffusion";
			for (size_t i = 0; i != number_of_species_; ++i)
			{
				species_indexes_map_.insert(make_pair(species_name_list[i], i));
			}
		};
		/** Constructor for material with diffusion and reaction. */
		template <typename... ConstructorArgs>
		DiffusionReaction(BaseReactionModel &species_reaction,
						  StdVec<std::string> species_name_list, ConstructorArgs &&...args)
			: DiffusionReaction(species_name_list, std::forward<ConstructorArgs>(args)...)
		{
			species_reaction_ = &species_reaction;
			BaseMaterialType::material_type_ = "DiffusionReaction";
		};
		virtual ~DiffusionReaction(){};

		size_t NumberOfSpecies() { return number_of_species_; };
		size_t NumberOfSpeciesDiffusion() { return species_diffusion_.size(); };
		StdVec<BaseDiffusion *> SpeciesDiffusion() { return species_diffusion_; };
		BaseReactionModel *SpeciesReaction() { return species_reaction_; };
		std::map<std::string, size_t> SpeciesIndexMap() { return species_indexes_map_; };
		StdVec<std::string> getSpeciesNameList() { return species_name_list_; };
		void assignDiffusionReactionParticles(DiffusionReactionParticles<BaseParticlesType,
																		 BaseMaterialType> *diffusion_reaction_particles)
		{
			diffusion_reaction_particles_ = diffusion_reaction_particles;
			for (size_t k = 0; k < species_diffusion_.size(); ++k)
				species_diffusion_[k]->assignBaseParticles(diffusion_reaction_particles);
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

		/** Initialize a diffusion material. */
		template <class DiffusionType, typename... ConstructorArgs>
		void initializeAnDiffusion(const std::string &diffusion_species_name,
								   const std::string &gradient_species_name, ConstructorArgs &&...args)
		{
			species_diffusion_.push_back(
				diffusion_ptr_keeper_.createPtr<DiffusionType>(
					species_indexes_map_[diffusion_species_name],
					species_indexes_map_[diffusion_species_name], std::forward<ConstructorArgs>(args)...));
		};

		virtual DiffusionReaction<BaseParticlesType, BaseMaterialType> *
		ThisObjectPtr() override { return this; };
	};

	/**
	 * @class MonoFieldElectroPhysiology
	 * @brief material class for electro_physiology.
	 */
	class MonoFieldElectroPhysiology
		: public DiffusionReaction<SolidParticles, Solid>
	{
	public:
		explicit MonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
											Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
		virtual ~MonoFieldElectroPhysiology(){};
	};

	/**
	 * @class LocalMonoFieldElectroPhysiology
	 * @brief material class for electro_physiology with locally oriented fibers.
	 */
	class LocalMonoFieldElectroPhysiology
		: public DiffusionReaction<SolidParticles, Solid>
	{
	public:
		explicit LocalMonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
												 Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
		virtual ~LocalMonoFieldElectroPhysiology(){};

		virtual void readFromXmlForLocalParameters(const std::string &filefullpath) override;
	};
}
#endif //DIFFUSION_REACTION_H