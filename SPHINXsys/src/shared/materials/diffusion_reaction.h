/**
 * @file 	difussion_reaction.h
 * @brief 	Desrcibe the dffusive and reaction in which 
 *          the dynamcis is characterized duffusion equation and reactive source terms.
 *			Typical physical processes are diffusion, heat condution 
 *			and chemical and bioligical reactions. 
 * @author  Xiangyu Hu, Chi Zhang
 * @version 0.2.0
 * 			Here, Version 0.2.0 is set, as reaction-diffusion equation is solved. 
 */
#pragma once

#include "base_material.h"
#include "solid_particles.h"
#include <map>
#include <functional>
using namespace std::placeholders;

namespace SPH 
{
	/** preclaimed classes */
	template<class BaseParticlesType, class BaseMaterialType>
	class DiffusionReactionParticles;
	class ElectroPhysiologyParticles;
	/**
	 * @class BaseDiffusion
	 * @brief diffusision property abstract base class.
	 */
	class BaseDiffusion
	{
	public:
		/** Default constructor. */
		BaseDiffusion(size_t diffusion_species_index, size_t gradient_species_index)
		: diffusion_species_index_(diffusion_species_index), gradient_species_index_(gradient_species_index) {};
		virtual ~BaseDiffusion() {};

		/** Index of the species with diffusion transport. */
		size_t diffusion_species_index_;
		/** Index of the species whose gradient guides the diffusion. */
		size_t gradient_species_index_;

		/** Get the reference diffusivity for time step criterion. */
		virtual Real getReferenceDiffusivity() = 0;
		/** Get diffusion coefficient along the interacting particle direction. */
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd& direction_from_j_to_i) = 0;
		/** initialize the local property. */
		virtual void initializeLocalDiffusionProperties(BaseParticles* base_particles) {};
		/** Setup the local property after initialization. */
		virtual void setupLocalDiffusionProperties(StdVec<Vecd>& material_fiber) {};
	};

	/**
	 * @class IsotropicDiffusion
	 * @brief isotropic diffusision property.
	 */
	class IsotropicDiffusion : public BaseDiffusion
	{
	protected:
		Real diff_cf_; /*> diffusion coefficient. */

	public:
		/** Constructor. */
		IsotropicDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf) : BaseDiffusion(diffusion_species_index, gradient_species_index),
			diff_cf_(diff_cf) {};
		virtual ~IsotropicDiffusion() {};

		/** Get the reference diffusivity. */
		virtual Real getReferenceDiffusivity() override { return diff_cf_; };
		/** Get diffusion coefficient along the interacting particle direction. */
		virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd& direction_from_j_to_i) override
		{
			return diff_cf_;
		};
	};

	/**
	 * @class DirectionalDiffusion
	 * @brief Diffusision is biased along a specific direction.
	 */
	class DirectionalDiffusion : public IsotropicDiffusion
	{
	protected:
		Vecd bias_direction_; /*> Reference bias direction. */
		/*> The bias diffusion coefficient along the fiber direction. */
		Real bias_diff_cf_;
		/*> The transformed diffusivity with inverse Cholesky decomposition. */
		Matd transf_diffusivity_;
		/** Intialize directional diffusivity. */
		void intializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
	public:
		/** Constructor*/
		DirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf, Real bias_diff_cf, Vecd bias_direction) 
			: IsotropicDiffusion(diffusion_species_index, gradient_species_index,
				diff_cf), bias_diff_cf_(bias_diff_cf), bias_direction_(bias_direction),
			transf_diffusivity_(1.0) 
		{
			intializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
		};
		virtual ~DirectionalDiffusion() {};

		/** Get the reference diffusivity. */
		virtual Real getReferenceDiffusivity() override
		{
			return SMAX(diff_cf_, diff_cf_ + bias_diff_cf_);
		};
		/** 
		 * @brief Obtain diffusion coefficient along the interacting particle direction. 
		 * @param[in] particle_index_i Index of particle i;
		 * @param[in] particle_index_j Index of particle j;
		 * @param[in] inter_particle_direction Normal direction pointing from i to j;
		 */
		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, size_t particle_index_j, Vecd& inter_particle_direction) override
		{
			Vecd grad_ij = transf_diffusivity_ * inter_particle_direction;
			return 1.0 / grad_ij.scalarNormSqr();
		};
	};

	/**
	 * @class DirectionBiasedDiffusionMaterial
	 * @brief Diffusision is biased along a specific direction.
	 */
	class LocalDirectionalDiffusion : public DirectionalDiffusion
	{
	protected:
		/* Local bias (usually due to fiber orientation) direction. */
		StdVec<Vecd> local_bias_direction_;
		/* Local transformed diffusivity with inverse Cholesky decomposition. */
		StdVec<Matd> local_transf_diffusivity_;
	public:
		/** Constructor*/
		LocalDirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
			Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
			: DirectionalDiffusion(diffusion_species_index, gradient_species_index, diff_cf, bias_diff_cf, bias_direction) {};
		virtual ~LocalDirectionalDiffusion() {};
		/**
		 * @brief Obtain diffusion coefficient along the interacting particle direction.
		 * @param[in] particle_index_i Index of particle i;
		 * @param[in] particle_index_j Index of particle j;
		 * @param[in] inter_particle_direction Normal direction pointing from i to j;
		 */
		virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, size_t particle_index_j, Vecd& inter_particle_direction) override
		{
			Matd trans_diffusivity = getAverageValue(local_transf_diffusivity_[particle_index_i], local_transf_diffusivity_[particle_index_j]);
			Vecd grad_ij = trans_diffusivity * inter_particle_direction;
			return 1.0 / grad_ij.scalarNormSqr();
		};
		/** initialize the local property. */
		virtual void initializeLocalDiffusionProperties(BaseParticles* base_particles) override;
		/** Setup the local property after initialization. */
		virtual void setupLocalDiffusionProperties(StdVec<Vecd>& material_fiber) override;
	};

	/** Reaction functor . */
	typedef std::function<Real(StdVec<Real>&)> ReactionFunctor;
	/**
	 * @class BaseReactionModel
	 * @brief Base class for all reaction models.
	 */
	class BaseReactionModel
	{
	protected:
		/** assign derived material properties*/
		virtual void assignDerivedReactionParameters() = 0;
	public:
		/** Default constructor. */
		BaseReactionModel() {};
		virtual ~BaseReactionModel() {};

		/** All reactive species. */
		IndexVector reactive_species_;
		/** Get the prodcution rates. */
		StdVec<ReactionFunctor> get_production_rates_;
		/** Get the loss rates. */
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

		virtual Real getProductionRateIonicCurrent(StdVec<Real>& species) = 0;
		virtual Real getLossRateIonicCurrent(StdVec<Real>& species) = 0;
		virtual Real getProductionRateGateVariable(StdVec<Real>& species) = 0;
		virtual Real getLossRateGateVariable(StdVec<Real>& species) = 0;
		virtual Real getProductionActiveContractionStress(StdVec<Real>& species);
		virtual Real getLossRateActiveContractionStress(StdVec<Real>& species);

		/** assign derived material properties*/
		virtual void assignDerivedReactionParameters() override {};
	public:
		ElectroPhysiologyReaction() : BaseReactionModel(), k_a_(1.0),
			voltage_(0), gate_variable_(1), active_contraction_stress_(2) {};
		virtual ~ElectroPhysiologyReaction() {};
		/** Initialize reaction model. */
		void initilaizeElectroPhysiologyReaction(size_t voltage,
			size_t gate_variable, size_t active_contraction_stress);
	};

	class AlievPanfilowModel : public ElectroPhysiologyReaction
	{
	protected:
		/** Parameters for two variable cell model. */
		Real k_, a_, b_, mu_1_, mu_2_, epsilon_, c_m_;

		virtual Real getProductionRateIonicCurrent(StdVec<Real>& species) override;
		virtual Real getLossRateIonicCurrent(StdVec<Real>& species) override;
		virtual Real getProductionRateGateVariable(StdVec<Real>& species) override;
		virtual Real getLossRateGateVariable(StdVec<Real>& species) override;

		/** assign derived material properties*/
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

		/** Total number of species. */
		size_t number_of_species_;
		/** Map from species names to indexes. */
		map<string, size_t> species_indexes_map_;
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
		/** The reaction model for all reactive species. */
		BaseReactionModel* species_reaction_;
		/** particles for this material */
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* diffusion_reaction_particles_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override 
		{
			BaseMaterialType::assignDerivedMaterialParameters();
		};
		/** 
		 * @brief Insert a species.
		 * @param[in] species_name Name of species
		 */
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

		/** Get number of species. */
		size_t getNumberOfSpecies() { return number_of_species_; };
		/** Get diffusion species */
		StdVec<BaseDiffusion*> getDiffusionSpecies() { return species_diffusion_; };
		/** Get reaction model */
		BaseReactionModel* getReactionModel() { return species_reaction_; };
		/** Get species to index map. */
		map<string, size_t> getSpeciesIndexMap() { return  species_indexes_map_; };
		/** assign particles to this material */
		void assignDiffusionReactionParticles(DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* diffusion_reaction_particles) {
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

		/** the interface for dynamical cast*/
		virtual DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>* PointToThisObject() override {
			return this;
		};
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

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override
		{
			DiffusionReactionMaterial<SolidParticles, Solid>::assignDerivedMaterialParameters();
		};
	public:
		/** Constructor*/
		MonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction);
		virtual ~MonoFieldElectroPhysiology() {};

		/** Initialize diffusion reaction material. */
		virtual void initializeDiffusion() override;

		/** the interface for dynamical cast*/
		virtual MonoFieldElectroPhysiology* PointToThisObject() override 
		{
			return this;
		};
	};

	/**
	 * @class LocalMonoFieldElectroPhysiology
	 * @brief material class for electro_physiology with locally oriented fibers.
	 */
	class LocalMonoFieldElectroPhysiology 
		: public MonoFieldElectroPhysiology
	{
	protected:
		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override
		{
			MonoFieldElectroPhysiology::assignDerivedMaterialParameters();
		};
	public:
		/** Constructor*/
		LocalMonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction)
			:MonoFieldElectroPhysiology(electro_physiology_reaction) {
			MonoFieldElectroPhysiology::material_name_ = "LocalMonoFieldElectroPhysiology";
		};
		virtual ~LocalMonoFieldElectroPhysiology() {};
		/** Initialize diffusion reaction material. */
		virtual void initializeDiffusion() override;
		/** Assign the fiber property into the diffusion material. */
		void assignFiberProperties(StdVec<Vecd> &material_fiber);
		/** the interface for dynamical cast*/
		virtual LocalMonoFieldElectroPhysiology* PointToThisObject() override 
		{
			return this;
		};
	};
}
