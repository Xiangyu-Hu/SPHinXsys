/**
 * @file 	diffusion_reaction_particles_with_boundary.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H
#define DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H

#include "base_particles.h"
#include "base_body.h"
#include "base_material.h"
#include "diffusion_reaction.h"
//#include "diffusion_reaction_particles.h"

namespace SPH
{

	/**
	 * @class DiffusionReactionParticlesWithBoundary
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template <class BaseParticlesType, class DiffusionReactionMaterialType>
	class DiffusionReactionParticlesWithBoundary : public BaseParticlesType
	{
	protected:
		StdVec<std::string> all_species_names_;
		std::map<std::string, size_t> species_indexes_map_;
		StdVec<StdLargeVec<Real>*> diffusion_species_;
		StdVec<StdLargeVec<Real>*> gradient_species_;
		StdVec<StdLargeVec<Real>*> reactive_species_;

	public:
		StdVec<StdLargeVec<Real>> all_species_; /**< array of diffusion/reaction scalars */
		DiffusionReactionMaterialType& diffusion_reaction_material_;
		static constexpr int NumReactiveSpecies = DiffusionReactionMaterialType::NumReactiveSpecies;
		typedef DiffusionReactionMaterialType DiffusionReactionMaterial;

		StdLargeVec<Vecd> normal_vector_; /**< normal vector for Neumann or Robin boundary conditions */
		StdLargeVec<Real> heat_flux_; /**< heat flux for Neumann boundary condition */
		StdLargeVec<Real> convection_; /**< convection for Robin boundary condition */
		StdLargeVec<Real> T_infinity_; /**< T_infinity for Robin boundary condition */

		DiffusionReactionParticlesWithBoundary(SPHBody& sph_body, DiffusionReactionMaterialType* diffusion_reaction_material)
			: BaseParticlesType(sph_body, diffusion_reaction_material),
			all_species_names_(diffusion_reaction_material->AllSpeciesNames()),
			species_indexes_map_(diffusion_reaction_material->AllSpeciesIndexMap()),
			diffusion_reaction_material_(*diffusion_reaction_material)
		{
			all_species_.resize(all_species_names_.size());
			const IndexVector& diffusion_species_indexes = diffusion_reaction_material_.DiffusionSpeciesIndexes();
			const IndexVector& gradient_species_indexes = diffusion_reaction_material_.GradientSpeciesIndexes();
			for (size_t i = 0; i != diffusion_species_indexes.size(); ++i)
			{
				diffusion_species_.push_back(&all_species_[diffusion_species_indexes[i]]);
				gradient_species_.push_back(&all_species_[gradient_species_indexes[i]]);
			}

			const IndexVector& reactive_species_indexes = diffusion_reaction_material_.ReactiveSpeciesIndexes();
			for (size_t i = 0; i != reactive_species_indexes.size(); ++i)
			{
				reactive_species_.push_back(&all_species_[reactive_species_indexes[i]]);
			}

			normal_vector_.resize(this->all_species_.size());
			heat_flux_.resize(this->all_species_.size());
			convection_.resize(this->all_species_.size());
			T_infinity_.resize(this->all_species_.size());
		};
		virtual ~DiffusionReactionParticlesWithBoundary() {};

		StdVec<std::string>& AllSpeciesNames() { return all_species_names_; };
		std::map<std::string, size_t> AllSpeciesIndexMap() { return species_indexes_map_; };
		StdVec<StdLargeVec<Real>*>& DiffusionSpecies() { return diffusion_species_; };
		StdVec<StdLargeVec<Real>*>& GradientSpecies() { return gradient_species_; };
		StdVec<StdLargeVec<Real>*>& ReactiveSpecies() { return reactive_species_; };

		virtual void initializeOtherVariables() override
		{
			BaseParticlesType::initializeOtherVariables();

			std::map<std::string, size_t>::iterator itr;
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr)
			{
				// Register a specie.
				this->registerVariable(all_species_[itr->second], itr->first);
				/** the scalars will be sorted if particle sorting is called, Note that we call a template function from a template class. */
				this->template registerSortableVariable<Real>(itr->first);
				/** add species to basic output particle data. */
				this->template addVariableToWrite<Real>(itr->first);
			}

			this->registerVariable(normal_vector_, "UnitNormalVector");
			this->template addVariableToWrite<Vecd>("UnitNormalVector");
			
			this->registerVariable(heat_flux_, "HeatFlux");
			this->template addVariableToWrite<Real>("HeatFlux");
			
			this->registerVariable(convection_, "Convection");
			this->template addVariableToWrite<Real>("Convection");

			this->registerVariable(T_infinity_, "T_infinity");
			this->template addVariableToWrite<Real>("T_infinity");

		};

		virtual DiffusionReactionParticlesWithBoundary<BaseParticlesType, DiffusionReactionMaterialType>* ThisObjectPtr() override { return this; };
	};
}

#endif //  DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H