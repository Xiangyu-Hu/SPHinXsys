/**
* @file 	diffusion_reaction_particles.h
* @brief 	This is the derived class of diffusion reaction particles.
* @author	Xiangyu Huand Chi Zhang
* @version	0.1
* @version 0.2.1
* add muscle particles and muscle data.
*/

#pragma once

#include "base_particles.h"
#include "base_material.h"
#include "base_body.h"
#include "diffusion_reaction.h"
#include "solid_particles.h"
#include "xml_engine.h"
#include <fstream>

using namespace std;
namespace SPH {

	/**
	 * @class DiffusionReactionData
	 * @brief Particles data for general diffusion/reaction problems.
	 */
	class DiffusionReactionData
	{
	public:
		/** Constructor. */
		DiffusionReactionData(size_t number_of_species);
		virtual ~DiffusionReactionData() {};

		/** The array of diffusion/reaction scalars. */
		StdVec<Real> species_n_;
		/** intermediate state for multi-step time integration */
		StdVec<Real> species_s_;
		/** The array of the time drivative. */
		StdVec<Real>  dspecies_dt_;
	};

	/**
	 * @class DiffusionReactionParticles
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template<class BaseParticlesType = BaseParticles, class BaseMaterialType = BaseMaterial>
	class DiffusionReactionParticles : public BaseParticlesType
	{
	protected:
		/** Total number of diffusion and reaction species . */
		size_t number_of_species_;
		map<string, size_t> species_indexes_map_;
	public:
		/** Vector of diffusion data. */
		StdLargeVec<DiffusionReactionData> diffusion_reaction_data_;

		/** Constructor. */
		DiffusionReactionParticles(SPHBody* body, 
			DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>* diffusion_reaction_material)
			: BaseParticlesType(body, diffusion_reaction_material)
		{
			diffusion_reaction_material->assignDiffusionReactionParticles(this);
			number_of_species_ 		= diffusion_reaction_material->getNumberOfSpecies();
			species_indexes_map_ 	= diffusion_reaction_material->getSpeciesIndexMap();
			for (size_t i = 0; i < this->base_particle_data_.size(); ++i)
			{
				diffusion_reaction_data_.push_back(DiffusionReactionData(number_of_species_));
			}
		};
		/** Destructor. */
		virtual ~DiffusionReactionParticles() {};

		/** Get species index map. */
		map<string, size_t> getSpeciesIndexMap() { return  species_indexes_map_; };
		/** add buffer particles which latter may be realized for particle dynamics*/
		virtual void AddABufferParticle() override {
			BaseParticlesType::AddABufferParticle();
			diffusion_reaction_data_.push_back(DiffusionReactionData(number_of_species_));
		};
		/** Copy state from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index,
			size_t another_particle_index) override 
		{
			BaseParticlesType::CopyFromAnotherParticle(this_particle_index, another_particle_index);
			diffusion_reaction_data_[this_particle_index] = diffusion_reaction_data_[another_particle_index];
		};
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index) override {
			BaseParticlesType::swapParticles(this_particle_index, that_particle_index);
			std::swap(diffusion_reaction_data_[this_particle_index], diffusion_reaction_data_[that_particle_index]);
		};
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream& output_file) override {
			BaseParticlesType::WriteParticlesToVtuFile(output_file);

			map<string, size_t>::iterator itr;
			size_t number_of_particles = this->body_->number_of_particles_;
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) {
				output_file << "    <DataArray Name=\" "<< itr->first <<" \" type=\"Float32\" Format=\"ascii\">\n";
				output_file << "    ";
				size_t k = itr->second;
				for (size_t i = 0; i != number_of_particles; ++i) {
					output_file << diffusion_reaction_data_[i].species_n_[k] << " ";
				}
				output_file << std::endl;
				output_file << "    </DataArray>\n";
			}
		};
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream& output_file) override
		{
			//BaseParticlesType::WriteParticlesToPltFile(output_file);
			size_t number_of_particles = this->body_->number_of_particles_;
			map<string, size_t>::iterator itr;
			
			if(Vecd(0).size()==2)
			{
				output_file << " VARIABLES = \" x \", \"y\", \"ID \" " << " ";
			}
			if(Vecd(0).size()==3)
			{
				output_file << " VARIABLES = \" x \", \"y\" , \"z\", \" ID \" " << " ";
			}

			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) 
			{
				output_file << ","<< "\"" << itr->first << " \" " << "";
			}
			output_file << "\n";

			for (size_t i = 0; i != number_of_particles; ++i)
			{
				for(size_t j = 0; j != Vecd(0).size(); ++j)
				{
					output_file << this->body_->base_particles_->base_particle_data_[i].pos_n_[j] << "  ";
				}

				output_file << i << " ";

				for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) 
				{
					size_t k = itr->second;
					output_file << diffusion_reaction_data_[i].species_n_[k] << " ";
				}

				output_file << "\n";
			}
		};
		/** Write particle data in XML format. */
		virtual void WriteParticlesToXmlFile(std::string& filefullpath) override{};
		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string& filefullpath) override {};
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string& filefullpath) override {};
		/** Pointer to this object. */
		virtual DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* 
			PointToThisObject() override { return this; };
	};

	/**
	 * @class ElectroPhysiologyParticles
	 * @brief A group of particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyParticles 
		: public DiffusionReactionParticles<SolidParticles, Solid>
	{
	public:
		/** Constructor. */
		ElectroPhysiologyParticles(SPHBody* body, 
			DiffusionReactionMaterial<SolidParticles, Solid>* diffusion_reaction_material);
		/** Destructor. */
		virtual ~ElectroPhysiologyParticles() {};

		/** Pointer to this object. */
		virtual ElectroPhysiologyParticles* PointToThisObject() override { return this; };
	};
}