/**
 * @file diffusion_reaction.cpp
 * @brief These are classes for diffusion and reaction properties
 * @author Chi Zhang and Xiangyu Hu
 */

#include "diffusion_reaction.h"

#include "diffusion_reaction_particles.h"

namespace SPH 
{
//=================================================================================================//
	void DirectionalDiffusion::initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction) 
	{
		bias_diff_cf_ 		= bias_diff_cf;
		bias_direction_ 	= bias_direction;
		Matd diff_i 		= diff_cf_* Matd(1.0)
							+ bias_diff_cf_ * SimTK::outer(bias_direction_, bias_direction_);
		transformed_diffusivity_ = inverseCholeskyDecomposition(diff_i);
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::assignBaseParticles(BaseParticles* base_particles)
	{
		DirectionalDiffusion::assignBaseParticles(base_particles);
		initializeFiberDirection();
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::initializeFiberDirection()
	{
		base_particles_->registerAVariable<indexVector, Vecd>(local_bias_direction_, "Fiber");
		base_particles_->addAVariableNameToList<indexVector, Vecd>(reload_local_parameters_, "Fiber");
	}
	//=================================================================================================//
	void LocalDirectionalDiffusion::readFromXmlForLocalParameters(std::string& filefullpath)
	{
		BaseMaterial::readFromXmlForLocalParameters(filefullpath);
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t i = 0; i != total_real_particles; i++)
		{
			Matd diff_i = diff_cf_ * Matd(1.0) 
				+ bias_diff_cf_ * SimTK::outer(local_bias_direction_[i], local_bias_direction_[i]);
			local_transformed_diffusivity_.push_back(inverseCholeskyDecomposition(diff_i));
		}
		std::cout << "\n Local diffusion parameters setup finished " << std::endl;
	};
	//=================================================================================================//
	void ElectroPhysiologyReaction::initializeElectroPhysiologyReaction(size_t voltage, size_t gate_variable, 
			size_t active_contraction_stress)
	{
		voltage_ 					= voltage;
		gate_variable_ 				= gate_variable;
		active_contraction_stress_ 	= active_contraction_stress;

		reactive_species_.push_back(voltage);
		reactive_species_.push_back(gate_variable);		
		reactive_species_.push_back(active_contraction_stress);

		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateIonicCurrent, this, _1, _2));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateGateVariable, this, _1, _2));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionActiveContractionStress, this, _1, _2));

		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateIonicCurrent, this, _1, _2));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateGateVariable, this, _1, _2));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateActiveContractionStress, this, _1, _2));
	};
//=================================================================================================//
	Real ElectroPhysiologyReaction::
		getProductionActiveContractionStress(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real voltage_dim = species[voltage_][particle_i] * 100.0 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
	}
//=================================================================================================//
	Real ElectroPhysiologyReaction::
		getLossRateActiveContractionStress(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real voltage_dim = species[voltage_][particle_i] * 100.0 - 80.0;
		return 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
	}
//=================================================================================================//
	Real AlievPanfilowModel::
		getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real voltage = species[voltage_][particle_i];
		return - k_ * voltage * (voltage * voltage - a_ * voltage - voltage) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::
		getLossRateIonicCurrent(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real gate_variable = species[gate_variable_][particle_i];
		return  (k_ * a_ + gate_variable) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::
		getProductionRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real voltage = species[voltage_][particle_i];
		Real gate_variable = species[gate_variable_][particle_i];
		Real temp = epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
		return - temp * k_ * voltage * (voltage - b_ - 1.0);
	}
//=================================================================================================//
	Real AlievPanfilowModel::
		getLossRateGateVariable(StdVec<StdLargeVec<Real>>& species, size_t particle_i)
	{
		Real voltage = species[voltage_][particle_i];
		Real gate_variable = species[gate_variable_][particle_i];
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
	}
//=================================================================================================//
	MonoFieldElectroPhysiology::MonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction)
		: DiffusionReactionMaterial<SolidParticles, Solid>(electro_physiology_reaction), diff_cf_(1.0), bias_diff_cf_(0.0), 
		bias_direction_(FirstAxisVector(Vecd(0)))
	{
		material_name_ = "MonoFieldElectroPhysiology";
		insertASpecies("Voltage");
		insertASpecies("GateVariable");
		insertASpecies("ActiveContractionStress");

		electro_physiology_reaction->initializeElectroPhysiologyReaction(species_indexes_map_["Voltage"],
			species_indexes_map_["GateVariable"], species_indexes_map_["ActiveContractionStress"]);
	};
//=================================================================================================//
	void MonoFieldElectroPhysiology::initializeDiffusion()
	{
		DirectionalDiffusion* voltage_diffusion
			= new DirectionalDiffusion(species_indexes_map_["Voltage"], species_indexes_map_["Voltage"],
				diff_cf_, bias_diff_cf_, bias_direction_);
		species_diffusion_.push_back(voltage_diffusion);
	}
	//=================================================================================================//
	void LocalMonoFieldElectroPhysiology::initializeDiffusion()
	{
		LocalDirectionalDiffusion* voltage_diffusion
			= new LocalDirectionalDiffusion(species_indexes_map_["Voltage"], species_indexes_map_["Voltage"],
				diff_cf_, bias_diff_cf_, bias_direction_);
		species_diffusion_.push_back(voltage_diffusion);
	}
//=================================================================================================//
	void LocalMonoFieldElectroPhysiology::readFromXmlForLocalParameters(std::string& filefullpath)
	{
		species_diffusion_[0]->readFromXmlForLocalParameters(filefullpath);
	}
//=================================================================================================//
}
//=================================================================================================//
