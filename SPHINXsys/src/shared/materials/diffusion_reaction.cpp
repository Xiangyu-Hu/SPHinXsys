/**
 * @file diffusion_reaction.cpp
 * @brief These are classes for diffusion and reaction properties
 * @author Chi Zhang and Xiangyu Hu
 * @version 0.1.0
 */

#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"

namespace SPH 
{
//=================================================================================================//
	void DirectionalDiffusion::intializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction) 
	{
		bias_diff_cf_ 		= bias_diff_cf;
		bias_direction_ 	= bias_direction;
		Matd diff_i 		= diff_cf_* Matd(1.0)
							+ bias_diff_cf_ * SimTK::outer(bias_direction_, bias_direction_);
		transf_diffusivity_ = inverseCholeskyDecomposition(diff_i);
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::initializeLocalProperties(BaseParticles* base_particles)
	{
		size_t number_of_particles = base_particles->getSPHBody()->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; i++)
		{
			local_bias_direction_.push_back(Vecd(0));
			local_transf_diffusivity_.push_back(Matd(0));
		}
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::setupLocalProperties(StdVec<Vecd>& material_fiber)
	{
		if (material_fiber.size() != local_bias_direction_.size()) {
			std::cout << "\n Error:  material properties does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		for (size_t i = 0; i < material_fiber.size(); i++)
		{
			local_bias_direction_[i] = material_fiber[i];
			Matd diff_i = diff_cf_ * Matd(1.0) + bias_diff_cf_ * SimTK::outer(material_fiber[i], material_fiber[i]);
			local_transf_diffusivity_[i] = inverseCholeskyDecomposition(diff_i);
		}
		std::cout << "\n Local diffusion properties setup finished " << std::endl;
	};
	//=================================================================================================//
	void ElectroPhysiologyReaction::initilaizeElectroPhysiologyReaction(size_t voltage, size_t gate_variable, 
			size_t active_contraction_stress)
	{
		voltage_ 					= voltage;
		gate_variable_ 				= gate_variable;
		active_contraction_stress_ 	= active_contraction_stress;

		reactive_species_.push_back(voltage);
		reactive_species_.push_back(gate_variable);		
		reactive_species_.push_back(active_contraction_stress);

		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateIonicCurrent, this, _1));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateGateVariable, this, _1));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionActiveContractionStress, this, _1));

		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateIonicCurrent, this, _1));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateGateVariable, this, _1));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateActiveContractionStress, this, _1));
	};
//=================================================================================================//
	Real ElectroPhysiologyReaction::getProductionActiveContractionStress(StdVec<Real>& species)
	{
		Real voltage_dim = species[voltage_] * 100.0 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
	}
//=================================================================================================//
	Real ElectroPhysiologyReaction::getLossRateActiveContractionStress(StdVec<Real>& species)
	{
		Real voltage_dim = species[voltage_] * 100.0 - 80.0;
		return 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
	}
//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateIonicCurrent(StdVec<Real>& species)
	{
		Real voltage = species[voltage_];
		return - k_ * voltage * (voltage * voltage - a_ * voltage - voltage) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::getLossRateIonicCurrent(StdVec<Real>& species)
	{
		Real gate_variable = species[gate_variable_];
		return  (k_ * a_ + gate_variable) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateGateVariable(StdVec<Real>& species)
	{
		Real voltage = species[voltage_];
		Real gate_variable = species[gate_variable_];
		Real temp = epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + 1.0e-6);
		return - temp * k_ * voltage * (voltage - b_ - 1.0);
	}
//=================================================================================================//
	Real AlievPanfilowModel::getLossRateGateVariable(StdVec<Real>& species)
	{
		Real voltage = species[voltage_];
		Real gate_variable = species[gate_variable_];
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + 1.0e-6);
	}
//=================================================================================================//
	MonoFieldElectroPhysiology::MonoFieldElectroPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction)
		: DiffusionReactionMaterial<SolidParticles, Solid>(electro_physiology_reaction), diff_cf_(1.0), bias_diff_cf_(0.0), 
		bias_direction_(FisrtAxisVector(Vecd(0)))
	{
		material_name_ = "MonoFieldElectroPhysiology";
		insertASpecies("Voltage");
		insertASpecies("GateVariable");
		insertASpecies("ActiveContractionStress");

		electro_physiology_reaction->initilaizeElectroPhysiologyReaction(species_indexes_map_["Voltage"],
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
	void LocalMonoFieldElectroPhysiology::initializeLocalProperties(BaseParticles* base_particles)
	{
		species_diffusion_[0]->initializeLocalProperties(base_particles);
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
	void LocalMonoFieldElectroPhysiology::assignFiberProperties(StdVec<Vecd> &material_fiber)
	{
		species_diffusion_[0]->setupLocalProperties(material_fiber);
	}
//=================================================================================================//
}
//=================================================================================================//
