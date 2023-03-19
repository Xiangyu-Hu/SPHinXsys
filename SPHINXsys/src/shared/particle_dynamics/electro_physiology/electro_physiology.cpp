#include "electro_physiology.h"

namespace SPH
{
	//=================================================================================================//
	void ElectroPhysiologyReaction::initializeElectroPhysiologyReaction()
	{
		reactive_species_.push_back(voltage_);
		reactive_species_.push_back(gate_variable_);
		reactive_species_.push_back(active_contraction_stress_);

		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateIonicCurrent, this, _1));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionRateGateVariable, this, _1));
		get_production_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getProductionActiveContractionStress, this, _1));

		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateIonicCurrent, this, _1));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateGateVariable, this, _1));
		get_loss_rates_.push_back(std::bind(&ElectroPhysiologyReaction::getLossRateActiveContractionStress, this, _1));
	};
	//=================================================================================================//
	Real ElectroPhysiologyReaction::getProductionActiveContractionStress(LocalSpecies &species)
	{
		Real voltage_dim = species[voltage_] * 100.0 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
	}
	//=================================================================================================//
	Real ElectroPhysiologyReaction::getLossRateActiveContractionStress(LocalSpecies &species)
	{
		Real voltage_dim = species[voltage_] * 100.0 - 80.0;
		return 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
	}
	//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateIonicCurrent(LocalSpecies &species)
	{
		Real voltage = species[voltage_];
		return -k_ * voltage * (voltage * voltage - a_ * voltage - voltage) / c_m_;
	}
	//=================================================================================================//
	Real AlievPanfilowModel::getLossRateIonicCurrent(LocalSpecies &species)
	{
		Real gate_variable = species[gate_variable_];
		return (k_ * a_ + gate_variable) / c_m_;
	}
	//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateGateVariable(LocalSpecies &species)
	{
		Real voltage = species[voltage_];
		Real gate_variable = species[gate_variable_];
		Real temp = epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
		return -temp * k_ * voltage * (voltage - b_ - 1.0);
	}
	//=================================================================================================//
	Real AlievPanfilowModel::
		getLossRateGateVariable(LocalSpecies &species)
	{
		Real voltage = species[voltage_];
		Real gate_variable = species[gate_variable_];
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + voltage + Eps);
	}
	//=================================================================================================//
	MonoFieldElectroPhysiology::
		MonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
								   Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
		: DiffusionReaction<Solid, 3>(electro_physiology_reaction, electro_physiology_reaction.getSpeciesNameList())
	{
		material_type_name_ = "MonoFieldElectroPhysiology";
		initializeAnDiffusion<DirectionalDiffusion>("Voltage", "Voltage", diff_cf, bias_diff_cf, bias_direction);
	};
	//=================================================================================================//
	LocalMonoFieldElectroPhysiology::
		LocalMonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
										Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
		: DiffusionReaction<Solid, 3>(electro_physiology_reaction, electro_physiology_reaction.getSpeciesNameList())
	{
		material_type_name_ = "LocalMonoFieldElectroPhysiology";
		initializeAnDiffusion<LocalDirectionalDiffusion>("Voltage", "Voltage", diff_cf, bias_diff_cf, bias_direction);
	}
	//=================================================================================================//
	ElectroPhysiologyParticles::ElectroPhysiologyParticles(
		SPHBody &sph_body, DiffusionReaction<Solid, 3> *diffusion_reaction_material)
		: DiffusionReactionParticles<SolidParticles, Solid, 3>(sph_body, diffusion_reaction_material) {}
	//=================================================================================================//
	ElectroPhysiologyReducedParticles::ElectroPhysiologyReducedParticles(
		SPHBody &sph_body, DiffusionReaction<Solid, 3> *diffusion_reaction_material)
		: DiffusionReactionParticles<SolidParticles, Solid, 3>(sph_body, diffusion_reaction_material) {}
	//=================================================================================================//
}
