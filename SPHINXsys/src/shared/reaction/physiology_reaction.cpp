/** 
 * @file 	physiology_reaction.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */
#include "physiology_reaction.h"
//=================================================================================================//
namespace SPH 
{
//=================================================================================================//
	AlievPanfilowModel::AlievPanfilowModel(string reaction_model, SPHBody *body, Real c_m, Real k, Real a, 
			Real mu_1, Real mu_2, Real epsilon, Real k_a) : ElectroPhysiology(reaction_model, body), 
			c_m_(c_m), k_(k), a_(a), mu_1_(mu_1), mu_2_(mu_2), epsilon_(epsilon), k_a_(k_a)
	{	
		//
	}
//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return -k_ * volatage * (volatage * volatage - a_ * volatage - volatage) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::getLoseRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return  (k_ * a_ + gate_variable) / c_m_;
	}
//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateOfGateVarible(Real volatage, Real gate_variable)
	{
		Real epsilon = epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
		return -epsilon * k_ * volatage * (volatage - a_ - 1.0);
	}
//=================================================================================================//
	Real AlievPanfilowModel::getLoseRateOfGateVarible(Real volatage, Real gate_variable)
	{
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
	}
//=================================================================================================//
	Real AlievPanfilowModel::getProductionRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
		// Real factor = volatage >= 0.05 ? 1.0 : 10.0;
		// return factor * k_a_ * volatage;
	}
//=================================================================================================//
	Real AlievPanfilowModel::getLoseRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		return 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		// return volatage >= 0.05 ? 1.0 : 10.0;
	}
//=================================================================================================//
	FitzHughNagumoModel::FitzHughNagumoModel(string reaction_model, SPHBody *body, Real c_m, Real a, Real epsilon, Real beta, 
			Real gamma, Real sigma, Real k_a) : ElectroPhysiology(reaction_model, body), 
			c_m_(c_m), a_(a), epsilon_(epsilon), beta_(beta), gamma_(gamma), sigma_(sigma), k_a_(k_a)
	{	

	}
//=================================================================================================//
	Real FitzHughNagumoModel::getProductionRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return volatage * (-volatage * volatage + a_ * volatage + volatage) / c_m_ - gate_variable / c_m_ ;
	}
//=================================================================================================//
	Real FitzHughNagumoModel::getLoseRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return  a_ / c_m_;
	}
//=================================================================================================//
	Real FitzHughNagumoModel::getProductionRateOfGateVarible(Real volatage, Real gate_variable)
	{
		return epsilon_ * (beta_ * volatage - sigma_);
	}
//=================================================================================================//
	Real FitzHughNagumoModel::getLoseRateOfGateVarible(Real volatage, Real gate_variable)
	{
		return epsilon_ * gamma_ ;
	}
//=================================================================================================//
	Real FitzHughNagumoModel::getProductionRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
	}
//=================================================================================================//
	Real FitzHughNagumoModel::getLoseRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor;
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getProductionRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return -52 * volatage * (volatage * volatage - 0.05 * volatage - volatage) / c_m_;
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getLoseRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return  (52 * 0.05 + 8 * gate_variable) / c_m_;
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getProductionRateOfGateVarible(Real volatage, Real gate_variable)
	{
		Real epsilon = epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
		return -epsilon * 8 * volatage * (volatage - 0.35 - 1.0);
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getLoseRateOfGateVarible(Real volatage, Real gate_variable)
	{
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getProductionRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		Real factor = 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
		return factor * k_a_ * (voltage_dim + 80.0);
	}
//=================================================================================================//
	Real ModifiedAlievPanfilowModel::getLoseRateOfActiveContractionStress(Real volatage)
	{
		Real voltage_dim = volatage * 100 - 80.0;
		return 0.1 + (1.0 - 0.1) * exp(-exp(-voltage_dim));
	}
//=================================================================================================//
}
//=================================================================================================//