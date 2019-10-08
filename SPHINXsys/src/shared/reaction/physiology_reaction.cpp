/** 
 * @file 	physiology_reaction.cpp
 * @brief 	These are functions defined in physiology_reaction.h
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
#include "physiology_reaction.h"
//=================================================================================================//
namespace SPH 
{
//=================================================================================================//
	Real ElectrophysiologyReaction::getProductionRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return -k_ * volatage * (volatage * volatage - a_ * volatage - volatage) / c_m_;
	}
//=================================================================================================//
	Real ElectrophysiologyReaction::getLoseRateOfIonicCurrent(Real volatage, Real gate_variable)
	{
		return  (k_ * a_ + gate_variable) / c_m_;
	}
//=================================================================================================//
	Real ElectrophysiologyReaction::getProductionRateOfGateVarible(Real volatage, Real gate_variable)
	{
		Real epsilon = epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
		return -epsilon * k_ * volatage * (volatage - a_ - 1.0);
	}
//=================================================================================================//
	Real ElectrophysiologyReaction::getLoseRateOfGateVarible(Real volatage, Real gate_variable)
	{
		return epsilon_ + mu_1_ * gate_variable / (mu_2_ + volatage + 1.0e-6);
	}
//=================================================================================================//
}
//=================================================================================================//