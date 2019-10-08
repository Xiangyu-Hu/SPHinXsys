 /** 
 * @file 	physiology_reaction.h
 * @brief 	This is the class that presents the two variable model of celluar electrophysiology.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
#pragma once
#include "base_reaction.h"
#include "base_body.h"

using namespace std;

namespace SPH {
	/**
	 * @class ElectrophysiologyReaction
	 * @brief The simplest Electrophysiology Reaction model, which reduces the complex of array of ion currents to two variables that
	 * decribe excitation and recovery.
	 */
	class ElectrophysiologyReaction : public Reaction
	{
	public:
		ElectrophysiologyReaction(string reaction_model, SPHBody *body, Real c_m, Real k, Real a, 
			Real mu_1, Real mu_2, Real epsilon) : Reaction(reaction_model), 
			c_m_(c_m), k_(k), a_(a), mu_1_(mu_1), mu_2_(mu_2), epsilon_(epsilon)
			{	
        		std::cout << "The electrophysiology reaction model is applied !" << std::endl;
				delete body->base_reaction_;
				body->base_reaction_ = this;
			}
		virtual ~ElectrophysiologyReaction() {};
		/** the interface for dynamical cast*/
		virtual ElectrophysiologyReaction* PointToThisObject() override { return this; };
		/**
		 * @brief Get the Production Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getProductionRateOfIonicCurrent(Real volatage, Real gate_variable);
		/**
		 * @brief Get the Lose Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variableo
		 */
		virtual Real getLoseRateOfIonicCurrent(Real volatage, Real gate_variable);
		/**
		 * @brief Get the Production Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable.
		 */
		virtual Real getProductionRateOfGateVarible(Real volatage, Real gate_variable);
		/**
		 * @brief Get the Lose Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getLoseRateOfGateVarible(Real volatage, Real gate_variable);
		/** Parameters for two variable cell model. */
		Real k_, a_, mu_1_, mu_2_, epsilon_, c_m_;
	};
}