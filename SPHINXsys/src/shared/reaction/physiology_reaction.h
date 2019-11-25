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

namespace SPH 
{
	/**
	 * @class ElectroPhysiology
	 * @brief The Electrophysiology Reaction model
	 */
	class ElectroPhysiology : public Reaction
	{
	public:
		ElectroPhysiology(string reaction_model, SPHBody *body) : Reaction(reaction_model)
			{	
        		std::cout << "The electrophysiology reaction model is applied !" << std::endl;
				delete body->base_reaction_;
				body->base_reaction_ = this;
			}
		virtual ~ElectroPhysiology() {};
		/** the interface for dynamical cast*/
		virtual ElectroPhysiology* PointToThisObject() override { return this; };
		/**
		 * @brief Get the Production Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getProductionRateOfIonicCurrent(Real volatage, Real gate_variable) = 0;
		/**
		 * @brief Get the Lose Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variableo
		 */
		virtual Real getLoseRateOfIonicCurrent(Real volatage, Real gate_variable) = 0 ;
		/**
		 * @brief Get the Production Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable.
		 */
		virtual Real getProductionRateOfGateVarible(Real volatage, Real gate_variable) = 0;
		/**
		 * @brief Get the Lose Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getLoseRateOfGateVarible(Real volatage, Real gate_variable) = 0 ;
		/**
		 * @brief Get the production Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage
		 */
		virtual Real getProductionRateOfActiveContractionStress(Real volatage) = 0 ;
		/**
		 * @brief Get the lose Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage
		 */
		virtual Real getLoseRateOfActiveContractionStress(Real volatage) = 0 ;
	};
	/**
	 * @class AlievPanfilowModel
	 * @brief The simplest Electrophysiology Reaction model, which reduces the complex of array of ion currents to two variables that
	 * decribe excitation and recovery.
	 */
	class AlievPanfilowModel : public ElectroPhysiology
	{
	public:
		explicit AlievPanfilowModel(string reaction_model, SPHBody *body, Real c_m, Real k, Real a, 
			Real mu_1, Real mu_2, Real epsilon, Real k_a = 0.0);
		virtual ~AlievPanfilowModel() {};
		/** the interface for dynamical cast*/
		virtual AlievPanfilowModel* PointToThisObject() override { return this; };
		/**
		 * @brief Get the Production Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getProductionRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variableo
		 */
		virtual Real getLoseRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Production Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable.
		 */
		virtual Real getProductionRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getLoseRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the production Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage 
		 */
		virtual Real getProductionRateOfActiveContractionStress(Real volatage) override ;
		/**
		 * @brief Get the lose Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage
		 */
		virtual Real getLoseRateOfActiveContractionStress(Real volatage) override ;
		/** Parameters for two variable cell model. */
		Real k_, a_, mu_1_, mu_2_, epsilon_, c_m_,k_a_;
	};
	/**
	 * @class FitzHughNagumoModel
	 * @brief The simplest Electrophysiology Reaction model, which reduces the complex of array of ion currents to two variables that
	 * decribe excitation and recovery.
	 */
	class FitzHughNagumoModel : public ElectroPhysiology
	{
	public:
		explicit FitzHughNagumoModel(string reaction_model, SPHBody *body, Real c_m, Real a, Real epsilon, Real beta, 
			Real gamma, Real sigma, Real k_a);
		virtual ~FitzHughNagumoModel() {};
		/** the interface for dynamical cast*/
		virtual FitzHughNagumoModel* PointToThisObject() override { return this; };
		/**
		 * @brief Get the Production Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getProductionRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variableo
		 */
		virtual Real getLoseRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Production Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable.
		 */
		virtual Real getProductionRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getLoseRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the production Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage
		 */
		virtual Real getProductionRateOfActiveContractionStress(Real volatage) override ;
		/**
		 * @brief Get the lose Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage 
		 */
		virtual Real getLoseRateOfActiveContractionStress(Real volatage) override ;
		/** Parameters for two variable cell model. */
		Real c_m_, a_, epsilon_, beta_, gamma_, sigma_, k_a_;
	};
		/**
	 * @class AlievPanfilowModel
	 * @brief The modified AP model, https://doi.org/10.1016/j.cma.2017.03.015
	 */
	class ModifiedAlievPanfilowModel : public AlievPanfilowModel
	{
	public:
		explicit ModifiedAlievPanfilowModel(string reaction_model, SPHBody *body, Real c_m, Real k, Real a, 
			Real mu_1, Real mu_2, Real epsilon, Real k_a) : AlievPanfilowModel(reaction_model,body, c_m,k,a, 
			mu_1, mu_2, epsilon, k_a){};
		virtual ~ModifiedAlievPanfilowModel() {};
		/** the interface for dynamical cast*/
		virtual ModifiedAlievPanfilowModel* PointToThisObject() override { return this; };
		/**
		 * @brief Get the Production Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getProductionRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Ionic current.
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variableo
		 */
		virtual Real getLoseRateOfIonicCurrent(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Production Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable.
		 */
		virtual Real getProductionRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the Lose Rate of Gate variable
		 * @param[in] volatage Trans membrane volatage
		 * @param[in] gate_variable Gate variable
		 */
		virtual Real getLoseRateOfGateVarible(Real volatage, Real gate_variable) override;
		/**
		 * @brief Get the production Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage 
		 */
		virtual Real getProductionRateOfActiveContractionStress(Real volatage) override ;
		/**
		 * @brief Get the lose Rate of active contraction stress
		 * @param[in] volatage Trans membrane volatage
		 */
		virtual Real getLoseRateOfActiveContractionStress(Real volatage) override ;
	};
}