 /** 
 * @file 	chemical_reaction.h
 * @brief 	This is the class that presents chemical reaction model
 * @author 	Chi Zhang Jianhang Wang and Xiangyu Hu
 * @version 0.1.0
 * 			The starting point of a chemical reaction model.
 * 			Chi Zhang
 */
#pragma once
#include "base_reaction.h"
#include "base_body.h"

using namespace std;

namespace SPH 
{
	/**
	 * @class ChemicalSpecies
	 * @brief Species of chemical reaction.
	 */
	// class ChemicalSpecies
	// {

	// }
	/**
	 * @class ChemicalReaction
	 * @brief The Chemical Reaction model
	 */
	class ChemicalReaction : public Reaction
	{
	public:
		ChemicalReaction(string reaction_model, SPHBody *body) : Reaction(reaction_model)
		{	
        		std::cout << "The chemical reaction model is applied !" << std::endl;
				delete body->base_reaction_;
				body->base_reaction_ = this;
		}
		virtual ~ChemicalReaction() {};
		/** the interface for dynamical cast*/
		virtual ChemicalReaction* PointToThisObject() override { return this; };
    };
}