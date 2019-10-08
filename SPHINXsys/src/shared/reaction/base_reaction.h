 /** 
 * @file 	base_reaction.h
 * @brief 	This is the class that presents the base reaction class.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Maybe chemical reaction in the future.
 * 			Chi Zhang
 */
#pragma once

#include "base_data_package.h"

using namespace std;

namespace SPH 
{
	/**
	 * @class Reaction
	 * @brief This class define the models of celluar electrophysiology.
	 */
	class Reaction
	{
	protected:
		const string reaction_model_;
	public:
		Reaction(string reaction_model) : reaction_model_(reaction_model)
		{

		};

		virtual ~Reaction() {};
		virtual Reaction* PointToThisObject() { return this; };
	};

}