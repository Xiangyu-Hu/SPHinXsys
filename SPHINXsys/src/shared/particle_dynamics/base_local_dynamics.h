/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
* @file base_local_dynamics.h
* @brief This is for the base classes of local particle dynamics, which describe the
* dynamics of a particle.
* @author  Xiangyu Hu
*/

#ifndef BASE_LOCAL_DYNAMICS_H
#define BASE_LOCAL_DYNAMICS_H

#include "base_data_package.h"
#include "sph_data_containers.h"

namespace SPH
{

	/**
	 * @class BaseLocalDynamics
	 * @brief The new version of base class for all local particle dynamics.
	 */
	template <class ReturnType>
	class BaseLocalDynamics
	{
		SPHBody &sph_body_;

	public:
		explicit BaseLocalDynamics(SPHBody &sph_body) : sph_body_(sph_body){};
		virtual ~BaseLocalDynamics(){};

		typedef ReturnType LocalDynamicsReturnType;
		void setBodyUpdated() { sph_body_.setNewlyUpdated(); };
		virtual ReturnType setupDynamics(Real dt = 0.0) = 0; //setup global parameters
	};

	/**
	 * @class LocalDynamics
	 * @brief The new version of base class for all local particle dynamics,
	 * which loops along particles.
	 */
	class LocalDynamics : public BaseLocalDynamics<void>
	{
	public:
		explicit LocalDynamics(SPHBody &sph_body)
			: BaseLocalDynamics<void>(sph_body){};
		virtual ~LocalDynamics(){};

		/** the function for set global parameters for the particle dynamics */
		virtual void setupDynamics(Real dt = 0.0) override{};
	};

	/**
	 * @class LocalDynamicsReduce
	 * @brief The new version of base class for all local particle dynamics.
	 */
	template <typename ReturnType, typename ReduceOperation>
	class LocalDynamicsReduce : public LocalDynamics
	{
	public:
		explicit LocalDynamicsReduce(SPHBody &sph_body, ReturnType reference)
			: LocalDynamics(sph_body), reference_(reference),
			  quantity_name_("ReducedQuantity"){};
		virtual ~LocalDynamicsReduce(){};

		using ReduceReturnType = ReturnType;
		ReturnType InitialReference() { return reference_; };
		std::string QuantityName() { return quantity_name_; };
		ReduceOperation &getReduceOperation() { return operation_; };
		virtual ReturnType outputResult(ReturnType reduced_value) { return reduced_value; }

	protected:
		ReturnType reference_;
		ReduceOperation operation_;
		std::string quantity_name_;
	};    
}
#endif //BASE_LOCAL_DYNAMICS_H