/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	inelastic_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These classes are based on isotropic linear elastic solid.
* 			Several more complex materials, including neo-hookean, FENE noe-hookean
*			and anisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
*/
#pragma once

#include "elastic_solid.h"

namespace SPH {

	/**
	* @class PlasticSolid
	* @brief Abstract class for a generalized plastic solid
	*/
	class PlasticSolid : public NeoHookeanSolid
	{
	protected:
		Real yield_stress_;

		virtual void initializePlasticParameters() = 0;
	public:
		/** Constructor */
		PlasticSolid() :NeoHookeanSolid()
		{
			material_name_ = "PlasticSolid";
		};
		virtual ~PlasticSolid() {};

		Real YieldStress() { return yield_stress_; };
		/** compute the stress through defoemation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd& deformation, size_t index_i, Real dt = 0.0) = 0;

		virtual PlasticSolid* ThisObjectPtr() override { return this; };
	};

	/**
	 * @class HardeningPlasticSolid
	 * @brief Class for a generalized plastic solid
	 */
	class HardeningPlasticSolid : public PlasticSolid
	{
	protected:
		Real hardening_modulus_;
		const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;
		const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
		StdLargeVec<Matd>	inverse_plastic_strain_;	/**< inverse of plastic right cauchy green strain tensor */
		StdLargeVec<Real>	hardening_parameter_;		/**< hardening parameter */

		virtual void initializePlasticParameters() override;
	public:
		/** Constructor */
		HardeningPlasticSolid() :PlasticSolid()
		{
			material_name_ = "HardeningPlasticSolid";
		};
		virtual ~HardeningPlasticSolid() {};

		Real HardeningModulus() { return hardening_modulus_; };
		/** assign particles to this material */
		virtual void assignElasticSolidParticles(ElasticSolidParticles* elastic_particles) override;
		/** compute the stress through defoemation, and plastic relaxation. */
		virtual Matd PlasticConstitutiveRelation(const Matd& deformation, size_t index_i, Real dt = 0.0) override;

		virtual HardeningPlasticSolid* ThisObjectPtr() override { return this; };
	};
}
