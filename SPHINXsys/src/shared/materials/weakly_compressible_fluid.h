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
 * @file 	weakly_compressible_fluid.h
 * @brief 	Describe the weakly compressible fluid which is used 
 * 			model incompressible fluids. Here, we have included several equation of states.
 * 			Futhermore, A typical non-newtonian fluid model is included.  
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */

#pragma once

#include "base_material.h"

namespace SPH {

	class ViscoelasticFluidParticles;

	/**
	 * @class WeaklyCompressibleFluid
	 * @brief Linear equation of state (EOS).
	 */
	class WeaklyCompressibleFluid : public Fluid
	{
	protected:
		Real p0_; /**< reference pressure */

		virtual void assignDerivedMaterialParameters() override 
		{
			Fluid::assignDerivedMaterialParameters();
			p0_ = rho_0_ * c_0_ * c_0_;
		};
	public:
		explicit WeaklyCompressibleFluid() : Fluid(), p0_(1.0) {
			material_name_ = "WeaklyCompressibleFluid";
		};
		virtual ~WeaklyCompressibleFluid() {};

		virtual Real getPressure(Real rho) override;
		virtual Real DensityFromPressure(Real p) override;
		virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
		virtual WeaklyCompressibleFluid* ThisObjectPtr() override {return this;};
	};

	/**
	* @class WeaklyCompressibleFluidFreeSurface
	* @brief Equation of state (EOS) with cut-off pressure.
	*/
	template<class WeaklyCompressibleFluidType>
	class WeaklyCompressibleFluidFreeSurface : public WeaklyCompressibleFluid
	{
	protected:
		WeaklyCompressibleFluidType* fluid_;
		Real cutoff_pressure_, cutoff_density_;

		virtual void assignDerivedMaterialParameters() {
			WeaklyCompressibleFluid::assignDerivedMaterialParameters();
		};
	public:
		WeaklyCompressibleFluidFreeSurface(Real cutoff_pressure)
			: WeaklyCompressibleFluid(),
			cutoff_pressure_(cutoff_pressure) {
			fluid_ = new WeaklyCompressibleFluidType();
			material_name_ = fluid_->material_name_ + "FreeSurface";
			cutoff_density_ = fluid_->DensityFromPressure(cutoff_pressure);
		}; 
		virtual ~WeaklyCompressibleFluidFreeSurface() {};

		virtual Real getPressure(Real rho) override {
			return rho < cutoff_density_ ? cutoff_pressure_ : fluid_->getPressure(rho);
		};
	};

	/**
	 * @class SymmetricTaitFluid
	 * @brief Tait EOS for positive and negative pressure symmetrically.
	 */
	class SymmetricTaitFluid : public WeaklyCompressibleFluid
	{
	protected:
		
		int gamma_; /**< determine the stiffness of the fluid */

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override 
		{
			WeaklyCompressibleFluid::assignDerivedMaterialParameters();
		};

	public:
		SymmetricTaitFluid() : WeaklyCompressibleFluid(), gamma_(2) 
		{
			material_name_ = "SymmetricTaitFluid";
		};
		virtual ~SymmetricTaitFluid() {};

		virtual Real getPressure(Real rho) override;
		virtual Real DensityFromPressure(Real p) override;
		virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
	};

	/**
	 * @class Oldroyd_B_Fluid
	 * @brief linear EOS with relaxation time and polymetric viscosity.
	 */
	class Oldroyd_B_Fluid : public WeaklyCompressibleFluid
	{
	protected:
		Real lambda_; /**< relaxation time */
		Real mu_p_; /**< polymeric viscosity */
		ViscoelasticFluidParticles* viscoelastic_fluid_particles_;

		virtual void assignDerivedMaterialParameters() override 
		{
			WeaklyCompressibleFluid::assignDerivedMaterialParameters();
		};
	public:
		explicit Oldroyd_B_Fluid() : WeaklyCompressibleFluid(),
			lambda_(1.0), mu_p_(0.0) 
		{
			material_name_ = "Oldroyd_B_Fluid";
		};
		virtual ~Oldroyd_B_Fluid() {};

		void assignViscoelasticFluidParticles(ViscoelasticFluidParticles* viscoelastic_fluid_particles) 
		{
			viscoelastic_fluid_particles_ = viscoelastic_fluid_particles;
		};
		Real getReferenceRelaxationTime() { return lambda_; };
		Real ReferencePolymericViscosity() { return mu_p_; };
		virtual Oldroyd_B_Fluid* ThisObjectPtr() override {return this;};
	};
}
