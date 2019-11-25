/**
 * @file 	weakly_compressible_fluid.h
 * @brief 	Desrcibe the weakly compressible fluid which is used 
 * 			model incompressible fluids. Here, we have included serveral equation of states.
 * 			Futhermore, A typical non-newtonian fluid model is included.  
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 * @version 0.1.0
 */

#pragma once

#include "base_material.h"

namespace SPH {

	/**
	 * @class WeaklyCompressibleFluid
	 * @brief Linear equation of state (EOS).
	 */
	class WeaklyCompressibleFluid : public Fluid
	{
	protected:
		/** reference pressure */
		Real p0_;
	public:
		/** constructor */
		explicit WeaklyCompressibleFluid(string fluid_name, SPHBody *body, Real rho_0 = 1.0,
			Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0);
		/** Constructor without body information. */
		explicit WeaklyCompressibleFluid(Real rho_0 = 1.0, Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0);
		virtual ~WeaklyCompressibleFluid() {};

		/** the interface for dynamical cast*/
		virtual WeaklyCompressibleFluid* PointToThisObject() override { return this; };

		virtual Real GetPressure(Real rho) override;
		virtual Real ReinitializeRho(Real p) override;
		virtual Real GetSoundSpeed(Real p = 0.0, Real rho = 1.0) override;

		/** riemann soslver */
		virtual Real RiemannSolverForPressure(Real rhol, Real Rhor, Real pl, 
			Real pr, Real ul, Real ur) override;
		virtual Real RiemannSolverForVelocity(Real rhol, Real Rhor, Real pl, 
			Real pr, Real ul, Real ur) override;
	};

	/**
	* @class FluidFreeSurface
	* @brief Linear equation of state (EOS) with cut-off pressure.
	*/
	template<class WeaklyCompressibleFluidType>
	class WeaklyCompressibleFluidFreeSurface : public WeaklyCompressibleFluid
	{
	protected:
		WeaklyCompressibleFluidType *fluid_;
		Real cutoff_pressure_, cutoff_density_;
	public:
		/** constructor */
		explicit WeaklyCompressibleFluidFreeSurface(string fluid_name, SPHBody *body, 
			Real cutoff_pressure = 0.0, Real rho_0 = 1.0, Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0) 
			: WeaklyCompressibleFluid(fluid_name, body, rho_0, c_0, mu, k), 
			cutoff_pressure_(cutoff_pressure) {
			fluid_ = new WeaklyCompressibleFluidType(rho_0, c_0, mu, k);
				cutoff_density_ = fluid_->ReinitializeRho(cutoff_pressure);
		};
		virtual ~WeaklyCompressibleFluidFreeSurface() {};

		virtual Real GetPressure(Real rho) override { 
			return rho < cutoff_density_ ? cutoff_pressure_ : fluid_->GetPressure(rho);
		};
	};

	/**
	 * @class SymmetricTaitFluid
	 * @brief linear EOS for negative presssure and Tait EOS for positive pressure.
	 */
	class SymmetricTaitFluid : public WeaklyCompressibleFluid
	{
	protected:
		//determine the stiffness of the fluid
		int gamma_;

	public:
		///constructor
		explicit SymmetricTaitFluid(string fluid_name, SPHBody *body, Real rho_0 = 1.0,
			Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0);
		/** Constructor without body information. */
		explicit SymmetricTaitFluid(Real rho_0 = 1.0, Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0);
		virtual ~SymmetricTaitFluid() {};

		/** the interface for dynamical cast*/
		virtual SymmetricTaitFluid* PointToThisObject() override { return this; };

		virtual Real GetPressure(Real rho) override;
		virtual Real ReinitializeRho(Real p) override;
		virtual Real GetSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
	};

	/**
	 * @class Oldroyd_B_Fluid
	 * @brief linear EOS with relaxation time and polymetric viscosity.
	 */
	class Oldroyd_B_Fluid : public WeaklyCompressibleFluid
	{
	public:
		/** relaxation time */
		Real lambda_;
		/** polymeric viscosity */
		Real mu_p_;

		/** constructor */
		explicit Oldroyd_B_Fluid(string fluid_name, SPHBody *body, Real rho_0 = 1.0,
			Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0, 
			Real lambda = 1.0, Real mu_p = 0.0);
		virtual ~Oldroyd_B_Fluid() {};

		/** the interface for dynamical cast*/
		virtual Oldroyd_B_Fluid* PointToThisObject() override { return this; };

	};
}