/**
 * @file 	base_material.h
 * @brief 	This is the base classes of all materials. 
 *			Basically, it is a interface from which
 *			one can access devirved material by dynamic cast
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"
#include "base_body.h"
using namespace std;

namespace SPH {

	/** @class  Material
	 *  @brief Base of all materials
	 *  @details Note that the same material object can be shared by several
	 *  SPH bodies.
	*/
	class Material
	{
	protected:
		const string material_name_;

	public:
		explicit Material(string material_name)
			: material_name_(material_name) {};
		virtual ~Material() {};

		/** the interface for dynamical cast*/
		virtual Material* PointToThisObject() { return this; };
	};

	/** @class  Fluid
	 *  @brief Base calss  of all fluids
	*/
	class Fluid : public Material
	{
	protected:
		const string fluid_name_;

	public:
		/** reference density, sound speed, viscosity, thermal condution rate */
		const Real rho_0_, c_0_, mu_, k_;

		/** In the constructor, the base material is deleted. */
		Fluid(string fluid_name, SPHBody *body, Real rho_0 = 1.0, Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0)
			: Material(fluid_name), fluid_name_(fluid_name), rho_0_(rho_0), c_0_(c_0), mu_(mu), k_(k) 
		{
			delete body->base_material_;
			body->base_material_ = this;
		};
		virtual ~Fluid() {};

		/** the interface for dynamical cast*/
		virtual Fluid* PointToThisObject() override { return this; };

		Real GetReferenceSoundSpeed() { return c_0_; };
		Real GetReferenceDensity() { return rho_0_; };
		virtual Real GetPressure(Real rho) = 0; 
		virtual Real GetPressure(Real rho, Real rho_e) { return GetPressure(rho); };
		virtual Real ReinitializeRho(Real p) = 0;
		virtual Real GetSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;
		virtual Real RiemannSolverForPressure(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;
		virtual Real RiemannSolverForVelocity(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;
	};

	/** @class  Solid
	 *  @brief Base calss  of all solids
	*/
	class Solid : public Material
	{
	protected:
		const string solid_name_;

	public:
		/** reference density */
		const Real rho_0_;
		/** In the constructor, the base material is deleted. */
		Solid(string solid_name, SPHBody *body, Real rho_0 = 1.0)
			: Material(solid_name), solid_name_(solid_name), rho_0_(rho_0) 
		{
			delete body->base_material_;
			body->base_material_ = this;
		};
		virtual ~Solid() {};

		/** the interface for dynamical cast*/
		virtual Solid* PointToThisObject() override { return this; };
	};
}