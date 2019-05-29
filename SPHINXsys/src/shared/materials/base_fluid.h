#pragma once

#include "base_data_package.h"

using namespace std;

namespace SPH {

	/*
	Fluid suggests the properties of a simple fluid
	*/
	class Fluid
	{
	protected:
		//name of the fluid
		const string fluid_name_;

	public:

		//reference density, sound speed, 
		//viscosity, thermal condution rate
		const Real rho_0_, c_0_, mu_, k_;

		//constructor
		Fluid(string fluid_name, Real rho_0 = 1.0, 
			Real c_0 = 1.0, Real mu = 0.0, Real k = 0.0) : fluid_name_(fluid_name),
			rho_0_(rho_0), c_0_(c_0), mu_(mu), k_(k)
		{

		};

		virtual ~Fluid() {};

		Real GetReferenceSoundSpeed() {
			return c_0_;
		};

		Real GetReferenceDensity() {
			return rho_0_;
		};

		virtual Real GetPressure(Real rho) = 0;
		virtual Real GetPressure(Real rho, Real rho_e)
		{
			return GetPressure(rho);
		};
		virtual Real ReinitializeRho(Real p) = 0;
		virtual Real GetSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;

		virtual Real RiemannSolverForPressure(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;
		virtual Real RiemannSolverForVelocity(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;

	};

}