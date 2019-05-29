#pragma once

#include "base_data_package.h"

using namespace std;

namespace SPH {

	/*
	Fluid suggests the properties of a simple solid
	*/
	class Solid
	{
	protected:
		//name of the solid
		const string solid_name_;

	public:

		//reference density, sound speed, 
		//viscosity, thermal condution rate
		const Real rho_0_;

		//constructor
		Solid(string solid_name, Real rho_0 = 1.0) : solid_name_(solid_name),
			rho_0_(rho_0)
		{

		};

		virtual ~Solid() {};

	};

}