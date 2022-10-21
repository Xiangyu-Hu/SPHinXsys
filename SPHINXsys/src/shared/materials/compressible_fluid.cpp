#include "compressible_fluid.h"

using namespace std;

namespace SPH {
	//=============================================================================================//
	Real CompressibleFluid::getPressure(Real rho, Real rho_e)
	{
		return rho_e * (gamma_ - 1.0);
	}
	//=============================================================================================//
	Real CompressibleFluid::getSoundSpeed(Real p, Real rho)
	{
		return sqrt(gamma_ * p / rho);
	}
	//=============================================================================================//
}
