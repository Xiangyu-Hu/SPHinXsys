/**
 * @file 	weakly_compressible_fluid.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "weakly_compressible_fluid.h"

using namespace std;

namespace SPH {
	//===============================================================//
	Real WeaklyCompressibleFluid::getPressure(Real rho)
	{
		return p0_ * (rho / rho_0_ - 1.0);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::DensityFromPressure(Real p)
	{
		return rho_0_ * (p / p0_ + 1.0);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::getSoundSpeed(Real p, Real rho)
	{
		return c_0_;
	}
	//===============================================================//
	Real SymmetricTaitFluid::getPressure(Real rho)
	{
		Real rho_ratio = rho / rho_0_;
		return rho_ratio > 1.0
			? p0_ * (powern(rho_ratio, gamma_) - 1.0) / Real(gamma_)
			: -p0_ * (powern(1.0 / rho_ratio, gamma_) - 1.0) / Real(gamma_);
	}
	//===============================================================//
	Real SymmetricTaitFluid::DensityFromPressure(Real p)
	{
		return p > 0.0
			? rho_0_ * pow(1.0 + Real(gamma_) * p / p0_, 1.0 / Real(gamma_))
			: rho_0_ / pow(1.0 - Real(gamma_) * p / p0_, 1.0 / Real(gamma_));
	}
	//===============================================================//
	Real SymmetricTaitFluid::getSoundSpeed(Real p, Real rho)
	{
		Real rho_ratio = rho / rho_0_;
		return rho_ratio > 1.0
			? sqrt((p0_ + Real(gamma_) * p) / rho)
			: sqrt((p0_ - Real(gamma_) * p) / rho);
	}
	//===============================================================//
}
