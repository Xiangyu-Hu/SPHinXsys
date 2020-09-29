/**
 * @file 	weakly_compressible_fluid.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "weakly_compressible_fluid.h"

using namespace std;

namespace SPH {
	//===============================================================//
	Real WeaklyCompressibleFluid::GetPressure(Real rho)
	{
		return p0_ * (rho / rho_0_ - 1.0);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::DensityFromPressure(Real p)
	{
		return rho_0_ * (p / p0_ + 1.0);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::GetSoundSpeed(Real p, Real rho)
	{
		return c_0_;
	}
	//===============================================================//
	Real WeaklyCompressibleFluid
		::RiemannSolverForVelocity(Real rhol, Real rhor, Real pl,
			Real pr, Real ul, Real ur)
	{
		Real rhol_cl = GetSoundSpeed(pl, rhol) * rhol;
		Real rhor_cr = GetSoundSpeed(pr, rhor) * rhor;

		return (rhol_cl * ul + rhor_cr * ur + pl - pr) / (rhol_cl + rhor_cr);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid
		::RiemannSolverForPressure(Real rhol, Real rhor, Real pl,
			Real pr, Real ul, Real ur)
	{
		Real rhol_cl = GetSoundSpeed(pl, rhol) * rhol;
		Real rhor_cr = GetSoundSpeed(pr, rhor) * rhor;
		Real clr = (rhol_cl + rhor_cr) / (rhol + rhor);

		return (rhol_cl * pr + rhor_cr * pl + rhol_cl * rhor_cr * (ul - ur)
			* SMIN(3.0 * SMAX((ul - ur) / clr, 0.0), 1.0)) / (rhol_cl + rhor_cr);
	}
	//===============================================================//
	Real SymmetricTaitFluid::GetPressure(Real rho)
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
	Real SymmetricTaitFluid::GetSoundSpeed(Real p, Real rho)
	{
		Real rho_ratio = rho / rho_0_;
		return rho_ratio > 1.0
			? sqrt((p0_ + Real(gamma_) * p) / rho)
			: sqrt((p0_ - Real(gamma_) * p) / rho);
	}
	//===============================================================//
}