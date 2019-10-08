#include "weakly_compressible_fluid.h"

using namespace std;

namespace SPH {
	//===============================================================//
	WeaklyCompressibleFluid::WeaklyCompressibleFluid(string fluid_name, 
		SPHBody *body, Real rho0, Real c0, Real mu, Real k) 
		: Fluid(fluid_name, body, rho0, c0, mu, k)
	{
		p0_ = rho_0_* c_0_ * c_0_;
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::GetPressure(Real rho)
	{
		return p0_ *(rho/rho_0_ - 1.0);
	}
	//===============================================================//
	Real WeaklyCompressibleFluid::ReinitializeRho(Real p)
	{
		return rho_0_*(p/ p0_ + 1.0);
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
		Real rhol_cl = GetSoundSpeed(pl, rhol)*rhol;
		Real rhor_cr = GetSoundSpeed(pr, rhor)*rhor;
		Real clr =(rhol_cl + rhor_cr)/ (rhol + rhor);

		return (rhol_cl * pr + rhor_cr * pl + rhol_cl * rhor_cr * (ul - ur)
			*SMIN(3.0*SMAX((ul - ur) / clr, 0.0), 1.0)) / (rhol_cl + rhor_cr);
	}
	//===============================================================//
	SymmetricTaitFluid
		::SymmetricTaitFluid(string fluid_name, SPHBody *body, 
			Real rho_0, Real c_0, Real mu, Real k)
		:WeaklyCompressibleFluid(fluid_name, body, rho_0, c_0, mu, k),
		gamma_(2)
	{

	}
	//===============================================================//
	Real SymmetricTaitFluid::GetPressure(Real rho)
	{
		Real rho_ratio = rho / rho_0_;
		return rho_ratio > 1.0 
			? p0_ * (powern(rho_ratio, gamma_) - 1.0) / Real(gamma_)
			: - p0_ * (powern(1.0/rho_ratio, gamma_) - 1.0) / Real(gamma_);
	}
	//===============================================================//
	Real SymmetricTaitFluid::ReinitializeRho(Real p)
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
			? sqrt((p0_ + Real(gamma_)*p) / rho)
			: sqrt((p0_ - Real(gamma_) * p) / rho);
	}
	//===============================================================//
	Oldroyd_B_Fluid
		::Oldroyd_B_Fluid(string fluid_name, SPHBody *body, Real rho_0,
			Real c_0, Real mu, Real k,	Real lambda, Real mu_p)
		: WeaklyCompressibleFluid(fluid_name, body, rho_0, c_0, mu, k),
		lambda_(lambda), mu_p_(mu_p)
	{
	//===============================================================//
	}
}