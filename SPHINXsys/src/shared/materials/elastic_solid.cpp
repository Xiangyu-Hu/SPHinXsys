#include "elastic_solid.h"

namespace SPH {
	//===================================================================//
	ElasticSolid::ElasticSolid(string elastic_solid_name,
		Real rho_0, Real E_0, Real nu, Real eta_0)
		: Solid(elastic_solid_name, rho_0), 
		E_0_(E_0), nu_(nu), eta_0_(eta_0)
	{	

	}
	//===================================================================//
	Matd ElasticSolid::ConstitutiveRelation(Matd &F, Real local_G, Real local_lambda)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = local_lambda * strain.trace() *  Matd(1.0)
			+ 2.0*local_G * strain;
		return F * sigmaPK2;
	}
	//===================================================================//
	Matd ElasticSolid::DampingStress(Matd &F, Matd &dF_dt, Real local_eta)
	{
		Matd strain_rate = 0.5*(~dF_dt*F + ~F * dF_dt);
		Matd sigmaPK2 = local_eta *strain_rate;
		return F * sigmaPK2;
	}
	//===================================================================//
	Real ElasticSolid::GetArtificalViscosity(Real rho,
		Real sound_speed, Real smoothing_length)
	{
		return 0.5*rho*sound_speed*smoothing_length;
	}
	//===================================================================//
	Real ElasticSolid::GetSoundSpeed(Real rho_0, Real E, Real nu)
	{
		return  sqrt(E / 3.0 / (1.0 - 2.0 * nu) / rho_0);
	}
	//===================================================================//
	Real ElasticSolid::GetShearModulus(Real E, Real nu)
	{
		return 0.5 * E / (1.0 + nu);
	}
	//===================================================================//
	Real ElasticSolid::GetLambda(Real E, Real nu)
	{
		return nu * E / (1.0 + nu) / (1.0 - 2.0 * nu);
	}
	//===================================================================//
	NeoHookeanSolid::NeoHookeanSolid(string elastic_solid_name,
		Real rho_0, Real E_0, Real nu, Real eta_0)
		: ElasticSolid(elastic_solid_name, rho_0, E_0, nu, eta_0)
	{

	}
	//===================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd &F, Real local_G, Real local_lambda)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = local_G * Matd(1.0)
			+ (local_lambda*log(det(F)) - local_G)*inverse(right_cauchy);
		return F * sigmaPK2;
	}
	//===================================================================//
	FeneNeoHookeanSolid::FeneNeoHookeanSolid(string elastic_solid_name,
		Real rho_0, Real E_0, Real nu, Real eta_0, Real j1_m)
		: ElasticSolid(elastic_solid_name, rho_0, E_0, nu, eta_0),
		j1_m_(j1_m)
	{

	}
	//===================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd &F, Real local_G, Real local_lambda)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = local_G / (1.0 - 2.0*strain.trace() / j1_m_)*Matd(1.0)
			+ (local_lambda*log(det(F)) - local_G)*inverse(right_cauchy);
		return F * sigmaPK2;
	}
	//===================================================================//
	Muscle::Muscle(string elastic_solid_name, Real a_0[4], Real b_0[4],
		Real rho_0, Real poisson, Real eta_0)
		: ElasticSolid(elastic_solid_name, rho_0, 2.0*a_0[0]*b_0[0]*(1.0 + poisson), poisson, eta_0)
	{
		for (size_t i = 0; i != 4; ++i)
		{
			a_0_[i] = a_0[i];
			b_0_[i] = b_0[i];
		}
	}
	//===================================================================//
	Matd Muscle::ConstitutiveRelation(Matd &F,
		Real local_lambda, Real a[4], Real b[4], Vecd local_f0, Vecd local_s0,
		Matd local_f0f0, Matd local_s0s0, Matd local_f0s0)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy*local_f0, local_f0) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy*local_s0, local_s0) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy*local_f0, local_s0);
		Real ln_J = log(det(F));
		Real I_1_1 = right_cauchy.trace() - 3.0 - 2.0*ln_J;
		Matd sigmaPK2 = a[0]*exp(b[0] * I_1_1)*Matd(1.0)
			+ (local_lambda*ln_J - a[0] * exp(b[0] * I_1_1))*inverse(right_cauchy)
			+ 2.0* a[1]*I_ff_1*exp(b[1]*I_ff_1*I_ff_1)*local_f0f0
			+ 2.0* a[2]*I_ss_1*exp(b[2]*I_ss_1*I_ss_1)*local_s0s0
			+ a[3]* I_fs*exp(b[3]* I_fs*I_fs)*local_f0s0;
		return F * sigmaPK2;
	}
}