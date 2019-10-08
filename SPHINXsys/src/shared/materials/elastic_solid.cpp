/** 
 * @file elastic_solid.cpp
 * @brief These are functions defined in elastic_solid.h
 * @author Chi Zhang and Xiangyu Hu
 * @version  0.1
 * @version  0.2.1
 *			 add the electrophysiology to muscle body.
 */
#include "elastic_solid.h"
//=================================================================================================//
namespace SPH {
	//===================================================================//
	ElasticSolid::ElasticSolid(string elastic_solid_name, SPHBody *body,
		Real rho_0, Real E_0, Real nu, Real eta_0)
		: Solid(elastic_solid_name, body, rho_0), 
		E_0_(E_0), nu_(nu), eta_0_(eta_0)
	{	
		lambda_0_ = SetLambda(E_0_, nu_);
		G_0_ = SetShearModulus(E_0_, nu_);
		c_0_ = SetSoundSpeed(rho_0, E_0, nu);
	}
//=================================================================================================//
	Real ElasticSolid::GetSoundSpeed(size_t particle_index_i)
	{
		return c_0_;
	}
//=================================================================================================//
	Matd ElasticSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = lambda_0_ * strain.trace() *  Matd(1.0)
			+ 2.0*G_0_ * strain;
		return F * sigmaPK2;
	}
//=================================================================================================//
	Matd ElasticSolid::DampingStress(Matd &F, Matd &dF_dt, size_t particle_index_i)
	{
		Matd strain_rate = 0.5*(~dF_dt*F + ~F * dF_dt);
		Matd sigmaPK2 = eta_0_ *strain_rate;
		return F * sigmaPK2;
	}
//=================================================================================================//
	Matd ElasticSolid::NumericalDampingStress(Matd &F, Matd &dF_dt, Real numerical_viscoisty)
	{
		Matd strain_rate = 0.5*(~dF_dt*F + ~F * dF_dt);
		Matd sigmaPK2 = numerical_viscoisty * strain_rate;
		return F * sigmaPK2;
	}	
//=================================================================================================//
	Real ElasticSolid::GetArtificalViscosity(Real rho,
		Real sound_speed, Real smoothing_length)
	{
		return 0.5*rho*sound_speed*smoothing_length;
	}
//=================================================================================================//
	Real ElasticSolid::SetSoundSpeed(Real rho_0, Real E, Real nu)
	{
		return  sqrt(E / 3.0 / (1.0 - 2.0 * nu) / rho_0);
	}
//=================================================================================================//
	Real ElasticSolid::SetShearModulus(Real E, Real nu)
	{
		return 0.5 * E / (1.0 + nu);
	}
//=================================================================================================//
	Real ElasticSolid::SetLambda(Real E, Real nu)
	{
		return nu * E / (1.0 + nu) / (1.0 - 2.0 * nu);
	}
	//===================================================================//
	NeoHookeanSolid::NeoHookeanSolid(string elastic_solid_name, SPHBody *body,
		Real rho_0, Real E_0, Real nu, Real eta_0)
		: ElasticSolid(elastic_solid_name, body, rho_0, E_0, nu, eta_0)
	{

	}
//=================================================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = G_0_ * Matd(1.0)
			+ (lambda_0_*log(det(F)) - G_0_)*inverse(right_cauchy);
		return F * sigmaPK2;
	}
	//===================================================================//
	FeneNeoHookeanSolid::FeneNeoHookeanSolid(string elastic_solid_name, SPHBody *body,
		Real rho_0, Real E_0, Real nu, Real eta_0, Real j1_m)
		: ElasticSolid(elastic_solid_name, body, rho_0, E_0, nu, eta_0),
		j1_m_(j1_m)
	{

	}
//=================================================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = G_0_ / (1.0 - 2.0*strain.trace() / j1_m_)*Matd(1.0)
			+ (lambda_0_*log(det(F)) - G_0_)*inverse(right_cauchy);
		return F * sigmaPK2;
	}
//=================================================================================================//
	Muscle::Muscle(string elastic_solid_name,SPHBody *body, Real a_0[4], Real b_0[4], Vecd d_0, 
		Real rho_0, Real poisson, Real eta_0)
		: ElasticSolid(elastic_solid_name, body, 
			rho_0, 2.0*a_0[0]*b_0[0]*(1.0 + poisson), poisson, eta_0)
	{
		for (size_t i = 0; i != 4; ++i)
		{
			a_0_[i] = a_0[i];
			b_0_[i] = b_0[i];
		}
		d_0_ = d_0;
	}
//=================================================================================================//
	void Muscle::SetupLocalProperties(ElasticSolidParticles &elasticsolid_particles)
	{
		cout << "\n This function Muscle::SetupLocalProperties is not done yet. Exit the program! \n";
		exit(0);
	}
//=================================================================================================//
	Matd Muscle::ConstitutiveRelation(Matd &F, size_t i)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy*f0_[i], f0_[i]) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy*s0_[i], s0_[i]) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy*f0_[i], s0_[i]);
		Real ln_J = log(det(F));
		Real I_1_1 = right_cauchy.trace() - 3.0 - 2.0*ln_J;
		Matd sigmaPK2 = a_0_[0]*exp(b_0_[0] * I_1_1)*Matd(1.0)
			+ (lambda_0_*ln_J - a_0_[0] * exp(b_0_[0] * I_1_1))*inverse(right_cauchy)
			+ 2.0* a_0_[1]*I_ff_1*exp(b_0_[1]*I_ff_1*I_ff_1)*f0f0_[i]
			+ 2.0* a_0_[2]*I_ss_1*exp(b_0_[2]*I_ss_1*I_ss_1)*s0s0_[i]
			+ a_0_[3]* I_fs*exp(b_0_[3]* I_fs*I_fs)*f0s0_[i];
		return F * sigmaPK2;
	}
//=================================================================================================//
	Matd Muscle::getDiffussionTensor(Vecd pnt)
	{
		Matd matrix_unit(1.0);
		Vecd f0(0.0, 1.0);
		Matd diffusion_tensor = d_0_[0] * matrix_unit - d_0_[1] * SimTK::outer(f0, f0);

		return diffusion_tensor;
	}
//=================================================================================================//
	Real Muscle::getDiffusionTensorTrace(Vecd pnt)
	{
		Matd matrix_unit(1.0);
		Vecd f0(0.0, 1.0);
		Matd diffusion_tensor = d_0_[0] * matrix_unit + d_0_[1] * SimTK::outer(f0, f0);
		return diffusion_tensor.trace();
	}
//=================================================================================================//
}