#include "fluid_shell_interaction.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
    	//=================================================================================================//
		void DensitySummationComplexWithShell::interaction(size_t index_i, Real dt)
		{
			BaseDensitySummationComplex<DensitySummationInner>::interaction(index_i, dt);
			Real sigma = spacing_ref_ * BaseDensitySummationComplex<DensitySummationInner>::ContactSummation(index_i);
			rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
		}
		//=================================================================================================//
		void DensitySummationComplexAdaptiveWithShell::interaction(size_t index_i, Real dt)
		{
			BaseDensitySummationComplex<DensitySummationInnerAdaptive>::interaction(index_i, dt);
			Real sigma = spacing_ref_ *  BaseDensitySummationComplex<DensitySummationInnerAdaptive>::ContactSummation(index_i);
			rho_sum_[index_i] += sigma * rho0_ * rho0_ / mass_[index_i] /
								 sph_adaptation_.ReferenceNumberDensity(h_ratio_[index_i]);
		}
		//=================================================================================================//
    }
//=================================================================================================//
}