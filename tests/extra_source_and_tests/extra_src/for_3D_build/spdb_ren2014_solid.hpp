#ifndef SPDB_REN2014_SOLID_HPP
#define SPDB_REN2014_SOLID_HPP

#include "spdb_ren2014_solid.h"

namespace SPH
{
//=================================================================================================//
template <typename ExecutionPolicy>
Ren2014X20Cr13Solid::ConstituteKernel::ConstituteKernel(const ExecutionPolicy &ex_policy, Ren2014X20Cr13Solid &encloser)
    : NeoHookeanSolid::ConstituteKernel(ex_policy, encloser), temperature_K_(encloser.temperature_K_),
      inelastic_enabled_(encloser.inelastic_enabled_),
      inverse_plastic_strain_(encloser.dv_inverse_plastic_strain_->DelegatedData(ex_policy)),
      inelastic_eq_(encloser.dv_inelastic_eq_->DelegatedData(ex_policy)),
      range_flag_(encloser.dv_range_flag_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Matd Ren2014X20Cr13Solid::ConstituteKernel::ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt)
{
    return relax_ren2014_elastic_left_cauchy(deformation, index_i, dt, G0_, temperature_K_, inelastic_enabled_,
                                             inverse_plastic_strain_, inelastic_eq_, range_flag_);
}
//=================================================================================================//
} // namespace SPH

#endif
