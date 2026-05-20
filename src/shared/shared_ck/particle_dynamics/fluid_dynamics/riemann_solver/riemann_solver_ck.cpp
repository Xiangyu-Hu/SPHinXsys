#include "riemann_solver_ck.h"

namespace SPH
{
//=================================================================================================//
ImpedanceModel<WeaklyCompressibleFluid, WeaklyCompressibleFluid>::ImpedanceModel(
    WeaklyCompressibleFluid &fluid_i, WeaklyCompressibleFluid &fluid_j)
    : rho0c0_i_(fluid_i.ReferenceDensity() * fluid_i.ReferenceSoundSpeed()),
      rho0c0_j_(fluid_j.ReferenceDensity() * fluid_j.ReferenceSoundSpeed()),
      inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)),
      inv_rho0c0_ave_((rho0c0_i_ + rho0c0_j_) / (math::pow(rho0c0_i_, 2) + math::pow(rho0c0_j_, 2))),
      rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
      c0_ave_(0.5 * (fluid_i.ReferenceDensity() + fluid_j.ReferenceDensity()) * inv_rho0c0_ave_) {}
//=================================================================================================//
} // namespace SPH
