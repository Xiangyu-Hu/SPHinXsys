#include "riemann_solver_ck.h"

namespace SPH
{
//=================================================================================================//
ImpedanceModel<WeaklyCompressibleFluid, WeaklyCompressibleFluid>::ImpedanceModel(
    WeaklyCompressibleFluid &fluid_i, WeaklyCompressibleFluid &fluid_j)
    : rho0_i_(fluid_i.ReferenceDensity()), rho0_j_(fluid_j.ReferenceDensity()),
      rho0c0_i_(rho0_i_ * fluid_i.ReferenceSoundSpeed()),
      rho0c0_j_(rho0_j_ * fluid_j.ReferenceSoundSpeed()),
      inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)),
      inv_rho0c0_ave_((rho0c0_i_ + rho0c0_j_) / (math::pow(rho0c0_i_, 2) + math::pow(rho0c0_j_, 2))),
      rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
      inv_c0_ave_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_) {}
//=================================================================================================//
} // namespace SPH
