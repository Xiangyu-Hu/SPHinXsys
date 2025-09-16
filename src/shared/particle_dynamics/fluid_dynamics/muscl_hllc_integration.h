#pragma once
#include "sphinxsys.h"
#include "muscl_reconstruction.hpp"
#include "eulerian_fluid_dynamics/eulerian_riemann_solver.h"

namespace SPH {
namespace fluid_dynamics {

struct MUSCLHLLCBridgeConfig {
    SecondOrderConfig muscl_cfg;
    bool use_hllc_dissipation_limiter = false;
    Real hllc_limiter_parameter = 0.0;
    Real gamma = 1.4;
    Real small = 1e-12;
};

// Bridge: reconstruct L/R interface primitives with MUSCL, then call HLLC.
class MUSCL_HLLC_Bridge {
  public:
    MUSCL_HLLC_Bridge(CompressibleFluid &fluid_i,
                      CompressibleFluid &fluid_j,
                      const MUSCLHLLCBridgeConfig &cfg);

    // Produce star-state at interface using MUSCL + HLLC.
#if SPH_NDIM == 2
    CompressibleFluidStarState getInterfaceState(
        const CompressibleFluidState &cell_i,
        const CompressibleFluidState &cell_j,
        const Vecd &xi, const Vecd &xj,
        const Vecd &xf, const Vecd &e_ij,
        const Vecd &grad_rho_i, const Vecd &grad_rho_j,
        const Vecd &grad_u_i,   const Vecd &grad_u_j,
        const Vecd &grad_v_i,   const Vecd &grad_v_j,
        const Vecd &grad_p_i,   const Vecd &grad_p_j) const;
#elif SPH_NDIM == 3
    CompressibleFluidStarState getInterfaceState(
        const CompressibleFluidState &cell_i,
        const CompressibleFluidState &cell_j,
        const Vecd &xi, const Vecd &xj,
        const Vecd &xf, const Vecd &e_ij,
        const Vecd &grad_rho_i, const Vecd &grad_rho_j,
        const Vecd &grad_u_i,   const Vecd &grad_u_j,
        const Vecd &grad_v_i,   const Vecd &grad_v_j,
        const Vecd &grad_w_i,   const Vecd &grad_w_j,
        const Vecd &grad_p_i,   const Vecd &grad_p_j) const;
#endif

  private:
    CompressibleFluid &fluid_i_;
    CompressibleFluid &fluid_j_;
    MUSCLHLLCBridgeConfig cfg_;

    // helpers
    inline Real total_energy_from_prim(Real rho, const Vecd &vel, Real p) const;
    inline void enforce_positivity(Real &rho, Real &p) const;
    inline Vecd unit_normal(const Vecd &n) const;
};

} // namespace fluid_dynamics
} // namespace SPH
