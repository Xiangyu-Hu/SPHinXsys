#include "muscl_hllc_integration.h"
#include "eulerian_fluid_dynamics/eulerian_riemann_solver.h"
#include <algorithm>
#include <cmath>

namespace SPH {
namespace fluid_dynamics {

using namespace SPH;

MUSCL_HLLC_Bridge::MUSCL_HLLC_Bridge(CompressibleFluid &fluid_i,
                                     CompressibleFluid &fluid_j,
                                     const MUSCLHLLCBridgeConfig &cfg)
: fluid_i_(fluid_i), fluid_j_(fluid_j), cfg_(cfg) {}

inline void MUSCL_HLLC_Bridge::enforce_positivity(Real &rho, Real &p) const {
    rho = std::max(rho, cfg_.small);
    p   = std::max(p,   cfg_.small);
}

inline Real MUSCL_HLLC_Bridge::total_energy_from_prim(Real rho, const Vecd &vel, Real p) const {
    Real u2 = vel.squaredNorm();
    // Total energy density (J/m^3)
    return p/(cfg_.gamma - 1.0) + 0.5 * rho * u2;
}

inline Vecd MUSCL_HLLC_Bridge::unit_normal(const Vecd &n) const {
    Real nn = n.norm();
    return (nn > 0) ? (n / nn) : n;
}

#if SPH_NDIM == 2
CompressibleFluidStarState MUSCL_HLLC_Bridge::getInterfaceState(
    const CompressibleFluidState &cell_i,
    const CompressibleFluidState &cell_j,
    const Vecd &xi, const Vecd &xj,
    const Vecd &xf, const Vecd &e_ij_in,
    const Vecd &grad_rho_i, const Vecd &grad_rho_j,
    const Vecd &grad_u_i,   const Vecd &grad_u_j,
    const Vecd &grad_v_i,   const Vecd &grad_v_j,
    const Vecd &grad_p_i,   const Vecd &grad_p_j) const
{
    const Vecd e_ij = unit_normal(e_ij_in);

    // 1) Pack MUSCL primitives from cell-centered states
    Primitives Pi{cell_i.rho_, cell_i.vel_, cell_i.p_, cell_i.E_};
    Primitives Pj{cell_j.rho_, cell_j.vel_, cell_j.p_, cell_j.E_};

    // 2) MUSCL reconstruction of interface primitives (L/R)
    LR lr = reconstruct_primitives_muscl(
        Pi, Pj,
        grad_rho_i, grad_rho_j,
        grad_u_i,   grad_u_j,
        grad_v_i,   grad_v_j,
        grad_p_i,   grad_p_j,
        xi, xj, xf, e_ij,
        cfg_.muscl_cfg
    );
    
    // Debug: MUSCL reconstruction results
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "MUSCL L: rho=" << lr.L.rho << " p=" << lr.L.p << " vel=" << lr.L.vel.transpose() << std::endl;
    std::cout << "MUSCL R: rho=" << lr.R.rho << " p=" << lr.R.p << " vel=" << lr.R.vel.transpose() << std::endl;
    #endif

    // 3) Positivity guard at interface primitives
    Real rhoL = lr.L.rho, rhoR = lr.R.rho;
    Real  pL  = lr.L.p,    pR  = lr.R.p;
    enforce_positivity(rhoL, pL);
    enforce_positivity(rhoR, pR);

    Vecd vL = lr.L.vel;
    Vecd vR = lr.R.vel;

    // 4) Build CompressibleFluidState for HLLC with EOS-consistent energy
    Real EL = total_energy_from_prim(rhoL, vL, pL);
    Real ER = total_energy_from_prim(rhoR, vR, pR);

    // Locals provide lifetime for reference-based state wrappers during the call
    Real rhoL_local = rhoL, pL_local = pL, EL_local = EL;
    Real rhoR_local = rhoR, pR_local = pR, ER_local = ER;
    Vecd vL_local = vL, vR_local = vR;

    CompressibleFluidState Ls(rhoL_local, vL_local, pL_local, EL_local);
    CompressibleFluidState Rs(rhoR_local, vR_local, pR_local, ER_local);
    
    // Debug output
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "Bridge L: rho=" << Ls.rho_ << " p=" << Ls.p_ << " vel=" << Ls.vel_.transpose() << " E=" << Ls.E_ << std::endl;
    std::cout << "Bridge R: rho=" << Rs.rho_ << " p=" << Rs.p_ << " vel=" << Rs.vel_.transpose() << " E=" << Rs.E_ << std::endl;
    std::cout << "Interface normal e_ij: " << e_ij.transpose() << std::endl;
    #endif

    // 5) Call HLLC (optionally with dissipation limiter)
    CompressibleFluidStarState star = [&]() {
        if (cfg_.use_hllc_dissipation_limiter) {
            HLLCWithLimiterRiemannSolver solver(fluid_i_, fluid_j_, cfg_.hllc_limiter_parameter);
            return solver.getInterfaceState(Ls, Rs, e_ij);
        } else {
            HLLCRiemannSolver solver(fluid_i_, fluid_j_);
            return solver.getInterfaceState(Ls, Rs, e_ij);
        }
    }();

    // Debug: HLLC solver results
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "HLLC Star: rho=" << star.rho_ << " p=" << star.p_ << " vel=" << star.vel_.transpose() << " E=" << star.E_ << std::endl;
    #endif
    
    return star;
}
#elif SPH_NDIM == 3
CompressibleFluidStarState MUSCL_HLLC_Bridge::getInterfaceState(
    const CompressibleFluidState &cell_i,
    const CompressibleFluidState &cell_j,
    const Vecd &xi, const Vecd &xj,
    const Vecd &xf, const Vecd &e_ij_in,
    const Vecd &grad_rho_i, const Vecd &grad_rho_j,
    const Vecd &grad_u_i,   const Vecd &grad_u_j,
    const Vecd &grad_v_i,   const Vecd &grad_v_j,
    const Vecd &grad_w_i,   const Vecd &grad_w_j,
    const Vecd &grad_p_i,   const Vecd &grad_p_j) const
{
    const Vecd e_ij = unit_normal(e_ij_in);

    Primitives Pi{cell_i.rho_, cell_i.vel_, cell_i.p_, cell_i.E_};
    Primitives Pj{cell_j.rho_, cell_j.vel_, cell_j.p_, cell_j.E_};

 
    LR lr = reconstruct_primitives_muscl(
        Pi, Pj,
        grad_rho_i, grad_rho_j,
        grad_u_i,   grad_u_j,
        grad_v_i,   grad_v_j,
        grad_w_i,   grad_w_j,
        grad_p_i,   grad_p_j,
        xi, xj, xf, e_ij,
        cfg_.muscl_cfg
    );
    
    // Debug: MUSCL reconstruction results
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "MUSCL L: rho=" << lr.L.rho << " p=" << lr.L.p << " vel=" << lr.L.vel.transpose() << std::endl;
    std::cout << "MUSCL R: rho=" << lr.R.rho << " p=" << lr.R.p << " vel=" << lr.R.vel.transpose() << std::endl;
    #endif

    Real rhoL = lr.L.rho, rhoR = lr.R.rho;
    Real  pL  = lr.L.p,    pR  = lr.R.p;
    enforce_positivity(rhoL, pL);
    enforce_positivity(rhoR, pR);

    Vecd vL = lr.L.vel;
    Vecd vR = lr.R.vel;

    Real EL = total_energy_from_prim(rhoL, vL, pL);
    Real ER = total_energy_from_prim(rhoR, vR, pR);

    Real rhoL_local = rhoL, pL_local = pL, EL_local = EL;
    Real rhoR_local = rhoR, pR_local = pR, ER_local = ER;
    Vecd vL_local = vL, vR_local = vR;

    CompressibleFluidState Ls(rhoL_local, vL_local, pL_local, EL_local);
    CompressibleFluidState Rs(rhoR_local, vR_local, pR_local, ER_local);

    // Debug output
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "Bridge L: rho=" << Ls.rho_ << " p=" << Ls.p_ << " vel=" << Ls.vel_.transpose() << " E=" << Ls.E_ << std::endl;
    std::cout << "Bridge R: rho=" << Rs.rho_ << " p=" << Rs.p_ << " vel=" << Rs.vel_.transpose() << " E=" << Rs.E_ << std::endl;
    std::cout << "Interface normal e_ij: " << e_ij.transpose() << std::endl;
    #endif

    // Call HLLC (optionally with dissipation limiter)
    CompressibleFluidStarState star = [&]() {
        if (cfg_.use_hllc_dissipation_limiter) {
            HLLCWithLimiterRiemannSolver solver(fluid_i_, fluid_j_, cfg_.hllc_limiter_parameter);
            return solver.getInterfaceState(Ls, Rs, e_ij);
        } else {
            HLLCRiemannSolver solver(fluid_i_, fluid_j_);
            return solver.getInterfaceState(Ls, Rs, e_ij);
        }
    }();

    // Debug: HLLC solver results
    #ifdef SPH_DEBUG_MUSCL
    std::cout << "HLLC Star: rho=" << star.rho_ << " p=" << star.p_ << " vel=" << star.vel_.transpose() << " E=" << star.E_ << std::endl;
    #endif
    
    return star;
}
#endif

} // namespace fluid_dynamics
} // namespace SPH
