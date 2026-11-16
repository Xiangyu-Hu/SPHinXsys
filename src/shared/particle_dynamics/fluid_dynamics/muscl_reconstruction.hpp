#pragma once
 
#include "data_type.h"
 #include <algorithm>
 #include <cmath>
 #include <string>
 #include <Eigen/Core> 
 
 namespace SPH {
 namespace fluid_dynamics {

struct Primitives
{
    Real rho;     // density
    Vecd vel;     // velocity (u,v[,w])
    Real p;       // pressure
    Real E;       // total energy
};

/// Left/Right pair (interface states)
struct LR
{
    Primitives L;
    Primitives R;
};

// Slope limiter enum
enum class SlopeLimiter : int {
    None = 0,
    Minmod,
    MC,
    VanLeer
};

/// Config for reconstruction
struct SecondOrderConfig
{
    SlopeLimiter limiter          = SlopeLimiter::Minmod;
    bool         positivity       = true;
    // Clamp face values to [min(Ui,Uj), max(Ui,Uj)]. Only valid when x_iface lies on the
    // segment between xi and xj (e.g. classic collinear FV faces) -- for off-line interface
    // points (wall/ghost reconstructions, general SPH neighbor geometry) this bound does not
    // hold and would corrupt an otherwise-exact reconstruction, so it defaults to off.
    bool         monotonicity     = false;
    bool         coupled_limiting = true;  // one limiter factor for all primitives (recommended)
    Real         small            = 1e-12;
    Real         gamma            = 1.4;   // ideal-gas gamma for EOS-based energy
};

/// Scalar limiters (standalone, inlined)
inline Real limiter_minmod(Real a, Real b)
{
    if (a * b <= (Real)0) return (Real)0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

inline Real limiter_mc(Real a, Real b)
{
    if (a * b <= (Real)0) return (Real)0;
    const Real s = (a > 0) ? (Real)1 : (Real)(-1);
    const Real aa = std::abs(a), bb = std::abs(b);
    return s * std::min({ (Real)2*aa, (Real)2*bb, (Real)0.5*(aa+bb) });
}

inline Real limiter_vanleer(Real a, Real b)
{
    if (a * b <= (Real)0) return (Real)0;
    return (Real)2.0 * a * b / (a + b);
}

/// Pick limiter
inline Real apply_limiter(SlopeLimiter type, Real a, Real b)
{
    switch (type) {
        case SlopeLimiter::Minmod: return limiter_minmod(a, b);
        case SlopeLimiter::MC:     return limiter_mc(a, b);
        case SlopeLimiter::VanLeer:return limiter_vanleer(a, b);
        case SlopeLimiter::None:
        default:                   return a; // no limiting (use left slope)
    }
}

/// Compute left/right limiter factors phi_i, phi_j for one scalar stencil.
/// i-side uses forward jump (Uj-Ui); j-side uses backward jump (Ui-Uj) and slope toward interface.
template <typename Vec>
inline std::pair<Real, Real>
compute_muscl_limiter_factors(Real Ui, Real Uj,
                              const Vec& gradUi, const Vec& gradUj,
                              const Vec& xi, const Vec& xj,
                              const SecondOrderConfig& cfg)
{
    const Vec dx_vec   = xj - xi;
    const Vec dx_back  = xi - xj;          // from j toward i (interface direction for j)
    const Real du_fwd  = Uj - Ui;
    const Real du_back = Ui - Uj;

    const Real si = gradUi.dot(dx_vec);
    const Real sj = gradUj.dot(dx_back);   // backward slope at j

    auto safe_div = [&](Real num, Real den){ return (std::abs(den) > 1e-14) ? (num/den) : 0.0; };

    const Real phi_i_raw = apply_limiter(cfg.limiter, si, du_fwd);
    const Real phi_j_raw = apply_limiter(cfg.limiter, sj, du_back);

    const Real phi_i = safe_div(phi_i_raw, (std::abs(si) > 1e-14 ? si : (Real)1));
    const Real phi_j = safe_div(phi_j_raw, (std::abs(sj) > 1e-14 ? sj : (Real)1));

    return {phi_i, phi_j};
}

/// Clamp reconstructed L/R face values to [min(Ui,Uj), max(Ui,Uj)] when requested.
/// See SecondOrderConfig::monotonicity for when this bound is valid.
inline void clamp_to_bounds(Real& UL, Real& UR, Real Ui, Real Uj, const SecondOrderConfig& cfg)
{
    if (!cfg.monotonicity) return;
    const Real Umin = std::min(Ui, Uj);
    const Real Umax = std::max(Ui, Uj);
    UL = std::clamp(UL, Umin, Umax);
    UR = std::clamp(UR, Umin, Umax);
}

template <typename Vec>
inline std::pair<Real, Real>
reconstruct_scalar_muscl(Real Ui, const Vec& gradUi,
                         Real Uj, const Vec& gradUj,
                         const Vec& xi, const Vec& xj,
                         const Vec& x_iface,
                         const SecondOrderConfig& cfg)
{
    const Vec di_vec = x_iface - xi;
    const Vec dj_vec = x_iface - xj;

    if (cfg.limiter == SlopeLimiter::None) {
        Real UL = gradUi.dot(di_vec) + Ui;
        Real UR = gradUj.dot(dj_vec) + Uj;
        clamp_to_bounds(UL, UR, Ui, Uj, cfg);
        return {UL, UR};
    }

    const auto phis = compute_muscl_limiter_factors(Ui, Uj, gradUi, gradUj, xi, xj, cfg);
    Real UL = Ui + phis.first  * gradUi.dot(di_vec);
    Real UR = Uj + phis.second * gradUj.dot(dj_vec);

    clamp_to_bounds(UL, UR, Ui, Uj, cfg);

    return {UL, UR};
}


#if SPH_NDIM == 2
template <typename Vec>
inline LR reconstruct_primitives_muscl(const Primitives& Pi, const Primitives& Pj,
                                       const Vec& grad_rho_i, const Vec& grad_rho_j,
                                       const Vec& grad_u_i,   const Vec& grad_u_j,
                                       const Vec& grad_v_i,   const Vec& grad_v_j,
                                       const Vec& grad_p_i,   const Vec& grad_p_j,
                                       const Vec& xi, const Vec& xj,
                                       const Vec& x_iface,
                                       const SecondOrderConfig& cfg)
#elif SPH_NDIM == 3
template <typename Vec>
inline LR reconstruct_primitives_muscl(const Primitives& Pi, const Primitives& Pj,
                                       const Vec& grad_rho_i, const Vec& grad_rho_j,
                                       const Vec& grad_u_i,   const Vec& grad_u_j,
                                       const Vec& grad_v_i,   const Vec& grad_v_j,
                                       const Vec& grad_w_i,   const Vec& grad_w_j,
                                       const Vec& grad_p_i,   const Vec& grad_p_j,
                                       const Vec& xi, const Vec& xj,
                                       const Vec& x_iface,
                                       const SecondOrderConfig& cfg)
#endif
{
    LR out;
    const Vec di_vec = x_iface - xi;
    const Vec dj_vec = x_iface - xj;

    Real phi_i = (Real)1, phi_j = (Real)1;
    if (cfg.limiter != SlopeLimiter::None) {
        if (cfg.coupled_limiting) {
            // One limiter factor (from density) for all primitives avoids inconsistent face states.
            const auto phis = compute_muscl_limiter_factors(
                Pi.rho, Pj.rho, grad_rho_i, grad_rho_j, xi, xj, cfg);
            phi_i = phis.first;
            phi_j = phis.second;
        }
    }

    auto reconstruct_component = [&](Real Ui, Real Uj, const Vec& gradUi, const Vec& gradUj) {
        Real local_phi_i = phi_i, local_phi_j = phi_j;
        if (cfg.limiter != SlopeLimiter::None && !cfg.coupled_limiting) {
            const auto phis = compute_muscl_limiter_factors(Ui, Uj, gradUi, gradUj, xi, xj, cfg);
            local_phi_i = phis.first;
            local_phi_j = phis.second;
        } else if (cfg.limiter == SlopeLimiter::None) {
            local_phi_i = (Real)1;
            local_phi_j = (Real)1;
        }
        Real UL = Ui + local_phi_i * gradUi.dot(di_vec);
        Real UR = Uj + local_phi_j * gradUj.dot(dj_vec);
        clamp_to_bounds(UL, UR, Ui, Uj, cfg);
        return std::pair<Real, Real>{UL, UR};
    };

    // density
    {
        auto lr = reconstruct_component(Pi.rho, Pj.rho, grad_rho_i, grad_rho_j);
        out.L.rho = lr.first;  out.R.rho = lr.second;
    }

    // velocity components
    {
        auto lr_u = reconstruct_component(Pi.vel[0], Pj.vel[0], grad_u_i, grad_u_j);
        auto lr_v = reconstruct_component(Pi.vel[1], Pj.vel[1], grad_v_i, grad_v_j);

#if SPH_NDIM == 2
        out.L.vel = Vecd(lr_u.first,  lr_v.first);
        out.R.vel = Vecd(lr_u.second, lr_v.second);
#elif SPH_NDIM == 3
        auto lr_w = reconstruct_component(Pi.vel[2], Pj.vel[2], grad_w_i, grad_w_j);
        out.L.vel = Vecd(lr_u.first,  lr_v.first,  lr_w.first);
        out.R.vel = Vecd(lr_u.second, lr_v.second, lr_w.second);
#endif
    }

    // pressure
    {
        auto lr = reconstruct_component(Pi.p, Pj.p, grad_p_i, grad_p_j);
        out.L.p = lr.first;  out.R.p = lr.second;
    }

    // positivity safeguard
    if (cfg.positivity) {
        out.L.rho = std::max(out.L.rho, cfg.small);
        out.R.rho = std::max(out.R.rho, cfg.small);
        out.L.p   = std::max(out.L.p,   cfg.small);
        out.R.p   = std::max(out.R.p,   cfg.small);
    }

    // total energy reconstructed from (rho, v, p) via EOS
    {
        const Real inv_gamma_minus_one = 1.0 / (cfg.gamma - 1.0);
        const Real kinL = (Real)0.5 * out.L.vel.dot(out.L.vel);
        const Real kinR = (Real)0.5 * out.R.vel.dot(out.R.vel);
        const Real eL   = out.L.p * inv_gamma_minus_one / out.L.rho;
        const Real eR   = out.R.p * inv_gamma_minus_one / out.R.rho;
        out.L.E = eL + kinL;
        out.R.E = eR + kinR;
    }

    return out;
}

} // namespace fluid_dynamics
} // namespace SPH
