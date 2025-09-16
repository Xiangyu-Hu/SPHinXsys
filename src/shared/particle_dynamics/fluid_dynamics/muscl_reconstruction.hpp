#pragma once
 
 #include "sphinxsys.h"
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
    SlopeLimiter limiter      = SlopeLimiter::Minmod;
    Real         theta        = 2.0;     
    bool         positivity   = true;    
    Real         small        = 1e-12;   
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

template <typename Vec>
inline std::pair<Real, Real>
reconstruct_scalar_muscl(Real Ui, const Vec& gradUi,
                         Real Uj, const Vec& gradUj,
                         const Vec& xi, const Vec& xj,
                         const Vec& x_iface, const Vec& nhat,
                         const SecondOrderConfig& cfg)
{
    const Vec di_vec = x_iface - xi;
    const Vec dj_vec = x_iface - xj;

    if (cfg.limiter == SlopeLimiter::None) {
        Real UL = gradUi.dot(di_vec) + Ui;
        Real UR = gradUj.dot(dj_vec) + Uj;
        return {UL, UR};
    }

    const Vec dx_vec = xj - xi;
    const Real du    = Uj - Ui;

    const Real si = gradUi.dot(dx_vec);
    const Real sj = gradUj.dot(-dx_vec);

    auto safe_div = [&](Real num, Real den){ return (std::abs(den) > 1e-14) ? (num/den) : 0.0; };

    const Real phi_i_raw = apply_limiter(cfg.limiter, si, du);
    const Real phi_j_raw = apply_limiter(cfg.limiter, sj, Ui - Uj);
 
    const Real phi_i = safe_div(phi_i_raw, (std::abs(si) > 1e-14 ? si : (Real)1)); // ∈[0,1] civarı
    const Real phi_j = safe_div(phi_j_raw, (std::abs(sj) > 1e-14 ? sj : (Real)1));
 
    Real UL = Ui + phi_i * gradUi.dot(di_vec);
    Real UR = Uj + phi_j * gradUj.dot(dj_vec);
 
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
                                       const Vec& x_iface, const Vec& nhat,
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
                                       const Vec& x_iface, const Vec& nhat,
                                       const SecondOrderConfig& cfg)
#endif
{
    LR out;

    // density
    {
        auto lr = reconstruct_scalar_muscl(Pi.rho, grad_rho_i,
                                           Pj.rho, grad_rho_j,
                                           xi, xj, x_iface, nhat, cfg);
        out.L.rho = lr.first;  out.R.rho = lr.second;
    }

    // velocity components
    {
        // u
        auto lr_u = reconstruct_scalar_muscl(Pi.vel[0], grad_u_i,
                                             Pj.vel[0], grad_u_j,
                                             xi, xj, x_iface, nhat, cfg);
        // v
        auto lr_v = reconstruct_scalar_muscl(Pi.vel[1], grad_v_i,
                                             Pj.vel[1], grad_v_j,
                                             xi, xj, x_iface, nhat, cfg);

#if SPH_NDIM == 2
        out.L.vel = Vecd(lr_u.first,  lr_v.first);
        out.R.vel = Vecd(lr_u.second, lr_v.second);
#elif SPH_NDIM == 3
        auto lr_w = reconstruct_scalar_muscl(Pi.vel[2], grad_w_i,
                                             Pj.vel[2], grad_w_j,
                                             xi, xj, x_iface, nhat, cfg);
        out.L.vel = Vecd(lr_u.first,  lr_v.first,  lr_w.first);
        out.R.vel = Vecd(lr_u.second, lr_v.second, lr_w.second);
#endif
    }

    // pressure
    {
        auto lr = reconstruct_scalar_muscl(Pi.p, grad_p_i,
                                           Pj.p, grad_p_j,
                                           xi, xj, x_iface, nhat, cfg);
        out.L.p = lr.first;  out.R.p = lr.second;
    }

    // total energy 
    out.L.E = Pi.E;
    out.R.E = Pj.E;

    // positivity safeguard
    if (cfg.positivity) {
        out.L.rho = std::max(out.L.rho, cfg.small);
        out.R.rho = std::max(out.R.rho, cfg.small);
        out.L.p   = std::max(out.L.p,   cfg.small);
        out.R.p   = std::max(out.R.p,   cfg.small);
    }

    return out;
}

} // namespace fluid_dynamics
} // namespace SPH
