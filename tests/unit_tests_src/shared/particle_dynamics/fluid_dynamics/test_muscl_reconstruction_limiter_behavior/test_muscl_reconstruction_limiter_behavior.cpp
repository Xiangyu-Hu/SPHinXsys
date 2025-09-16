// test_muscl_reconstruction_limiter_behavior.cpp
#include "gtest/gtest.h"
#include "muscl_reconstruction.hpp"
#include "sphinxsys.h"

using namespace SPH;
using namespace SPH::fluid_dynamics;

#if SPH_NDIM == 2
static inline Vecd V(Real x, Real y) { return Vecd(x, y); }
static inline Vecd Z()               { return Vecd(0.0, 0.0); }
#elif SPH_NDIM == 3
static inline Vecd V(Real x, Real y, Real z=0.0) { return Vecd(x, y, z); }
static inline Vecd Z()                            { return Vecd(0.0, 0.0, 0.0); }
#else
#error "SPH_NDIM must be 2 or 3"
#endif

// Helpers
static inline Vecd unit(const Vecd& a)
{
    Real n = a.norm();
    return (n > (Real)0) ? a / n : a;
}

static inline Real dot(const Vecd& a, const Vecd& b)
{
    return a.dot(b);
}

// Linear field: U(x) = a + b·x
static inline Real U_lin(Real a, const Vecd& b, const Vecd& x)
{
    return a + dot(b, x);
}

// 1) Linear reproduction: limiter OFF
TEST(MUSCL_LimiterSuite, LinearReproductionScalar_Unlimited)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::None;
    cfg.positivity = false;

    const Real a = 1.0;
    const Vecd b = V(0.3, -0.2);

    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.5);
    const Vecd xf = V(0.6, 0.2);
    const Vecd nh = unit(xj - xi);

    const Real Ui = U_lin(a, b, xi);
    const Real Uj = U_lin(a, b, xj);

    auto lr = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg);
    const Real Uface = U_lin(a, b, xf);

    EXPECT_NEAR(lr.first,  Uface, 1e-12);
    EXPECT_NEAR(lr.second, Uface, 1e-12);
}

// 2) Linear reproduction: limiter ON
TEST(MUSCL_LimiterSuite, LinearReproductionScalar_Minmod)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::Minmod;
    cfg.positivity = false;

    const Real a = 1.0;
    const Vecd b = V(0.3, -0.2);

    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.5);
    const Vecd xf = V(0.6, 0.2);
    const Vecd nh = unit(xj - xi);

    const Real Ui = U_lin(a, b, xi);
    const Real Uj = U_lin(a, b, xj);

    auto lr = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg);
    const Real Uface = U_lin(a, b, xf);

    EXPECT_NEAR(lr.first,  Uface, 1e-12);
    EXPECT_NEAR(lr.second, Uface, 1e-12);
}

TEST(MUSCL_LimiterSuite, LinearReproductionScalar_MC)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::MC;
    cfg.positivity = false;

    const Real a = 1.0;
    const Vecd b = V(0.3, -0.2);

    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.5);
    const Vecd xf = V(0.6, 0.2);
    const Vecd nh = unit(xj - xi);

    const Real Ui = U_lin(a, b, xi);
    const Real Uj = U_lin(a, b, xj);

    auto lr = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg);
    const Real Uface = U_lin(a, b, xf);

    EXPECT_NEAR(lr.first,  Uface, 1e-12);
    EXPECT_NEAR(lr.second, Uface, 1e-12);
}

TEST(MUSCL_LimiterSuite, LinearReproductionScalar_VanLeer)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::VanLeer;
    cfg.positivity = false;

    const Real a = 1.0;
    const Vecd b = V(0.3, -0.2);

    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.5);
    const Vecd xf = V(0.6, 0.2);
    const Vecd nh = unit(xj - xi);

    const Real Ui = U_lin(a, b, xi);
    const Real Uj = U_lin(a, b, xj);

    auto lr = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg);
    const Real Uface = U_lin(a, b, xf);

    EXPECT_NEAR(lr.first,  Uface, 1e-12);
    EXPECT_NEAR(lr.second, Uface, 1e-12);
}

// 3) Discontinuity: limiter OFF vs ON (overshoot control)
TEST(MUSCL_LimiterSuite, DiscontinuityScalar_OvershootControl_Minmod)
{
    // Synthetic profile with opposite gradients; limiter should suppress overshoot.
    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.0);
    const Vecd xf = V(0.5, 0.0);
    const Vecd nh = unit(xj - xi);

    const Real Ui = 10.0;
    const Real Uj = 1.0;

    // Nonzero gradients can overshoot without limiter
    const Vecd gi = V(6.0, 0.0);
    const Vecd gj = V(-6.0, 0.0);

    // Unlimited (reference)
    SecondOrderConfig cfg_unlim; cfg_unlim.limiter = SlopeLimiter::None; cfg_unlim.positivity = false;
    auto lr_unlim = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg_unlim);

    // Minmod
    SecondOrderConfig cfg_mm; cfg_mm.limiter = SlopeLimiter::Minmod; cfg_mm.positivity = false;
    auto lr_mm = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg_mm);

    // TVD/monotonicity: limited results must stay within [Ui, Uj]
    Real Umin = std::min(Ui, Uj), Umax = std::max(Ui, Uj);
    EXPECT_GE(lr_mm.first,  Umin);
    EXPECT_LE(lr_mm.first,  Umax);
    EXPECT_GE(lr_mm.second, Umin);
    EXPECT_LE(lr_mm.second, Umax);

    // Slope shrinking: limited slope should not exceed unlimited
    EXPECT_LE(std::abs(lr_mm.first  - Ui), std::abs(lr_unlim.first  - Ui));
    EXPECT_LE(std::abs(lr_mm.second - Uj), std::abs(lr_unlim.second - Uj));
}

TEST(MUSCL_LimiterSuite, DiscontinuityScalar_OvershootControl_MC_VanLeer)
{
    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.0);
    const Vecd xf = V(0.5, 0.0);
    const Vecd nh = unit(xj - xi);

    const Real Ui = 10.0, Uj = 1.0;
    const Vecd gi = V(6.0, 0.0), gj = V(-6.0, 0.0);

    SecondOrderConfig cfg_unlim; cfg_unlim.limiter = SlopeLimiter::None; cfg_unlim.positivity = false;
    auto lr_unlim = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg_unlim);

    // MC
    {
        SecondOrderConfig cfg; cfg.limiter = SlopeLimiter::MC; cfg.positivity = false;
        auto lr = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg);
        Real Umin = std::min(Ui, Uj), Umax = std::max(Ui, Uj);
        EXPECT_GE(lr.first,  Umin); EXPECT_LE(lr.first,  Umax);
        EXPECT_GE(lr.second, Umin); EXPECT_LE(lr.second, Umax);
        EXPECT_LE(std::abs(lr.second - lr.first), std::abs(Uj - Ui) + 1e-12);
        EXPECT_LE(std::abs(lr.second - lr.first), std::abs(lr_unlim.second - lr_unlim.first) + 1e-12);
    }
    // VanLeer
    {
        SecondOrderConfig cfg; cfg.limiter = SlopeLimiter::VanLeer; cfg.positivity = false;
        auto lr = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg);
        Real Umin = std::min(Ui, Uj), Umax = std::max(Ui, Uj);
        EXPECT_GE(lr.first,  Umin); EXPECT_LE(lr.first,  Umax);
        EXPECT_GE(lr.second, Umin); EXPECT_LE(lr.second, Umax);
    }
}

// ---------- 4) Monoton smooth profile----------
TEST(MUSCL_LimiterSuite, MonotoneSmooth_Profile_NoExcessCut)
{
    const Real a = 2.0;
#if SPH_NDIM == 2
    const Vecd b = V(0.7, 0.2);
    const Vecd xi = V(-0.2, 0.4);
    const Vecd xj = V( 1.1, 0.7);
    const Vecd xf = V( 0.3, 0.6);
#else
    const Vecd b = V(0.7, 0.2, -0.1);
    const Vecd xi = V(-0.2, 0.4, 0.0);
    const Vecd xj = V( 1.1, 0.7, 0.2);
    const Vecd xf = V( 0.3, 0.6, 0.1);
#endif
    const Vecd nh = unit(xj - xi);

    const Real Ui = U_lin(a, b, xi);
    const Real Uj = U_lin(a, b, xj);
    const Real Uface = U_lin(a, b, xf);

    // Unlimited referans
    SecondOrderConfig cfg_ref; cfg_ref.limiter = SlopeLimiter::None; cfg_ref.positivity = false;
    auto lr_ref = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg_ref);

  
    for (auto lim : {SlopeLimiter::Minmod, SlopeLimiter::MC, SlopeLimiter::VanLeer})
    {
        SecondOrderConfig cfg; cfg.limiter = lim; cfg.positivity = false;
        auto lr = reconstruct_scalar_muscl(Ui, b, Uj, b, xi, xj, xf, nh, cfg);
        EXPECT_NEAR(lr.first,  Uface, 1e-12);
        EXPECT_NEAR(lr.second, Uface, 1e-12);
        EXPECT_NEAR(lr.first,  lr_ref.first,  1e-12);
        EXPECT_NEAR(lr.second, lr_ref.second, 1e-12);
    }
}

// ---------- 5) Primitives: lineer reproduction (limiter is open) ----------
TEST(MUSCL_LimiterSuite, PrimitivesLinear_WithLimiter)
{
    SecondOrderConfig cfg; cfg.limiter = SlopeLimiter::MC; cfg.positivity = false;

    // rho, u, v, p all are linear
    const Vecd xi = V(0.1, -0.2);
    const Vecd xj = V(1.1, 0.7);
    const Vecd xf = V(0.4, 0.1);
    const Vecd nh = unit(xj - xi);

    //  coefficients
    const Real ar = 1.2, ap = 4.0;      Vecd br = V(0.1, -0.05), bp = V(0.2, 0.15);
    const Real au = 0.3, av = -0.7;     Vecd bu = V(0.05, 0.01),  bv = V(-0.02, 0.04);

    auto Lf = [&](Real a0, const Vecd& b0, const Vecd& x){ return a0 + dot(b0, x); };

    Primitives Pi{}, Pj{};
    Pi.rho = Lf(ar, br, xi);   Pj.rho = Lf(ar, br, xj);
    Pi.vel = V(Lf(au, bu, xi), Lf(av, bv, xi));
    Pj.vel = V(Lf(au, bu, xj), Lf(av, bv, xj));
    Pi.p   = Lf(ap, bp, xi);   Pj.p   = Lf(ap, bp, xj);
    Pi.E = 0.0; Pj.E = 0.0;

    const Vecd gr_i = br, gr_j = br;
    const Vecd gu_i = bu, gu_j = bu;
    const Vecd gv_i = bv, gv_j = bv;
    const Vecd gp_i = bp, gp_j = bp;

    auto lr = reconstruct_primitives_muscl(
        Pi, Pj, gr_i, gr_j, gu_i, gu_j, gv_i, gv_j, gp_i, gp_j,
        xi, xj, xf, nh, cfg
    );

    const Real rho_f = Lf(ar, br, xf);
    const Real  u_f  = Lf(au, bu, xf);
    const Real  v_f  = Lf(av, bv, xf);
    const Real  p_f  = Lf(ap, bp, xf);

    EXPECT_NEAR(lr.L.rho,   rho_f, 1e-12);
    EXPECT_NEAR(lr.R.rho,   rho_f, 1e-12);
    EXPECT_NEAR(lr.L.vel[0], u_f,  1e-12);
    EXPECT_NEAR(lr.R.vel[0], u_f,  1e-12);
    EXPECT_NEAR(lr.L.vel[1], v_f,  1e-12);
    EXPECT_NEAR(lr.R.vel[1], v_f,  1e-12);
    EXPECT_NEAR(lr.L.p,      p_f,  1e-12);
    EXPECT_NEAR(lr.R.p,      p_f,  1e-12);
}

// ---------- 6) posivity guard ----------
TEST(MUSCL_LimiterSuite, PositivityClipping_RhoAndP)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::Minmod;
    cfg.positivity = true;
    cfg.small      = 1e-10;

    Primitives Pi{}, Pj{};
    Pi.rho = 1e-12;  Pj.rho = -5.0;          
    Pi.vel = Z();    Pj.vel = Z();
    Pi.p   = 1e-13;  Pj.p   = -2.0;          
    Pi.E   = 0.0;    Pj.E   = 0.0;

    const Vecd g0 = Z();
    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.0);
    const Vecd xf = V(0.5, 0.0);
    const Vecd nh = unit(xj - xi);

    LR lr = reconstruct_primitives_muscl(
        Pi, Pj, g0, g0, g0, g0, g0, g0, g0, g0, xi, xj, xf, nh, cfg
    );

    EXPECT_GE(lr.L.rho, cfg.small);
    EXPECT_GE(lr.R.rho, cfg.small);
    EXPECT_GE(lr.L.p,   cfg.small);
    EXPECT_GE(lr.R.p,   cfg.small);

    EXPECT_DOUBLE_EQ(lr.L.vel[0], 0.0);
    EXPECT_DOUBLE_EQ(lr.L.vel[1], 0.0);
    EXPECT_DOUBLE_EQ(lr.R.vel[0], 0.0);
    EXPECT_DOUBLE_EQ(lr.R.vel[1], 0.0);
}

// ---------- 7) Discontinuity----------
TEST(MUSCL_LimiterSuite, DiscontinuityScalar_ZeroGrad_DropsToFirstOrder)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::Minmod;
    cfg.positivity = false;

    const Real Ui = 10.0, Uj = 1.0;
    const Vecd gi = Z(),   gj = Z();

    const Vecd xi = V(0.0, 0.0);
    const Vecd xj = V(1.0, 0.0);
    const Vecd xf = V(0.5, 0.0);
    const Vecd nh = unit(xj - xi);

    auto lr = reconstruct_scalar_muscl(Ui, gi, Uj, gj, xi, xj, xf, nh, cfg);

    EXPECT_NEAR(lr.first,  Ui, 1e-12);
    EXPECT_NEAR(lr.second, Uj, 1e-12);

    Real Umin = std::min(Ui, Uj), Umax = std::max(Ui, Uj);
    EXPECT_GE(lr.first,  Umin); EXPECT_LE(lr.first,  Umax);
    EXPECT_GE(lr.second, Umin); EXPECT_LE(lr.second, Umax);
}
