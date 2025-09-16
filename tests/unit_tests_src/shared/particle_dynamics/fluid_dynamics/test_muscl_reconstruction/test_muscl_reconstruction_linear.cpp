// test_muscl_reconstruction_linear.cpp
#include "gtest/gtest.h"
#include "muscl_reconstruction.hpp"
#include "sphinxsys.h"

using namespace SPH;
using namespace SPH::fluid_dynamics;

TEST(MUSCL, LinearReproductionScalar)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::None;
    cfg.positivity = false;

    // 2D example: U(x,y) = a + b·x
    const Real a  = 1.0;
    const Vecd b  = Vecd(0.3, -0.2);

    auto U = [&](const Vecd& x){ return a + b.dot(x); };

    const Vecd xi(0.0, 0.0);
    const Vecd xj(1.0, 0.5);
    const Vecd x_iface(0.6, 0.2);
    Vecd nhat = (xj - xi);
    nhat /= nhat.norm();

    const Real Ui = U(xi), Uj = U(xj);

    // In a linear field, grad is constant: ∇U = b
    auto lr = reconstruct_scalar_muscl(
        Ui, b, Uj, b, xi, xj, x_iface, nhat, cfg
    );

    const Real Uface = U(x_iface);
    EXPECT_NEAR(lr.first,  Uface, 1e-12);
    EXPECT_NEAR(lr.second, Uface, 1e-12);
}

TEST(MUSCL, LinearReproductionPrimitives)
{
    SecondOrderConfig cfg;
    cfg.limiter    = SlopeLimiter::None;
    cfg.positivity = false;

    // rho(x,y), u(x,y), v(x,y), p(x,y) are linear
    const Vecd xi(0.1, -0.2), xj(1.1, 0.7), x_iface(0.4, 0.1);
    Vecd nhat = (xj - xi); nhat /= nhat.norm();

    auto Lf = [&](const Real a, const Vecd& b, const Vecd& x){ return a + b.dot(x); };

    // coefficients
    const Real ar=1.2, ap=4.0;  Vecd br(0.1, -0.05), bp(0.2, 0.15);
    const Real au=0.3, av=-0.7; Vecd bu(0.05, 0.01),  bv(-0.02, 0.04);

    Primitives Pi, Pj;
    Pi.rho = Lf(ar, br, xi);   Pj.rho = Lf(ar, br, xj);
    Pi.vel = Vecd(Lf(au, bu, xi), Lf(av, bv, xi));
    Pj.vel = Vecd(Lf(au, bu, xj), Lf(av, bv, xj));
    Pi.p   = Lf(ap, bp, xi);   Pj.p   = Lf(ap, bp, xj);
    Pi.E = 0.0; Pj.E = 0.0; // energy not tested here

    const Vecd gr_i = br, gr_j = br;
    const Vecd gu_i = bu, gu_j = bu;
    const Vecd gv_i = bv, gv_j = bv;
    const Vecd gp_i = bp, gp_j = bp;

    auto lr = reconstruct_primitives_muscl(
        Pi, Pj, gr_i, gr_j, gu_i, gu_j, gv_i, gv_j, gp_i, gp_j,
        xi, xj, x_iface, nhat, cfg
    );

    const Real rho_face = Lf(ar, br, x_iface);
    const Real  u_face  = Lf(au, bu, x_iface);
    const Real  v_face  = Lf(av, bv, x_iface);
    const Real  p_face  = Lf(ap, bp, x_iface);

    EXPECT_NEAR(lr.L.rho, rho_face, 1e-12);
    EXPECT_NEAR(lr.R.rho, rho_face, 1e-12);
    EXPECT_NEAR(lr.L.vel[0], u_face, 1e-12);
    EXPECT_NEAR(lr.R.vel[0], u_face, 1e-12);
    EXPECT_NEAR(lr.L.vel[1], v_face, 1e-12);
    EXPECT_NEAR(lr.R.vel[1], v_face, 1e-12);
    EXPECT_NEAR(lr.L.p,   p_face, 1e-12);
    EXPECT_NEAR(lr.R.p,   p_face, 1e-12);
}
