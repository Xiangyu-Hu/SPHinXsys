#include <gtest/gtest.h>
#include "sphinxsys.h"
#include "muscl_reconstruction.hpp"
#include "muscl_hllc_integration.h"
#include "eulerian_fluid_dynamics/eulerian_riemann_solver.h"
#include "materials/compressible_fluid.h"

using namespace SPH;
using namespace SPH::fluid_dynamics;

namespace {

// Helpers
#if SPH_NDIM == 2
inline Vecd V(Real u, Real v) { return Vecd(u, v); }
#elif SPH_NDIM == 3
inline Vecd V(Real u, Real v, Real w) { return Vecd(u, v, w); }
#endif

struct StateStorage
{
    Real rho;
    Vecd vel;
    Real p;
    Real E;
    CompressibleFluidState state;

    StateStorage(Real rho_in, const Vecd &vel_in, Real p_in, Real gamma,
                 bool energy_is_density = true)
        : rho(rho_in),
          vel(vel_in),
          p(p_in),
          E(energy_is_density ? p_in / (gamma - 1.0) + 0.5 * rho_in * vel_in.squaredNorm()
                               : p_in / ((gamma - 1.0) * rho_in) + 0.5 * vel_in.squaredNorm()),
          state(rho, vel, p, E)
    {
    }

    StateStorage(const StateStorage &) = delete;
    StateStorage &operator=(const StateStorage &) = delete;
};

struct BridgeFixture : public ::testing::Test {
    // Fluid materials (left/right)
    CompressibleFluid *fluidL_;
    CompressibleFluid *fluidR_;
    Real gamma_ = 1.4;
    Real rho0_ = 1.0; // reference density

    BridgeFixture()
    : fluidL_(nullptr), fluidR_(nullptr)
    {}

    void SetUp() override {
        // Material initialization with proper constructor
        fluidL_ = new CompressibleFluid(rho0_, gamma_);
        fluidR_ = new CompressibleFluid(rho0_, gamma_);
    }

    void TearDown() override {
        delete fluidL_;
        delete fluidR_;
    }

    MUSCL_HLLC_Bridge make_bridge(bool use_hllc_diss_limiter=false, Real limiter_param=0.0) {
    MUSCLHLLCBridgeConfig cfg;
    cfg.gamma = gamma_;
    cfg.muscl_cfg.limiter = SlopeLimiter::None; 
    cfg.muscl_cfg.positivity = true;
    cfg.muscl_cfg.small = 1e-12;
    cfg.use_hllc_dissipation_limiter = use_hllc_diss_limiter;
    cfg.hllc_limiter_parameter = limiter_param;
        return MUSCL_HLLC_Bridge(*fluidL_, *fluidR_, cfg);
    }

    void zero_grad(Vecd &gr) const { gr.setZero(); }
};

} // namespace

// ============ TESTS ============

// 1) Uniform state invariance: Pi == Pj
TEST_F(BridgeFixture, UniformState_NoChange) {
    auto bridge = make_bridge(false, 0.0);

#if SPH_NDIM == 2
    Vecd xi(0.0, 0.0), xj(1.0, 0.0), xf(0.5, 0.0), n(1.0, 0.0);
    Vecd vel = V(50.0, 0.0);
#elif SPH_NDIM == 3
    Vecd xi(0.0, 0.0, 0.0), xj(1.0, 0.0, 0.0), xf(0.5, 0.0, 0.0), n(1.0, 0.0, 0.0);
    Vecd vel = V(50.0, 0.0, 0.0);
#endif

    Real rho = 1.2, p = 1.0e5;
    StateStorage left(rho, vel, p, gamma_);
    StateStorage right(rho, vel, p, gamma_);
    CompressibleFluidState &Pi = left.state;
    CompressibleFluidState &Pj = right.state;

    Vecd gr_rho_i, gr_rho_j, gr_u_i, gr_u_j, gr_v_i, gr_v_j, gr_p_i, gr_p_j;
    zero_grad(gr_rho_i); zero_grad(gr_rho_j);
    zero_grad(gr_u_i);   zero_grad(gr_u_j);
    zero_grad(gr_v_i);   zero_grad(gr_v_j);
    zero_grad(gr_p_i);   zero_grad(gr_p_j);

#if SPH_NDIM == 2
    auto star = bridge.getInterfaceState(Pi, Pj, xi, xj, xf, n,
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_p_i, gr_p_j);
#elif SPH_NDIM == 3
    Vecd gr_w_i, gr_w_j; zero_grad(gr_w_i); zero_grad(gr_w_j);
    auto star = bridge.getInterfaceState(Pi, Pj, xi, xj, xf, n,
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_w_i, gr_w_j,
                                         gr_p_i, gr_p_j);
#endif

    std::cout << "UniformState Debug:" << std::endl;
    std::cout << "  Input: rho=" << rho << ", p=" << p << ", vel=" << vel.transpose() << std::endl;
    std::cout << "  Output: rho=" << star.rho_ << ", p=" << star.p_ << ", vel=" << star.vel_.transpose() << std::endl;
    
    EXPECT_NEAR(star.rho_, rho, 1e-10);
    EXPECT_NEAR((star.vel_ - vel).norm(), 0.0, 1e-10);
    EXPECT_NEAR(star.p_, p, 1e-6);
}

// 2) Linear reproduction (smooth)
TEST_F(BridgeFixture, LinearReproduction_SmoothProfile) {
    auto bridge = make_bridge(false, 0.0);

#if SPH_NDIM == 2
    Vecd xi(0.0, 0.0), xj(1.0, 0.0), xf(0.5, 0.0), n(1.0, 0.0);
#elif SPH_NDIM == 3
    Vecd xi(0.0, 0.0, 0.0), xj(1.0, 0.0, 0.0), xf(0.5, 0.0, 0.0), n(1.0, 0.0, 0.0);
#endif

    // ρ(x)=ρ0 + a x, u(x)=u0 + b x, p(x)=p0 + c x
    Real rho0=1.0, a=0.2;
    Real u0=10.0, b=1.0;
    Real v0=0.0,  bv=0.5;
    Real p0=1.0e5, c=3.0e4;

#if SPH_NDIM == 2
    StateStorage left(rho0 + a * xi[0], V(u0 + b * xi[0], v0 + bv * xi[0]), p0 + c * xi[0], gamma_);
    StateStorage right(rho0 + a * xj[0], V(u0 + b * xj[0], v0 + bv * xj[0]), p0 + c * xj[0], gamma_);
#elif SPH_NDIM == 3
    StateStorage left(rho0 + a * xi[0], V(u0 + b * xi[0], v0 + bv * xi[0], 0.0), p0 + c * xi[0], gamma_);
    StateStorage right(rho0 + a * xj[0], V(u0 + b * xj[0], v0 + bv * xj[0], 0.0), p0 + c * xj[0], gamma_);
#endif
    CompressibleFluidState &Pi = left.state;
    CompressibleFluidState &Pj = right.state;

 
    Vecd gr_rho_i=Vecd::Zero(), gr_rho_j=Vecd::Zero();
    Vecd gr_u_i=Vecd::Zero(),   gr_u_j=Vecd::Zero();
    Vecd gr_v_i=Vecd::Zero(),   gr_v_j=Vecd::Zero();
    Vecd gr_p_i=Vecd::Zero(),   gr_p_j=Vecd::Zero();
    gr_rho_i[0]=gr_rho_j[0]=a;
    gr_u_i[0]=gr_u_j[0]=b;
    gr_v_i[0]=gr_v_j[0]=bv;
    gr_p_i[0]=gr_p_j[0]=c;

#if SPH_NDIM == 2
    auto star = bridge.getInterfaceState(Pi, Pj, xi, xj, xf, n,
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_p_i, gr_p_j);
#elif SPH_NDIM == 3
    Vecd gr_w_i=Vecd::Zero(), gr_w_j=Vecd::Zero();
    auto star = bridge.getInterfaceState(Pi, Pj, xi, xj, xf, n,
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_w_i, gr_w_j,
                                         gr_p_i, gr_p_j);
#endif

    // Reference: MUSCL reconstruction of L/R and face value
    Primitives Pr_i{Pi.rho_, Pi.vel_, Pi.p_, Pi.E_};
    Primitives Pr_j{Pj.rho_, Pj.vel_, Pj.p_, Pj.E_};
    SecondOrderConfig cfg_ref; cfg_ref.limiter = SlopeLimiter::None; cfg_ref.positivity = false;
    auto lr = reconstruct_primitives_muscl(
        Pr_i, Pr_j,
        gr_rho_i, gr_rho_j, gr_u_i, gr_u_j, gr_v_i, gr_v_j, gr_p_i, gr_p_j,
        xi, xj, xf, n, cfg_ref
    );

    auto within = [](Real x, Real a, Real b){ Real lo=std::min(a,b)-1e-10, hi=std::max(a,b)+1e-10; return x>=lo && x<=hi; };

    // 1) MUSCL linear reproduction (face within L/R)
    Real rho_face = rho0 + a*xf[0];
    Real u_face   = u0  + b*xf[0];
    Real v_face   = v0  + bv*xf[0];
    Real p_face   = p0  + c*xf[0];
    EXPECT_TRUE(within(rho_face, lr.L.rho, lr.R.rho));
    EXPECT_TRUE(within(u_face,   lr.L.vel[0], lr.R.vel[0]));
    EXPECT_TRUE(within(v_face,   lr.L.vel[1], lr.R.vel[1]));
    EXPECT_TRUE(within(p_face,   lr.L.p,      lr.R.p));

    // 2) HLLC star-state controls
    // 2a) rho and p boundedness (no overshoot)
    EXPECT_TRUE(within(star.rho_,    lr.L.rho,    lr.R.rho));
    EXPECT_TRUE(within(star.p_,      lr.L.p,      lr.R.p));

    // 2b) normal velocity bounded by wave speeds: -n·v* ∈ [s_l, s_r]
    auto sound = [&](Real rho, Real p){ return std::sqrt(gamma_ * p / std::max(rho, (Real)1e-12)); };
    Real ul = -n.dot(lr.L.vel);
    Real ur = -n.dot(lr.R.vel);
    Real aL = sound(lr.L.rho, lr.L.p);
    Real aR = sound(lr.R.rho, lr.R.p);
    Real s_l = std::min(ul - aL, ur - aR);
    Real s_r = std::max(ul + aL, ur + aR);
    Real un_star = -n.dot(star.vel_);
    EXPECT_LE(un_star, s_r + 1e-8);
    EXPECT_GE(un_star, s_l - 1e-8);

    // 2c) Tangential velocity between L and R
    EXPECT_TRUE(within(star.vel_[1], lr.L.vel[1], lr.R.vel[1]));
}

// 3) Sod shock (1D embedded)
TEST_F(BridgeFixture, SodShock_1DNormal) {
    auto bridge = make_bridge(true, 0.5);

#if SPH_NDIM == 2
    Vecd xi(0.0, 0.0), xj(1.0, 0.0), xf(0.5, 0.0), n(1.0, 0.0);
    Vecd zeroV = V(0.0, 0.0);
#elif SPH_NDIM == 3
    Vecd xi(0.0, 0.0, 0.0), xj(1.0, 0.0, 0.0), xf(0.5, 0.0, 0.0), n(1.0, 0.0, 0.0);
    Vecd zeroV = V(0.0, 0.0, 0.0);
#endif

    StateStorage left(1.0, zeroV, 1.0e5, gamma_);
    StateStorage right(0.125, zeroV, 1.0e4, gamma_);
    CompressibleFluidState &PL = left.state;
    CompressibleFluidState &PR = right.state;

    // Add some gradients for Sod shock test
    Vecd gr_rho_i=Vecd::Zero(), gr_rho_j=Vecd::Zero();
    Vecd gr_u_i=Vecd::Zero(),   gr_u_j=Vecd::Zero();
    Vecd gr_v_i=Vecd::Zero(),   gr_v_j=Vecd::Zero();
    Vecd gr_p_i=Vecd::Zero(),   gr_p_j=Vecd::Zero();
    
    // Small gradients to create interface difference
    gr_rho_i[0] = -0.1; gr_rho_j[0] = 0.1;
    gr_u_i[0] = 0.0;    gr_u_j[0] = 0.0;
    gr_v_i[0] = 0.0;    gr_v_j[0] = 0.0;
    gr_p_i[0] = -1000;  gr_p_j[0] = 1000;
#if SPH_NDIM == 2
    auto star = bridge.getInterfaceState(PL, PR, xi, xj, xf, n, 
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_p_i, gr_p_j);
#elif SPH_NDIM == 3
    Vecd gr_w_i=Vecd::Zero(), gr_w_j=Vecd::Zero();
    auto star = bridge.getInterfaceState(PL, PR, xi, xj, xf, n, 
                                         gr_rho_i, gr_rho_j,
                                         gr_u_i, gr_u_j,
                                         gr_v_i, gr_v_j,
                                         gr_w_i, gr_w_j,
                                         gr_p_i, gr_p_j);
#endif

    EXPECT_GT(star.rho_, 0.0);
    EXPECT_GT(star.p_,   0.0);

    EXPECT_LT(std::abs(star.vel_[0]), 1e3);
    EXPECT_LT(std::abs(star.vel_[1]), 1e3);
}

// 4) Positivity guard
TEST_F(BridgeFixture, PositivityGuard_NoNegative) {
    auto bridge = make_bridge(false, 0.0);

#if SPH_NDIM == 2
    Vecd xi(0.0, 0.0), xj(1.0, 0.0), xf(0.5, 0.0), n(1.0, 0.0);
    Vecd uL = V(-300.0, 0.0), uR = V(300.0, 0.0);
#elif SPH_NDIM == 3
    Vecd xi(0.0, 0.0, 0.0), xj(1.0, 0.0, 0.0), xf(0.5, 0.0, 0.0), n(1.0, 0.0, 0.0);
    Vecd uL = V(-300.0, 0.0, 0.0), uR = V(300.0, 0.0, 0.0);
#endif

    StateStorage left(1.0, uL, 2.0e4, gamma_);
    StateStorage right(1.0, uR, 2.0e4, gamma_);
    CompressibleFluidState &PL = left.state;
    CompressibleFluidState &PR = right.state;

    Vecd z=Vecd::Zero();
#if SPH_NDIM == 2
    auto star = bridge.getInterfaceState(PL, PR, xi, xj, xf, n, z,z, z,z, z,z, z,z);
#elif SPH_NDIM == 3
    auto star = bridge.getInterfaceState(PL, PR, xi, xj, xf, n, z,z, z,z, z,z, z,z, z,z);
#endif

    EXPECT_GT(star.rho_, 0.0);
    EXPECT_GT(star.p_,   0.0);
}
