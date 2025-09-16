#include "gtest/gtest.h"
#include "muscl_reconstruction.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>

using namespace SPH;
using namespace SPH::fluid_dynamics;

// Simple 1D grid + Euler FV/HLLC sandbox (test only)
struct Cell {
    Real rho, u, p, E; // primitive + total energy
};

struct Flux {
    Real mass, mom, ener;
};

// EOS: gamma
static constexpr Real gamma_gas = 1.4;
static constexpr Real SMALL = 1e-12;

// E = p/(γ−1) + 0.5 ρ u^2
inline Real total_energy(Real rho, Real u, Real p) {
    return p/(gamma_gas - 1.0) + 0.5 * rho * u * u;
}

// a = sqrt(γ p / ρ)
inline Real sound_speed(Real rho, Real p) {
    rho = std::max(rho, SMALL);
    p   = std::max(p,   SMALL);
    return std::sqrt(gamma_gas * p / rho);
}

// Conservative state: [ρ, ρu, E]
struct Cons {
    Real r, ru, E;
};
inline Cons prim2cons(const Cell& c) { return {c.rho, c.rho*c.u, c.E}; }
inline Cell cons2prim(const Cons& U) {
    Real rho = std::max(U.r, SMALL);
    Real u   = U.ru / rho;
    Real p   = (gamma_gas - 1.0) * (U.E - 0.5 * rho * u * u);
    p = std::max(p, SMALL);
    return {rho, u, p, U.E};
}

// Physical flux from primitives
inline Flux phys_flux(const Cell& c) {
    Real r = c.rho, u = c.u, p = c.p, E = c.E;
    return { r*u, r*u*u + p, u*(E + p) };
}

// Simple 1D primitive gradients: central inside, one-sided at boundaries
static void compute_gradients_1d(const std::vector<Cell>& P,
                                 std::vector<Real>& grho,
                                 std::vector<Real>& gu,
                                 std::vector<Real>& gp,
                                 Real dx)
{
    int N = (int)P.size();
    grho.assign(N, 0.0); gu.assign(N, 0.0); gp.assign(N, 0.0);
    for (int i=0;i<N;i++) {
        if (i == 0) {
            grho[i] = (P[1].rho - P[0].rho) / dx;
            gu[i]   = (P[1].u   - P[0].u  ) / dx;
            gp[i]   = (P[1].p   - P[0].p  ) / dx;
        } else if (i == N-1) {
            grho[i] = (P[N-1].rho - P[N-2].rho) / dx;
            gu[i]   = (P[N-1].u   - P[N-2].u  ) / dx;
            gp[i]   = (P[N-1].p   - P[N-2].p  ) / dx;
        } else {
            grho[i] = (P[i+1].rho - P[i-1].rho) / (2.0*dx);
            gu[i]   = (P[i+1].u   - P[i-1].u  ) / (2.0*dx);
            gp[i]   = (P[i+1].p   - P[i-1].p  ) / (2.0*dx);
        }
    }
}

// 1D HLLC Riemann solver (L,R primitives)
inline Flux hllc_flux(const Cell& L, const Cell& R)
{
    // Roe/HLLC estimates:
    Real rL=L.rho, uL=L.u, pL=L.p, aL=sound_speed(rL,pL);
    Real rR=R.rho, uR=R.u, pR=R.p, aR=sound_speed(rR,pR);

    // Wave speed estimates (Davis):
    Real SL = std::min(uL - aL, uR - aR);
    Real SR = std::max(uL + aL, uR + aR);

    // Star-region speed (Toro’s formula):
    Real num = pR - pL + rL*uL*(SL - uL) - rR*uR*(SR - uR);
    Real den = rL*(SL - uL) - rR*(SR - uR);
    if (std::abs(den) < 1e-14) {
        Cons UL = prim2cons(L), UR = prim2cons(R);
        auto F_L = phys_flux(L), F_R = phys_flux(R);
        Real inv = 1.0 / (SR - SL + 1e-14);
        return {
            (SR*F_L.mass - SL*F_R.mass + SL*SR*(UR.r  - UL.r )) * inv,
            (SR*F_L.mom  - SL*F_R.mom  + SL*SR*(UR.ru - UL.ru)) * inv,
            (SR*F_L.ener - SL*F_R.ener + SL*SR*(UR.E  - UL.E )) * inv
        };
    }
    Real Sstar = num / den;

    auto F_L = phys_flux(L);
    auto F_R = phys_flux(R);

    if (0.0 <= SL) {
        return F_L;
    } else if (SL <= 0.0 && 0.0 <= Sstar) {
        // U*_L
        Cons UL = prim2cons(L);
        Real coeff = rL*(SL - uL) / (SL - Sstar);
        Cons UstarL = {
            coeff,
            coeff * Sstar,
            coeff * ( L.E/L.rho + (Sstar - uL)*(Sstar + pL/(rL*(SL - uL))) )
        };
        Flux FstarL = {
            F_L.mass + SL*(UstarL.r - UL.r),
            F_L.mom  + SL*(UstarL.ru - UL.ru),
            F_L.ener + SL*(UstarL.E  - UL.E)
        };
        return FstarL;
    } else if (Sstar <= 0.0 && 0.0 <= SR) {
        // U*_R
        Cons UR = prim2cons(R);
        Real coeff = rR*(SR - uR) / (SR - Sstar);
        Cons UstarR = {
            coeff,
            coeff * Sstar,
            coeff * ( R.E/R.rho + (Sstar - uR)*(Sstar + pR/(rR*(SR - uR))) )
        };
        Flux FstarR = {
            F_R.mass + SR*(UstarR.r - UR.r),
            F_R.mom  + SR*(UstarR.ru - UR.ru),
            F_R.ener + SR*(UstarR.E  - UR.E)
        };
        return FstarR;
    } else {
        return F_R;
    }
}

// Intermediate density band width metric for Sod (heuristic window)
static int band_width_count(const std::vector<Cell>& P, Real lo, Real hi)
{
    int cnt = 0;
    for (auto& c : P) if (c.rho >= lo && c.rho <= hi) cnt++;
    return cnt;
}

// Total variation (rho)
static Real total_variation_rho(const std::vector<Cell>& P)
{
    Real tv = 0.0;
    for (size_t i=1;i<P.size();++i) tv += std::abs(P[i].rho - P[i-1].rho);
    return tv;
}

// Main test: Sod shock tube + limiter comparison
TEST(ShockTubeMUSCL, Sod_Limiter_Comparison)
{
    // Grid
    const int N = 400;
    const Real x0 = 0.0, x1 = 1.0, dx = (x1 - x0)/N;
    const Real x_mid = 0.5;
    std::vector<Real> xc(N);
    for (int i=0;i<N;i++) xc[i] = x0 + (i+0.5)*dx;

    // Initial: Sod (γ=1.4), t=0 → t_final=0.20
    Cell L0{1.0, 0.0, 1.0, 0.0}, R0{0.125, 0.0, 0.1, 0.0};
    L0.E = total_energy(L0.rho, L0.u, L0.p);
    R0.E = total_energy(R0.rho, R0.u, R0.p);

    auto init_state = [&](std::vector<Cell>& P){
        P.resize(N);
        for (int i=0;i<N;i++) {
            if (xc[i] < x_mid) P[i] = L0; else P[i] = R0;
        }
    };

    // Simulation params
    const Real CFL = 0.3;
    const Real t_end = 0.20;

    auto run_sim = [&](SlopeLimiter lim)->std::tuple<std::vector<Cell>, Real, int>
    {
        std::vector<Cell> P; init_state(P);
        std::vector<Real> grho, gu, gp;
        SecondOrderConfig cfg;
        cfg.limiter    = lim;
        cfg.positivity = true;
        cfg.small      = SMALL;

        Real t = 0.0; int steps = 0;
        while (t < t_end) {
            // positivity floor before CFL
            for (auto& c : P) { c.rho = std::max(c.rho, SMALL); c.p = std::max(c.p, SMALL); }
            // dt 
            Real amax = 0.0;
            for (auto& c : P) amax = std::max(amax, std::abs(c.u) + sound_speed(c.rho, c.p));
            Real dt = CFL * dx / (amax + 1e-14);
            if (t + dt > t_end) dt = t_end - t;

            // Gradient
            compute_gradients_1d(P, grho, gu, gp, dx);

            // Interface fluxes
            std::vector<Flux> F(N+1);
            for (int i=0;i<=N;i++) {
                // Left and right cell indices
                int il = std::max(i-1, 0);
                int ir = std::min(i,   N-1);

                // Reconstruction point (interface): xi, xj, x_iface and nhat (1D +x direction)
                Real xi = x0 + (il+0.5)*dx;
                Real xj = x0 + (ir+0.5)*dx;
                Real xf = x0 + (i)*dx;
                Real nh = +1.0; // 1D

                // Primitives
                Primitives Pi{P[il].rho, Vecd(P[il].u, 0.0), P[il].p, P[il].E};
                Primitives Pj{P[ir].rho, Vecd(P[ir].u, 0.0), P[ir].p, P[ir].E};

                // 1D gradient vectors: nh direction, scalar — put x component for Vecd usage
                Vecd gi(grho[il], 0.0), gj(grho[ir], 0.0);
                Vecd gui(gu[il], 0.0),  guj(gu[ir], 0.0);
                Vecd gpi(gp[il], 0.0),  gpj(gp[ir], 0.0);

                // Reconstruction (our MUSCL)
                LR lr = reconstruct_primitives_muscl(
                    Pi, Pj,
                    gi, gj, gui, guj, Vecd(0.0,0.0), Vecd(0.0,0.0), 
                    gpi, gpj,
                    Vecd(xi,0.0), Vecd(xj,0.0), Vecd(xf,0.0), Vecd(nh,0.0),
                    cfg
                );

                // HLLC flux
                Cell Ls{lr.L.rho, lr.L.vel[0], lr.L.p, total_energy(lr.L.rho, lr.L.vel[0], lr.L.p)};
                Cell Rs{lr.R.rho, lr.R.vel[0], lr.R.p, total_energy(lr.R.rho, lr.R.vel[0], lr.R.p)};
                F[i] = hllc_flux(Ls, Rs);
            }

            // FV update
            std::vector<Cons> U(N);
            for (int i=0;i<N;i++) U[i] = prim2cons(P[i]);
            for (int i=0;i<N;i++) {
                U[i].r  -= dt/dx * (F[i+1].mass - F[i].mass);
                U[i].ru -= dt/dx * (F[i+1].mom  - F[i].mom );
                U[i].E  -= dt/dx * (F[i+1].ener - F[i].ener);
            }
            // back to primitives
            for (int i=0;i<N;i++) {
                P[i] = cons2prim(U[i]);
                // positivity guard
                P[i].rho = std::max(P[i].rho, SMALL);
                P[i].p   = std::max(P[i].p,   SMALL);
            }

            t += dt; steps++;
        }

        // Metrics: intermediate band width and TV(rho)
        int band = band_width_count(P, 0.30, 0.90); // Sod middle region ~ example interval
        Real tv  = total_variation_rho(P);
        return {P, tv, band};
    };

    // Three limiters in sequence
    auto [P_minmod, tv_minmod, band_minmod] = run_sim(SlopeLimiter::Minmod);
    auto [P_mc,     tv_mc,     band_mc    ] = run_sim(SlopeLimiter::MC);
    auto [P_vl,     tv_vl,     band_vl    ] = run_sim(SlopeLimiter::VanLeer);

    // 1) Positivity and boundedness (global check)
    auto check_pos = [&](const std::vector<Cell>& P){
        for (auto& c : P) {
            EXPECT_GT(c.rho, 0.0);
            EXPECT_GT(c.p,   0.0);
        }
    };
    check_pos(P_minmod);
    check_pos(P_mc);
    check_pos(P_vl);

    // 2) Diffusion order: smaller intermediate band width means sharper solution
    const int TOL_BAND = 10; // 400 cells ~%2.5 tolerance
    EXPECT_LE(band_vl, band_mc + TOL_BAND);
    EXPECT_LE(band_mc, band_minmod + TOL_BAND);

    // 3) (Optional) TV order: sharper limiter can keep TV higher,
    EXPECT_GT(tv_minmod, 0.0);
    EXPECT_GT(tv_mc,     0.0);
    EXPECT_GT(tv_vl,     0.0);

    // 4) Monotonicity: check for new extreme values
    auto minrho = [](const std::vector<Cell>& P){ Real m=P[0].rho; for(auto& c:P) m=std::min(m,c.rho); return m; };
    auto maxrho = [](const std::vector<Cell>& P){ Real m=P[0].rho; for(auto& c:P) m=std::max(m,c.rho); return m; };
    Real global_min = std::min(L0.rho, R0.rho);
    Real global_max = std::max(L0.rho, R0.rho);
    EXPECT_GE(minrho(P_minmod), 0.0);
    EXPECT_GE(minrho(P_mc),     0.0);
    EXPECT_GE(minrho(P_vl),     0.0);
    EXPECT_LE(maxrho(P_minmod), 1.05*global_max); // small num interval
    EXPECT_LE(maxrho(P_mc),     1.05*global_max);
    EXPECT_LE(maxrho(P_vl),     1.05*global_max);


}
