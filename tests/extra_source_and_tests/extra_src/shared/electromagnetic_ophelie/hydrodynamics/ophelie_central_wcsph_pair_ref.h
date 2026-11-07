/**
 * @file ophelie_central_wcsph_pair_ref.h
 * @brief Doc203 O0 — independent host reference for central WCSPH pair RHS.
 *        Does NOT call AcousticStep / RiemannSolver (forbidden as reference).
 */
#ifndef OPHELIE_CENTRAL_WCSPH_PAIR_REF_H
#define OPHELIE_CENTRAL_WCSPH_PAIR_REF_H

#include <cmath>
#include <cstddef>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{
namespace hydro
{

/** Wendland C2 in 2D: W(r,h), dW/dr (negative for r>0). Independent of framework Kernel. */
inline void wendlandC2_2D(double r, double h, double &W, double &dWdr)
{
    const double q = r / h;
    // Standard SPH Wendland C2 2D: α_d (1-q/2)^4 (2q+1) for q<=2, α_d = 7/(4π h^2)
    const double alpha = 7.0 / (4.0 * M_PI * h * h);
    if (q >= 2.0)
    {
        W = 0.0;
        dWdr = 0.0;
        return;
    }
    const double a = 1.0 - 0.5 * q;
    const double a2 = a * a;
    const double a3 = a2 * a;
    const double a4 = a2 * a2;
    W = alpha * a4 * (2.0 * q + 1.0);
    // dW/dq = alpha * [ -2 a^3 (2q+1) + 2 a^4 ] = alpha * 2 a^3 (-(2q+1) + a)
    //       = alpha * 2 a^3 (-2q -1 + 1 - q/2) = alpha * 2 a^3 (-2.5 q) = -5 alpha a^3 q
    const double dWdq = -5.0 * alpha * a3 * q;
    dWdr = (r > 0.0) ? dWdq / h : 0.0;
}

struct PairState2D
{
    double x = 0, y = 0;
    double u = 0, v = 0;
    double p = 0;
    double rho = 1.0;
    double Vol = 1.0;
    double mass = 1.0;
};

struct CentralRhs2D
{
    double ax = 0, ay = 0;     // acceleration a^p
    double force_x = 0, force_y = 0;
    double drho_dt = 0;
};

/**
 * Paper-form central RHS at particle i (full neighbor sum, including incomplete support).
 * a_i = - Σ_j m_j (p_i+p_j)/(ρ_i ρ_j) ∇_i W
 * ρ̇_i = ρ_i Σ_j V_j (u_i−u_j)·∇_i W
 * ∇_i W = (dW/dr) * e_ij, e_ij = (x_i−x_j)/r  (SPH convention toward j from i? )
 *
 * Framework AcousticStep uses e_ij = e_ij(i,j) and force -= AverageP * 2 * dW_ij V_j * e_ij
 * with AverageP=(p_i+p_j)/2 → Force_i = -V_i Σ (p_i+p_j) V_j (dW_ij e_ij).
 * Here dW_ij is |∇W| scalar from Neighbor: typically positive dW_ij * e_ij construction
 * matching framework: we use grad_i W = dWdr * (r_i - r_j)/r  with dWdr<=0 so grad points.
 */
inline CentralRhs2D centralRhsPaperForm(size_t i, const std::vector<PairState2D> &P, double h)
{
    CentralRhs2D out;
    const auto &pi = P[i];
    for (size_t j = 0; j < P.size(); ++j)
    {
        if (j == i)
            continue;
        const auto &pj = P[j];
        const double dx = pi.x - pj.x;
        const double dy = pi.y - pj.y;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r < 1.0e-30)
            continue;
        double W, dWdr;
        wendlandC2_2D(r, h, W, dWdr);
        (void)W;
        // ∇_i W = dW/dr * (r_i - r_j)/r
        const double gx = dWdr * dx / r;
        const double gy = dWdr * dy / r;
        const double psum = pi.p + pj.p;
        const double coef = pj.mass * psum / (pi.rho * pj.rho);
        out.ax -= coef * gx;
        out.ay -= coef * gy;
        const double du = pi.u - pj.u;
        const double dv = pi.v - pj.v;
        out.drho_dt += pi.rho * pj.Vol * (du * gx + dv * gy);
    }
    out.force_x = out.ax * pi.mass;
    out.force_y = out.ay * pi.mass;
    return out;
}

/**
 * Emulate NoRiemann AcousticStep1st force + 2nd continuity (algebra mirror, not calling solvers):
 * Force_i = - V_i Σ_j (p_i+p_j) V_j * (dW_ij * e_ij)
 * where dW_ij in framework is positive scalar |dW| used with e_ij; equivalent to
 * using ∇W = dWdr*(ri-rj)/r with dWdr<=0 and Force = -V_i Σ (p_i+p_j) V_j ∇W
 * when ∇W has the same sign as framework dW_ij*e_ij.
 *
 * Framework: dW_ij(i,j) returns Kernel::dW_ij which is typically positive * direction encoded in e_ij.
 * We mirror paper-form with mass=ρV → identical when ρ_i=ρ_j constant path.
 */
inline CentralRhs2D centralRhsNoRiemannMirror(size_t i, const std::vector<PairState2D> &P, double h)
{
    CentralRhs2D out;
    const auto &pi = P[i];
    for (size_t j = 0; j < P.size(); ++j)
    {
        if (j == i)
            continue;
        const auto &pj = P[j];
        const double dx = pi.x - pj.x;
        const double dy = pi.y - pj.y;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r < 1.0e-30)
            continue;
        double W, dWdr;
        wendlandC2_2D(r, h, W, dWdr);
        (void)W;
        const double gx = dWdr * dx / r;
        const double gy = dWdr * dy / r;
        // Force_i += - V_i (p_i+p_j) V_j ∇W
        out.force_x -= pi.Vol * (pi.p + pj.p) * pj.Vol * gx;
        out.force_y -= pi.Vol * (pi.p + pj.p) * pj.Vol * gy;
        const double du = pi.u - pj.u;
        const double dv = pi.v - pj.v;
        out.drho_dt += pi.rho * pj.Vol * (du * gx + dv * gy);
    }
    out.ax = out.force_x / pi.mass;
    out.ay = out.force_y / pi.mass;
    return out;
}

inline double relErr(double a, double b, double floor_abs = 1.0e-30)
{
    const double den = std::max(std::max(std::abs(a), std::abs(b)), floor_abs);
    return std::abs(a - b) / den;
}

} // namespace hydro
} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
