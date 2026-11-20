/**
 * @file electromagnetic_ophelie_phi_d2_local_poly.h
 * @brief D2 local quadratic Neumann collocation for φ (historically Doc217 E-G): host-side D2 local quadratic Neumann collocation for φ.
 *
 * Discrete family matches docs/ophelie/reference/doc217_runs/scripts/d2_cpu_local_poly_deepfix.py:
 *   interior: σΔφ + ∇σ·∇φ = −div(σF)
 *   boundary: n·∇φ = −n·F at true wall point (side r=R / caps z=±H)
 *   Doc218: F must be evaluated at the same wall point (wall_reconstructed);
 *           center_legacy / boundary_mean_legacy kept only as diagnostic modes.
 *   row scaling, adaptive neighborhood, pin+demean gauge, SparseLU.
 *
 * Not a full device Krylov solve — accuracy-first host SparseLU inside the GPU binary
 * (same honesty class as Doc216 P1 host Eigen). Reconstruct E/J/Q with the SAME local
 * polynomial gradient (plan §5.3.8), not SPH G.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_PHI_D2_LOCAL_POLY_H
#define ELECTROMAGNETIC_OPHELIE_PHI_D2_LOCAL_POLY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

namespace doc217_d2_detail
{
constexpr int kNBase = 10;
constexpr int kMinNbr = kNBase + 2;

inline void fillBasis(const Real *xi /*3*/, Real *B /*10*/)
{
    const Real x = xi[0], y = xi[1], z = xi[2];
    B[0] = Real(1);
    B[1] = x;
    B[2] = y;
    B[3] = z;
    B[4] = Real(0.5) * x * x;
    B[5] = Real(0.5) * y * y;
    B[6] = Real(0.5) * z * z;
    B[7] = x * y;
    B[8] = x * z;
    B[9] = y * z;
}

struct CellList
{
    Vecd mn{};
    Real inv = Real(1);
    std::vector<int> ijk; // 3*n
    std::unordered_map<long long, std::vector<int>> buckets;
    static long long key(int ix, int iy, int iz)
    {
        return ((long long)(ix + 100000) << 42) ^ ((long long)(iy + 100000) << 21) ^ (long long)(iz + 100000);
    }
};

inline CellList buildCellList(const Vecd *pos, size_t n, Real cutoff)
{
    CellList cl;
    cl.mn = pos[0];
    Vecd mx = pos[0];
    for (size_t i = 1; i < n; ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            cl.mn[k] = std::min(cl.mn[k], pos[i][k]);
            mx[k] = std::max(mx[k], pos[i][k]);
        }
    }
    cl.mn -= Vecd(Real(1e-12), Real(1e-12), Real(1e-12));
    cl.inv = Real(1) / cutoff;
    cl.ijk.assign(3 * n, 0);
    for (size_t i = 0; i < n; ++i)
    {
        const int ix = int(std::floor((pos[i][0] - cl.mn[0]) * cl.inv));
        const int iy = int(std::floor((pos[i][1] - cl.mn[1]) * cl.inv));
        const int iz = int(std::floor((pos[i][2] - cl.mn[2]) * cl.inv));
        cl.ijk[3 * i] = ix;
        cl.ijk[3 * i + 1] = iy;
        cl.ijk[3 * i + 2] = iz;
        cl.buckets[CellList::key(ix, iy, iz)].push_back(int(i));
    }
    return cl;
}

inline void neighborsAt(int i, const Vecd *pos, Real cutoff, const CellList &cl, std::vector<int> &out)
{
    out.clear();
    const int ix = cl.ijk[3 * i], iy = cl.ijk[3 * i + 1], iz = cl.ijk[3 * i + 2];
    const Real r2max = cutoff * cutoff + Real(1e-30);
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
            {
                auto it = cl.buckets.find(CellList::key(ix + dx, iy + dy, iz + dz));
                if (it == cl.buckets.end())
                    continue;
                for (int j : it->second)
                {
                    const Vecd d = pos[j] - pos[i];
                    if (d.squaredNorm() <= r2max)
                        out.push_back(j);
                }
            }
}

struct FitState
{
    Eigen::Matrix<Real, kNBase, kNBase> Minv{};
    Eigen::Matrix<Real, kNBase, Eigen::Dynamic> BtW{};
    std::vector<int> ids;
    Real h = 0;
    Real cond = 0;
};

inline bool fitMoment(const std::vector<int> &ids, const Vecd *pos, const Vecd &center, Real h, FitState &out)
{
    if (int(ids.size()) < kMinNbr)
        return false;
    const int m = int(ids.size());
    Eigen::Matrix<Real, Eigen::Dynamic, kNBase> B(m, kNBase);
    Eigen::Matrix<Real, Eigen::Dynamic, 1> w(m);
    for (int a = 0; a < m; ++a)
    {
        const Vecd xi = (pos[ids[a]] - center) / h;
        Real row[kNBase];
        const Real xiv[3] = {xi[0], xi[1], xi[2]};
        fillBasis(xiv, row);
        for (int k = 0; k < kNBase; ++k)
            B(a, k) = row[k];
        w(a) = std::exp(Real(-2) * xi.squaredNorm());
    }
    Eigen::Matrix<Real, kNBase, kNBase> Mt = Eigen::Matrix<Real, kNBase, kNBase>::Zero();
    for (int a = 0; a < m; ++a)
        Mt += w(a) * B.row(a).transpose() * B.row(a);
    const Real ridge = Real(1e-12) * (Mt.trace() / Real(kNBase) + Real(1));
    Mt += ridge * Eigen::Matrix<Real, kNBase, kNBase>::Identity();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, kNBase, kNBase>> es(Mt);
    if (es.info() != Eigen::Success)
        return false;
    const auto &evals = es.eigenvalues();
    const Real smax = evals.maxCoeff();
    const Real smin = evals.minCoeff();
    out.cond = Real(smax / (smin + Real(1e-30)));
    if (out.cond > Real(1e8) || smin < Real(1e-14) * smax)
        return false;
    Eigen::FullPivLU<Eigen::Matrix<Real, kNBase, kNBase>> lu(Mt);
    if (!lu.isInvertible())
        return false;
    out.Minv = lu.inverse();
    out.BtW.resize(kNBase, m);
    for (int a = 0; a < m; ++a)
        out.BtW.col(a) = w(a) * B.row(a).transpose();
    out.ids = ids;
    out.h = h;
    return true;
}

} // namespace doc217_d2_detail

struct OphelieDoc217D2Report
{
    Real linear_rel_residual = std::numeric_limits<Real>::quiet_NaN();
    Real linear_abs_residual = std::numeric_limits<Real>::quiet_NaN();
    Real b_norm = 0;
    Real true_rel_residual = std::numeric_limits<Real>::quiet_NaN(); // ||r||/||b|| (no floor)
    int rank_fail = 0;
    Real cond_med = 0;
    Real nbr_min = 0;
    Real compat_correction = 0;
    Real remaining_field_ratio = std::numeric_limits<Real>::quiet_NaN();
    bool solve_ok = false; // Doc219: SparseLU Success + finite φ; false → do not consume Q
    int sparselu_info = -1;
    const char *fail_stage = ""; // "", "grad_sigma_fit", "assembly_fit", "sparselu", "inject"
    /** Doc219 §5.2: pin-system vs unpinned physical residual (demeaned φ). */
    Real pin_rel_residual = std::numeric_limits<Real>::quiet_NaN();
    Real phys_rel_residual = std::numeric_limits<Real>::quiet_NaN();
    Real phys_abs_residual = std::numeric_limits<Real>::quiet_NaN();
    Real phys_b_norm = std::numeric_limits<Real>::quiet_NaN();
    Real phys_interior_abs = std::numeric_limits<Real>::quiet_NaN();
    Real phys_boundary_abs = std::numeric_limits<Real>::quiet_NaN();
    Real phys_ones_defect = std::numeric_limits<Real>::quiet_NaN(); // ||K_phys·1|| / (||K||scale)
};

struct OphelieDoc217D2ComplexBundle
{
    OphelieDoc217D2Report imag;
    OphelieDoc217D2Report realp;
    bool solved_real = false;
    bool ok = false; // imag (and real if solved) all solve_ok && rank_fail==0
};

/** Doc219: per-particle adaptive fits shared by assembly and recon (no silent h0 fallback). */
struct OphelieDoc217D2FitCache
{
    std::vector<doc217_d2_detail::FitState> fits;
    std::vector<char> ok;
    size_t n = 0;
    bool valid = false;
};

inline OphelieDoc217D2FitCache &ophelieDoc217D2FitCache()
{
    static OphelieDoc217D2FitCache cache;
    return cache;
}

/**
 * Classify near-wall particles; write wall points + outward normals.
 * side takes priority over caps (plan §5.3.6).
 */
inline void hostOphelieDoc217D2Classify(const Vecd *pos, size_t n, const Vecd &center, Real R, Real H, Real dp,
                                        Real wall_dp, std::vector<char> &surf, std::vector<Vecd> &nvec,
                                        std::vector<Vecd> &wall)
{
    surf.assign(n, 0);
    nvec.assign(n, Vecd::Zero());
    wall.assign(n, Vecd::Zero());
    for (size_t i = 0; i < n; ++i)
    {
        wall[i] = pos[i];
        const Vecd x = pos[i] - center;
        const Real r = std::hypot(x[0], x[1]);
        const Real z = x[2];
        const bool is_side = std::abs(r - R) <= wall_dp * dp;
        const bool is_top = (!is_side) && (std::abs(z - H) <= wall_dp * dp);
        const bool is_bot = (!is_side) && (std::abs(z + H) <= wall_dp * dp);
        if (!(is_side || is_top || is_bot))
            continue;
        surf[i] = 1;
        if (is_side)
        {
            const Real rh = std::max(r, Real(1e-30));
            wall[i][0] = center[0] + R * x[0] / rh;
            wall[i][1] = center[1] + R * x[1] / rh;
            nvec[i][0] = x[0] / rh;
            nvec[i][1] = x[1] / rh;
        }
        else if (is_top)
        {
            wall[i][2] = center[2] + H;
            nvec[i][2] = Real(1);
        }
        else
        {
            wall[i][2] = center[2] - H;
            nvec[i][2] = Real(-1);
        }
    }
}

using Doc217D2Fit = doc217_d2_detail::FitState;

inline Real hostOphelieDoc217D2EvalValue(const Doc217D2Fit &fit, const Real *vals /*on fit.ids*/,
                                         const Vecd &x_eval, const Vecd &pos_i)
{
    using namespace doc217_d2_detail;
    Eigen::Matrix<Real, kNBase, 1> rhs = Eigen::Matrix<Real, kNBase, 1>::Zero();
    for (int a = 0; a < int(fit.ids.size()); ++a)
        rhs += fit.BtW.col(a) * vals[a];
    const Eigen::Matrix<Real, kNBase, 1> c = fit.Minv * rhs;
    const Vecd xi = (x_eval - pos_i) / fit.h;
    Real B[kNBase];
    const Real xv[3] = {xi[0], xi[1], xi[2]};
    fillBasis(xv, B);
    Real v = Real(0);
    for (int k = 0; k < kNBase; ++k)
        v += B[k] * c(k);
    return v;
}

inline Vecd hostOphelieDoc217D2EvalGrad(const Doc217D2Fit &fit, const Real *phi_vals /*on fit.ids*/,
                                        const Vecd &x_eval, const Vecd &pos_i)
{
    using namespace doc217_d2_detail;
    Eigen::Matrix<Real, kNBase, 1> rhs = Eigen::Matrix<Real, kNBase, 1>::Zero();
    for (int a = 0; a < int(fit.ids.size()); ++a)
        rhs += fit.BtW.col(a) * phi_vals[a];
    const Eigen::Matrix<Real, kNBase, 1> c = fit.Minv * rhs;
    const Vecd xi = (x_eval - pos_i) / fit.h;
    const Real x = xi[0], y = xi[1], z = xi[2];
    // dB/dξ
    Real dBx[kNBase] = {0, 1, 0, 0, x, 0, 0, y, z, 0};
    Real dBy[kNBase] = {0, 0, 1, 0, 0, y, 0, x, 0, z};
    Real dBz[kNBase] = {0, 0, 0, 1, 0, 0, z, 0, x, y};
    Vecd g = Vecd::Zero();
    for (int k = 0; k < kNBase; ++k)
    {
        g[0] += dBx[k] * c(k);
        g[1] += dBy[k] * c(k);
        g[2] += dBz[k] * c(k);
    }
    return g / fit.h;
}

/** Doc219: same adaptive neighborhood for ∇σ / assembly / recon. */
inline bool hostOphelieDoc217D2AdaptiveFit(
    size_t i, const Vecd *pos, Real h0, Real h_max,
    const std::unordered_map<long long, doc217_d2_detail::CellList> &cells, doc217_d2_detail::FitState &fit,
    Real &h_used)
{
    using namespace doc217_d2_detail;
    std::vector<int> ids;
    h_used = h0;
    for (Real h = h0; h <= h_max * Real(1.0001);)
    {
        const long long key = (long long)std::llround(double(h) * 1e8);
        auto it = cells.find(key);
        if (it == cells.end())
            break;
        neighborsAt(int(i), pos, h, it->second, ids);
        if (fitMoment(ids, pos, pos[i], h, fit))
        {
            h_used = h;
            return true;
        }
        if (h >= h_max)
            break;
        h = std::min(h * Real(1.25), h_max);
    }
    return false;
}

inline OphelieDoc217D2Report hostOphelieDoc217D2SolveOneComponent(
    const Vecd *pos, const Real *Vol, const Real *sigma, const Vecd *A, Real *phi_out, size_t n, Real a_sign,
    Real omega, const OphelieParameters &params)
{
    using namespace doc217_d2_detail;
    OphelieDoc217D2Report report;
    const Vecd center = params.doc217_d2_center_;
    const Real R = params.doc217_d2_radius_;
    const Real H = params.doc217_d2_half_height_;
    const Real dp = params.doc217_d2_dp_;
    const Real h0 = params.doc217_d2_h_factor_ * dp;
    const Real h_max = params.doc217_d2_h_max_factor_ * dp;

    std::vector<char> surf;
    std::vector<Vecd> nvec, wall;
    hostOphelieDoc217D2Classify(pos, n, center, R, H, dp, Real(0.75), surf, nvec, wall);

    std::vector<Vecd> F(n), grad_sigma(n, Vecd::Zero());
    Real sigma_char = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        F[i] = a_sign * omega * A[i];
        sigma_char += sigma[i];
    }
    sigma_char /= Real(std::max<size_t>(n, 1));

    // Shared cell lists for adaptive h
    std::unordered_map<long long, CellList> cells;
    for (Real h = h0; h <= h_max * Real(1.0001);)
    {
        cells.emplace((long long)std::llround(double(h) * 1e8), buildCellList(pos, n, h));
        if (h >= h_max)
            break;
        h = std::min(h * Real(1.25), h_max);
    }
    cells.emplace((long long)std::llround(double(h_max) * 1e8), buildCellList(pos, n, h_max));

    // Doc219: ∇σ uses the SAME adaptive fit; failure is hard (no silent ∇σ=0).
    {
        FitState fit;
        Real h_used = h0;
        for (size_t i = 0; i < n; ++i)
        {
            if (params.doc219_inject_fit_fail_ && i == n / 2)
            {
                report.rank_fail = 1;
                report.fail_stage = "inject";
                report.solve_ok = false;
                std::cout << "[ophelie][doc219][d2] FAIL inject fit at i=" << i << std::endl;
                for (size_t j = 0; j < n; ++j)
                    phi_out[j] = Real(0);
                return report;
            }
            if (!hostOphelieDoc217D2AdaptiveFit(i, pos, h0, h_max, cells, fit, h_used))
            {
                report.rank_fail = 1;
                report.fail_stage = "grad_sigma_fit";
                report.solve_ok = false;
                std::cout << "[ophelie][doc219][d2] FAIL grad_sigma adaptiveFit at i=" << i
                          << " pos=(" << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << ")"
                          << " — abort (no silent ∇σ=0)" << std::endl;
                for (size_t j = 0; j < n; ++j)
                    phi_out[j] = Real(0);
                return report;
            }
            Eigen::Matrix<Real, kNBase, 1> rhs = Eigen::Matrix<Real, kNBase, 1>::Zero();
            for (int a = 0; a < int(fit.ids.size()); ++a)
                rhs += fit.BtW.col(a) * sigma[fit.ids[a]];
            const Eigen::Matrix<Real, kNBase, 1> c = fit.Minv * rhs;
            grad_sigma[i] = Vecd(c(1), c(2), c(3)) / fit.h;
        }
    }

    // Doc218 §5.3: host local matrix + SparseLU in double (GPU particle fields stay Real).
    using HostScalar = double;
    using HostVec = Eigen::Matrix<HostScalar, Eigen::Dynamic, 1>;
    std::vector<Eigen::Triplet<HostScalar>> trips;
    trips.reserve(n * 40);
    HostVec b = HostVec::Zero(static_cast<Eigen::Index>(n));
    std::vector<Real> conds;
    std::vector<int> nbrs;
    int rank_fail = 0;
    FitState fit;
    auto &fit_cache = ophelieDoc217D2FitCache();
    fit_cache.fits.assign(n, FitState{});
    fit_cache.ok.assign(n, 0);
    fit_cache.n = n;
    fit_cache.valid = false;

    for (size_t i = 0; i < n; ++i)
    {
        Real h_used = h0;
        if (!hostOphelieDoc217D2AdaptiveFit(i, pos, h0, h_max, cells, fit, h_used))
        {
            // Doc218/219: never silently replace a physics row.
            ++rank_fail;
            report.rank_fail = rank_fail;
            report.fail_stage = "assembly_fit";
            report.linear_rel_residual = Real(1);
            report.remaining_field_ratio = Real(1);
            report.solve_ok = false;
            std::cout << "[ophelie][doc219][d2] FAIL assembly adaptiveFit at i=" << i
                      << " pos=(" << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << ")"
                      << " nbr attempt logged — abort (no φ=0 soft pin)" << std::endl;
            for (size_t j = 0; j < n; ++j)
                phi_out[j] = Real(0);
            return report;
        }
        fit_cache.fits[i] = fit;
        fit_cache.ok[i] = 1;
        conds.push_back(fit.cond);
        nbrs.push_back(int(fit.ids.size()));
        const int m = int(fit.ids.size());
        Eigen::Matrix<Real, 1, Eigen::Dynamic> coeffs(m);

        if (surf[i])
        {
            const Real scale = h_used;
            const Vecd xi_w = (wall[i] - pos[i]) / h_used;
            const Real x = xi_w[0], y = xi_w[1], z = xi_w[2];
            Real dBx[kNBase] = {0, 1, 0, 0, x, 0, 0, y, z, 0};
            Real dBy[kNBase] = {0, 0, 1, 0, 0, y, 0, x, 0, z};
            Real dBz[kNBase] = {0, 0, 0, 1, 0, 0, z, 0, x, y};
            Eigen::Matrix<Real, kNBase, 1> row_c = Eigen::Matrix<Real, kNBase, 1>::Zero();
            for (int k = 0; k < kNBase; ++k)
                row_c(k) = (nvec[i][0] * dBx[k] + nvec[i][1] * dBy[k] + nvec[i][2] * dBz[k]) / h_used;
            Eigen::Matrix<Real, 1, kNBase> rc = scale * row_c.transpose();
            coeffs = rc * fit.Minv * fit.BtW;
            // Doc218: Neumann RHS must use F at the same wall point as n·∇φ.
            Vecd F_bc = F[i]; // center_legacy
            if (params.doc218_wall_F_mode_ == 1)
            {
                std::vector<Real> component(static_cast<size_t>(m));
                for (int k = 0; k < 3; ++k)
                {
                    for (int a = 0; a < m; ++a)
                        component[static_cast<size_t>(a)] = F[fit.ids[static_cast<size_t>(a)]][k];
                    F_bc[k] = hostOphelieDoc217D2EvalValue(fit, component.data(), wall[i], pos[i]);
                }
            }
            b(Eigen::Index(i)) = HostScalar(-nvec[i].dot(F_bc) * scale);
            for (int a = 0; a < m; ++a)
                trips.emplace_back(int(i), fit.ids[a], HostScalar(coeffs(a)));
        }
        else
        {
            const Real scale = (h_used * h_used) / (sigma_char + TinyReal);
            Eigen::Matrix<Real, kNBase, 1> row_c = Eigen::Matrix<Real, kNBase, 1>::Zero();
            row_c(4) = sigma[i] / (h_used * h_used);
            row_c(5) = sigma[i] / (h_used * h_used);
            row_c(6) = sigma[i] / (h_used * h_used);
            row_c(1) = grad_sigma[i][0] / h_used;
            row_c(2) = grad_sigma[i][1] / h_used;
            row_c(3) = grad_sigma[i][2] / h_used;
            Eigen::Matrix<Real, 1, kNBase> rc = scale * row_c.transpose();
            coeffs = rc * fit.Minv * fit.BtW;
            Real divF = Real(0);
            for (int k = 0; k < 3; ++k)
            {
                Eigen::Matrix<Real, kNBase, 1> rhs = Eigen::Matrix<Real, kNBase, 1>::Zero();
                for (int a = 0; a < m; ++a)
                    rhs += fit.BtW.col(a) * F[fit.ids[a]][k];
                const Eigen::Matrix<Real, kNBase, 1> ck = fit.Minv * rhs;
                divF += ck(1 + k) / h_used;
            }
            const Real div_sF = sigma[i] * divF + grad_sigma[i].dot(F[i]);
            b(Eigen::Index(i)) = HostScalar(-div_sF * scale);
            for (int a = 0; a < m; ++a)
                trips.emplace_back(int(i), fit.ids[a], HostScalar(coeffs(a)));
        }
    }

    HostScalar b_surf_sum = HostScalar(0);
    int n_surf = 0;
    for (size_t i = 0; i < n; ++i)
        if (surf[i])
        {
            b_surf_sum += b(Eigen::Index(i));
            ++n_surf;
        }
    const HostScalar corr = (n_surf > 0) ? b_surf_sum / HostScalar(n_surf) : HostScalar(0);
    report.compat_correction = Real(corr);
    // Doc218: boundary_mean_legacy (mode 0) only; formal path mode 1 = none.
    if (params.doc218_compat_mode_ == 0 && n_surf > 0)
        for (size_t i = 0; i < n; ++i)
            if (surf[i])
                b(Eigen::Index(i)) -= corr;

    // Doc219 §5.2: keep unpinned physical system before gauge pin.
    const std::vector<Eigen::Triplet<HostScalar>> trips_phys = trips;
    const HostVec b_phys = b;

    // Gauge pin: drop row-0 entries from triplets, then add identity.
    {
        std::vector<Eigen::Triplet<HostScalar>> pinned;
        pinned.reserve(trips.size() + 1);
        for (const auto &t : trips)
            if (t.row() != 0)
                pinned.push_back(t);
        pinned.emplace_back(0, 0, HostScalar(1));
        trips.swap(pinned);
        b(0) = HostScalar(0);
    }
    const Eigen::Index n_idx = static_cast<Eigen::Index>(n);
    Eigen::SparseMatrix<HostScalar> K(n_idx, n_idx);
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();
    Eigen::SparseMatrix<HostScalar> K_phys(n_idx, n_idx);
    K_phys.setFromTriplets(trips_phys.begin(), trips_phys.end());
    K_phys.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<HostScalar>> solver;
    solver.compute(K);
    report.sparselu_info = int(solver.info());
    HostVec phi = HostVec::Zero(n_idx);
    if (params.doc219_inject_solve_fail_)
    {
        report.solve_ok = false;
        report.fail_stage = "inject";
        report.sparselu_info = -99;
        std::cout << "[ophelie][doc219][d2] FAIL inject SparseLU failure" << std::endl;
        for (size_t j = 0; j < n; ++j)
            phi_out[j] = Real(0);
        fit_cache.valid = false;
        return report;
    }
    if (solver.info() != Eigen::Success)
    {
        report.solve_ok = false;
        report.fail_stage = "sparselu";
        std::cout << "[ophelie][doc219][d2] FAIL SparseLU compute info=" << report.sparselu_info
                  << " — abort (no WARN-and-continue)" << std::endl;
        for (size_t j = 0; j < n; ++j)
            phi_out[j] = Real(0);
        fit_cache.valid = false;
        return report;
    }
    phi = solver.solve(b);
    if (solver.info() != Eigen::Success || !phi.allFinite())
    {
        report.solve_ok = false;
        report.fail_stage = "sparselu";
        report.sparselu_info = int(solver.info());
        std::cout << "[ophelie][doc219][d2] FAIL SparseLU solve info=" << report.sparselu_info
                  << " finite=" << phi.allFinite() << std::endl;
        for (size_t j = 0; j < n; ++j)
            phi_out[j] = Real(0);
        fit_cache.valid = false;
        return report;
    }

    // Algebraic residual BEFORE demean (demean breaks the pin row and inflates
    // rel residual when ||b||~0, e.g. pure rotational A with φ≈0).
    const HostVec r_alg = K * phi - b;
    report.linear_abs_residual = Real(r_alg.norm());
    report.b_norm = Real(b.norm());
    report.true_rel_residual =
        report.b_norm > TinyReal ? report.linear_abs_residual / report.b_norm : report.linear_abs_residual;
    report.pin_rel_residual = report.true_rel_residual;
    // Doc217 near-zero RHS: ||r|| ≤ atol + rtol||b|| ; report effective rel for gates.
    const Real atol = Real(1.0e-8) * (Real(1) + report.b_norm);
    const Real rtol = Real(1.0e-6);
    const Real denom = std::max(report.b_norm, atol / std::max(rtol, Real(1.0e-30)));
    report.linear_rel_residual = report.linear_abs_residual / (denom + TinyReal);
    report.solve_ok = true;
    fit_cache.valid = true;

    Real vol_sum = Real(0), phiv = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        vol_sum += Vol[i];
        phiv += Vol[i] * Real(phi(Eigen::Index(i)));
    }
    const Real mean = phiv / (vol_sum + TinyReal);
    HostVec phi_demean = phi;
    for (size_t i = 0; i < n; ++i)
    {
        phi_out[i] = Real(phi(Eigen::Index(i))) - mean;
        phi_demean(Eigen::Index(i)) = HostScalar(phi_out[i]);
    }

    // Physical residual of demeaned φ on unpinned rows (includes former pin row).
    {
        const HostVec r_phys = K_phys * phi_demean - b_phys;
        report.phys_abs_residual = Real(r_phys.norm());
        report.phys_b_norm = Real(b_phys.norm());
        report.phys_rel_residual = report.phys_b_norm > TinyReal
                                       ? report.phys_abs_residual / report.phys_b_norm
                                       : report.phys_abs_residual;
        HostScalar int2 = HostScalar(0), bnd2 = HostScalar(0);
        for (size_t i = 0; i < n; ++i)
        {
            const HostScalar ri = r_phys(Eigen::Index(i));
            if (surf[i])
                bnd2 += ri * ri;
            else
                int2 += ri * ri;
        }
        report.phys_interior_abs = Real(std::sqrt(double(int2)));
        report.phys_boundary_abs = Real(std::sqrt(double(bnd2)));
        const HostVec ones = HostVec::Ones(n_idx);
        const HostVec k1 = K_phys * ones;
        report.phys_ones_defect = Real(k1.norm()) / (Real(1) + report.phys_b_norm);
    }
    std::cout << "[ophelie][doc219][d2] host_SparseLU_precision=double"
              << " fit_precision=Real assembly_precision=double_from_Real solve_precision=double"
              << " output_precision=Real"
              << " pin_true_rel=" << report.pin_rel_residual
              << " effective_rel=" << report.linear_rel_residual
              << " phys_rel=" << report.phys_rel_residual << " phys_abs=" << report.phys_abs_residual
              << " phys_b=" << report.phys_b_norm << " phys_int=" << report.phys_interior_abs
              << " phys_bnd=" << report.phys_boundary_abs << " ones_defect=" << report.phys_ones_defect
              << std::endl;

    report.rank_fail = rank_fail;
    if (!conds.empty())
    {
        std::vector<Real> tmp = conds;
        std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
        report.cond_med = tmp[tmp.size() / 2];
    }
    if (!nbrs.empty())
        report.nbr_min = Real(*std::min_element(nbrs.begin(), nbrs.end()));

    // remaining_field_ratio: reuse assembly fit cache (Doc219 — no silent h0 fallback).
    Real f2 = Real(0), r2 = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        if (!fit_cache.ok[i])
            continue;
        const FitState &fi = fit_cache.fits[i];
        std::vector<Real> phiv_local(fi.ids.size());
        for (size_t a = 0; a < fi.ids.size(); ++a)
            phiv_local[a] = phi_out[fi.ids[a]];
        const Vecd g = hostOphelieDoc217D2EvalGrad(fi, phiv_local.data(), pos[i], pos[i]);
        const Vecd rem = g + F[i];
        f2 += sigma[i] * F[i].squaredNorm();
        r2 += sigma[i] * rem.squaredNorm();
    }
    report.remaining_field_ratio = std::sqrt(r2 / (f2 + TinyReal));
    return report;
}

template <class ExecutionPolicy>
inline OphelieDoc217D2ComplexBundle execOphelieDoc217D2LocalPolySolveComplex(
    SolidBody &glass_body, Inner<> & /*inner*/, const OphelieGlassFieldNames &names, OphelieParameters &params,
    bool skip_real = false)
{
    OphelieDoc217D2ComplexBundle out;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    std::cout << "[ophelie][doc217][d2] local quadratic Neumann collocation ON (host SparseLU)"
              << " h0=" << params.doc217_d2_h_factor_ << "dp"
              << " h_max=" << params.doc217_d2_h_max_factor_ << "dp"
              << " R=" << params.doc217_d2_radius_ << " H=" << params.doc217_d2_half_height_
              << " dp=" << params.doc217_d2_dp_ << " n=" << n
              << " wall_F_mode=" << (params.doc218_wall_F_mode_ == 1 ? "wall_reconstructed" : "center_legacy")
              << " compat_mode=" << (params.doc218_compat_mode_ == 1 ? "none" : "boundary_mean_legacy")
              << std::endl;

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Vecd *Ar = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);

    out.imag = hostOphelieDoc217D2SolveOneComponent(pos, Vol, sigma, Ar, phi_i, n, Real(1), params.omega(), params);
    syncVariableToDevice<Real>(particles, names.phi_imag);
    std::cout << "[ophelie][doc217][d2] imag linear_rel=" << out.imag.linear_rel_residual
              << " ||r||=" << out.imag.linear_abs_residual << " ||b||=" << out.imag.b_norm
              << " rank_fail=" << out.imag.rank_fail << " nbr_min=" << out.imag.nbr_min
              << " compat_corr=" << out.imag.compat_correction
              << " remaining_field_ratio=" << out.imag.remaining_field_ratio
              << " solve_ok=" << (out.imag.solve_ok ? 1 : 0)
              << " fail_stage=" << (out.imag.fail_stage[0] ? out.imag.fail_stage : "none") << std::endl;

    if (!out.imag.solve_ok || out.imag.rank_fail != 0)
    {
        out.ok = false;
        return out;
    }

    if (!skip_real && params.edge_flux_complex_)
    {
        syncVariableToHost<Vecd>(particles, names.a_src_imag);
        Vecd *Ai = particles.getVariableDataByName<Vecd>(names.a_src_imag);
        Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
        out.realp =
            hostOphelieDoc217D2SolveOneComponent(pos, Vol, sigma, Ai, phi_r, n, Real(-1), params.omega(), params);
        out.solved_real = true;
        syncVariableToDevice<Real>(particles, names.phi_real);
        std::cout << "[ophelie][doc217][d2] real linear_rel=" << out.realp.linear_rel_residual
                  << " ||r||=" << out.realp.linear_abs_residual << " ||b||=" << out.realp.b_norm
                  << " rank_fail=" << out.realp.rank_fail
                  << " remaining_field_ratio=" << out.realp.remaining_field_ratio
                  << " solve_ok=" << (out.realp.solve_ok ? 1 : 0) << std::endl;
        if (!out.realp.solve_ok || out.realp.rank_fail != 0)
        {
            out.ok = false;
            return out;
        }
    }
    out.ok = true;
    return out;
}

/** Same-family recon: E=-(∇_poly φ + a_sign ω A). */
template <class ExecutionPolicy>
inline void execOphelieDoc217D2LocalPolyRecon(SolidBody &glass_body, Inner<> & /*inner*/,
                                              const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                              bool skip_real = false)
{
    using namespace doc217_d2_detail;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real h0 = params.doc217_d2_h_factor_ * params.doc217_d2_dp_;
    const Real h_max = params.doc217_d2_h_max_factor_ * params.doc217_d2_dp_;
    const Real omega = params.omega();

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
    Vecd *Ar = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Vecd *e_i = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    Vecd *j_i = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    Real *q_i = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);

    auto cl0 = buildCellList(pos, n, h0);
    auto clm = buildCellList(pos, n, h_max);
    std::vector<int> ids;
    FitState fit;
    std::cout << "[ophelie][doc217][d2] recon=local_poly E=-(grad_poly phi + a_sign omega A)" << std::endl;

    auto &fit_cache = ophelieDoc217D2FitCache();
    const bool use_cache = fit_cache.valid && fit_cache.n == n;
    if (!use_cache)
        std::cout << "[ophelie][doc219][d2] WARN recon without assembly fit cache — adaptive refit"
                  << std::endl;

    auto recon_one = [&](const Real *phi, const Vecd *A, Real a_sign, Vecd *e, Vecd *j, Real *q) {
        for (size_t i = 0; i < n; ++i)
        {
            Vecd g = Vecd::Zero();
            bool have = false;
            if (use_cache && fit_cache.ok[i])
            {
                fit = fit_cache.fits[i];
                have = true;
            }
            else
            {
                neighborsAt(int(i), pos, h0, cl0, ids);
                Real h = h0;
                if (int(ids.size()) < kMinNbr)
                {
                    neighborsAt(int(i), pos, h_max, clm, ids);
                    h = h_max;
                }
                have = fitMoment(ids, pos, pos[i], h, fit);
            }
            if (!have)
            {
                // Doc219: no silent gradφ=0 → zero E/J for a failed particle.
                std::cout << "[ophelie][doc219][d2] FAIL recon fit at i=" << i
                          << " — zeroing E/J/Q and marking cache invalid" << std::endl;
                fit_cache.valid = false;
                e[i] = Vecd::Zero();
                j[i] = Vecd::Zero();
                q[i] = Real(0);
                // Poison remaining so thermal cannot treat as valid heat source.
                for (size_t k = i; k < n; ++k)
                {
                    e[k] = Vecd::Zero();
                    j[k] = Vecd::Zero();
                    q[k] = std::numeric_limits<Real>::quiet_NaN();
                }
                return;
            }
            std::vector<Real> pv(fit.ids.size());
            for (size_t a = 0; a < fit.ids.size(); ++a)
                pv[a] = phi[fit.ids[a]];
            g = hostOphelieDoc217D2EvalGrad(fit, pv.data(), pos[i], pos[i]);
            const Vecd e_loc = -(g + a_sign * omega * A[i]);
            e[i] = e_loc;
            j[i] = sigma[i] * e_loc;
            q[i] = Real(0.5) * sigma[i] * e_loc.squaredNorm();
        }
    };

    recon_one(phi_i, Ar, Real(1), e_i, j_i, q_i);
    syncVariableToDevice<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToDevice<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToDevice<Real>(particles, names.joule_heat_edge_recon_imag);

    if (!skip_real && params.edge_flux_complex_)
    {
        syncVariableToHost<Real>(particles, names.phi_real);
        syncVariableToHost<Vecd>(particles, names.a_src_imag);
        Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
        Vecd *Ai = particles.getVariableDataByName<Vecd>(names.a_src_imag);
        Vecd *e_r = particles.getVariableDataByName<Vecd>(names.e_edge_recon_real);
        Vecd *j_r = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
        Real *q_r = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_real);
        recon_one(phi_r, Ai, Real(-1), e_r, j_r, q_r);
        syncVariableToDevice<Vecd>(particles, names.e_edge_recon_real);
        syncVariableToDevice<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToDevice<Real>(particles, names.joule_heat_edge_recon_real);

        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        Real *qc = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        for (size_t i = 0; i < n; ++i)
            qc[i] = q_i[i] + q_r[i];
        syncVariableToDevice<Real>(particles, names.joule_heat_edge_recon_complex);
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_D2_LOCAL_POLY_H
