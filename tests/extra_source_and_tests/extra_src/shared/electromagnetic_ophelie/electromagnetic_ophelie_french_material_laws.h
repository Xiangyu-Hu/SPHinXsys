#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "viscosity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Fixed-capacity, device-copyable temperature law.
 *
 * The digitized Jacoutot curves must be installed explicitly by the case
 * configuration. A one-point law is intentionally constant and is useful for
 * controlled regressions only; it is not a substitute for the paper's Fig. 2.
 */
class OphelieTemperatureLaw
{
  public:
    /** Thesis III.10 μ NUM has 23 nodes; keep headroom for endpoint hold points. */
    static constexpr size_t kMaxPoints = 32;

    std::array<Real, kMaxPoints> temperature_k_{};
    std::array<Real, kMaxPoints> value_{};
    size_t point_count_ = 0;
    bool logarithmic_value_ = false;

    void clear()
    {
        point_count_ = 0;
        logarithmic_value_ = false;
    }

    void setPoint(size_t index, Real temperature_k, Real value)
    {
        if (index >= kMaxPoints)
        {
            return;
        }
        temperature_k_[index] = temperature_k;
        value_[index] = value;
        point_count_ = std::max(point_count_, index + 1);
    }

    /** Doc156: scale all tabulated values (CEP1473K anchoring of thesis NUM shape). */
    void scaleValues(Real scale)
    {
        for (size_t i = 0; i < point_count_; ++i)
            value_[i] *= scale;
    }

    Real evaluate(Real temperature_k) const
    {
        if (point_count_ == 0)
        {
            return Real(0);
        }
        if (point_count_ == 1 || temperature_k <= temperature_k_[0])
        {
            return value_[0];
        }
        const size_t last = point_count_ - 1;
        if (temperature_k >= temperature_k_[last])
        {
            return value_[last];
        }
        for (size_t i = 0; i < last; ++i)
        {
            if (temperature_k <= temperature_k_[i + 1])
            {
                const Real t0 = temperature_k_[i];
                const Real t1 = temperature_k_[i + 1];
                const Real weight = (temperature_k - t0) / (t1 - t0 + TinyReal);
                if (logarithmic_value_ && value_[i] > TinyReal && value_[i + 1] > TinyReal)
                {
                    return std::exp((Real(1) - weight) * std::log(value_[i]) + weight * std::log(value_[i + 1]));
                }
                return (Real(1) - weight) * value_[i] + weight * value_[i + 1];
            }
        }
        return value_[last];
    }
};

struct OphelieFrenchGlassMaterialLaws
{
    OphelieTemperatureLaw rho;
    OphelieTemperatureLaw beta;
    OphelieTemperatureLaw mu;
    OphelieTemperatureLaw cp;
    OphelieTemperatureLaw k;
    OphelieTemperatureLaw sigma;
    bool paper_digitized_ = false;
    /** True when μ/cp/k come from CES 2008 Fig. 2 NUM digitization (not Table-1 constants). */
    bool fig2_thermal_rheology_digitized_ = false;
    const char *mu_source_ = "Jacoutot Table-1 @1473 K constant";
    const char *cp_source_ = "Jacoutot Table-1 @1473 K constant";
    const char *k_source_ = "Jacoutot Table-1 @1473 K constant";
    const char *sigma_source_ = "Jacoutot Table-1 @1473 K constant";
};

inline constexpr const char *kOphelieDynamicViscosityField = "OphelieDynamicViscosity";
inline constexpr const char *kOphelieHeatCapacityField = "OphelieHeatCapacity";

struct OphelieMaterialLawCheckResult
{
    bool passed = false;
    Real mu_1473 = 0.0;
    Real cp_1473 = 0.0;
    Real k_1473 = 0.0;
    Real sigma_1473 = 0.0;
    std::string message;
};

struct OphelieThermalRheologyStats
{
    /** Doc180: raw = curve evaluate before μ_cap; used = field read by viscous operators after clamp. */
    Real mu_raw_min = 0.0;
    Real mu_raw_max = 0.0;
    Real mu_raw_mean = 0.0;
    Real mu_used_min = 0.0;
    Real mu_used_max = 0.0;
    Real mu_used_mean = 0.0;
    /** Legacy aliases: after clamp helpers, equal to mu_used_*; after evaluate-only, equal to mu_raw_*. */
    Real mu_min = 0.0;
    Real mu_max = 0.0;
    Real mu_mean = 0.0;
    Real cp_min = 0.0;
    Real cp_max = 0.0;
    Real cp_mean = 0.0;
    Real k_min = 0.0;
    Real k_max = 0.0;
    Real k_mean = 0.0;
    Real mu_cap = 0.0;
    size_t n_mu_capped = 0;
    Real vol_frac_mu_capped = 0.0;
    std::string mu_stats_stage = "unset"; /**< e.g. after_evaluate_raw / after_cap_used */
};

/**
 * Exact Table-1 point from Jacoutot et al. (2008), T=1473 K.
 *
 * This is deliberately a constant regression law. Production French cases
 * should prefer thesis σ(T) laws (Fig. III.12 / IV.10) via
 * makeJacoutotThesisSigmaTemperatureLaw* and set paper_digitized_=true.
 *
 * Note: CEP paper Figure 2 is the coupling flowchart, not σ(T).
 */
inline OphelieFrenchGlassMaterialLaws makeJacoutot1473KReferenceMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws;
    laws.rho.setPoint(0, Real(1473), Real(2750));
    laws.beta.setPoint(0, Real(1473), Real(1.0e-5));
    laws.mu.setPoint(0, Real(1473), Real(4));
    laws.cp.setPoint(0, Real(1473), Real(1150));
    laws.k.setPoint(0, Real(1473), Real(4));
    laws.sigma.setPoint(0, Real(1473), Real(16));
    laws.mu.logarithmic_value_ = true;
    laws.sigma.logarithmic_value_ = true;
    laws.mu_source_ = "Jacoutot Table-1 @1473 K constant";
    laws.cp_source_ = "Jacoutot Table-1 @1473 K constant";
    laws.k_source_ = "Jacoutot Table-1 @1473 K constant";
    laws.sigma_source_ = "Jacoutot Table-1 @1473 K constant";
    return laws;
}

/**
 * OPEN_PARAMETER: visual digitization of Jacoutot et al., Chem. Eng. Sci. 63 (2008) Fig. 2 NUM
 * curves for μ(T), cp(T), k(T). Anchored at Table-1 values at 1473 K (μ=4, cp=1150, k=4).
 * Valid range for interpolation: [300, 1800] K; outside that range the table clamps.
 * Not an author-published closed-form; replace if an authoritative NUM fit becomes available.
 */
inline OphelieTemperatureLaw makeJacoutotCes2008Fig2ViscosityTemperatureLaw()
{
    OphelieTemperatureLaw mu;
    mu.logarithmic_value_ = true;
    // log10 μ roughly flat ~5 until ~800 K, then drops to Table-1 μ(1473)=4 and ~1 at 1800 K.
    const Real t[] = {Real(300), Real(550), Real(800), Real(1050), Real(1200), Real(1300),
                      Real(1400), Real(1473), Real(1550), Real(1650), Real(1800)};
    const Real v[] = {Real(1.0e5), Real(1.0e5), Real(8.0e4), Real(2.0e3), Real(80), Real(20),
                      Real(8), Real(4), Real(2.5), Real(1.5), Real(1)};
    for (size_t i = 0; i < 11; ++i)
    {
        mu.setPoint(i, t[i], v[i]);
    }
    return mu;
}

inline OphelieTemperatureLaw makeJacoutotCes2008Fig2HeatCapacityTemperatureLaw()
{
    OphelieTemperatureLaw cp;
    // Table-1 range 750–1200; NUM rises from ~800 toward a high-T plateau near 1200.
    const Real t[] = {Real(300), Real(550), Real(800), Real(1050), Real(1200), Real(1300),
                      Real(1400), Real(1473), Real(1550), Real(1650), Real(1800)};
    const Real v[] = {Real(800), Real(900), Real(1050), Real(1120), Real(1140), Real(1145),
                      Real(1148), Real(1150), Real(1160), Real(1180), Real(1200)};
    for (size_t i = 0; i < 11; ++i)
    {
        cp.setPoint(i, t[i], v[i]);
    }
    return cp;
}

inline OphelieTemperatureLaw makeJacoutotCes2008Fig2ThermalConductivityTemperatureLaw()
{
    OphelieTemperatureLaw k;
    // Near-linear rise; Table-1 k(1473)=4, journal range ~1–5.
    const Real t[] = {Real(300), Real(550), Real(800), Real(1050), Real(1200), Real(1300),
                      Real(1400), Real(1473), Real(1550), Real(1650), Real(1800)};
    const Real v[] = {Real(1.0), Real(1.5), Real(2.0), Real(2.6), Real(3.0), Real(3.3),
                      Real(3.6), Real(4.0), Real(4.2), Real(4.5), Real(5.0)};
    for (size_t i = 0; i < 11; ++i)
    {
        k.setPoint(i, t[i], v[i]);
    }
    return k;
}

inline void installJacoutotCes2008Fig2ThermalRheology(OphelieFrenchGlassMaterialLaws &laws)
{
    laws.mu = makeJacoutotCes2008Fig2ViscosityTemperatureLaw();
    laws.cp = makeJacoutotCes2008Fig2HeatCapacityTemperatureLaw();
    laws.k = makeJacoutotCes2008Fig2ThermalConductivityTemperatureLaw();
    laws.fig2_thermal_rheology_digitized_ = true;
    laws.mu_source_ = "CES 2008 Fig.2 NUM visual digitization (OPEN_PARAMETER), Table-1 anchor @1473 K";
    laws.cp_source_ = "CES 2008 Fig.2 NUM visual digitization (OPEN_PARAMETER), Table-1 anchor @1473 K";
    laws.k_source_ = "CES 2008 Fig.2 NUM visual digitization (OPEN_PARAMETER), Table-1 anchor @1473 K";
}

/**
 * Load digitized T,value CSV (header row; first two columns used).
 * Outside-node temperatures are handled by OphelieTemperatureLaw::evaluate (endpoint hold, no silent slope extrapolate).
 * For μ, set logarithmic_value=true (linear in T / log value).
 * If first T > cold_endpoint_T and prepend_cold_endpoint, insert (cold_endpoint_T, first_value) as explicit endpoint-hold.
 */
inline bool loadOphelieTemperatureLawFromCsv(OphelieTemperatureLaw &law, const std::string &path, bool logarithmic_value,
                                             std::string &message, bool prepend_cold_endpoint = false,
                                             Real cold_endpoint_T = Real(303))
{
    law.clear();
    law.logarithmic_value_ = logarithmic_value;
    std::ifstream in(path);
    if (!in)
    {
        message = "cannot_open:" + path;
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        message = "empty:" + path;
        return false;
    }
    struct Node
    {
        Real T;
        Real v;
    };
    std::vector<Node> nodes;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        std::string a, b;
        if (!std::getline(ss, a, ',') || !std::getline(ss, b, ','))
            continue;
        const Real T = static_cast<Real>(std::stod(a));
        const Real v = static_cast<Real>(std::stod(b));
        if (!(std::isfinite(T) && std::isfinite(v) && v > TinyReal))
            continue;
        nodes.push_back({T, v});
    }
    if (nodes.empty())
    {
        message = "no_nodes:" + path;
        return false;
    }
    std::sort(nodes.begin(), nodes.end(), [](const Node &x, const Node &y) { return x.T < y.T; });
    if (prepend_cold_endpoint && nodes.front().T > cold_endpoint_T + Real(1))
    {
        nodes.insert(nodes.begin(), Node{cold_endpoint_T, nodes.front().v});
    }
    if (nodes.size() > OphelieTemperatureLaw::kMaxPoints)
    {
        message = "too_many_nodes:" + path + " n=" + std::to_string(nodes.size());
        return false;
    }
    for (size_t i = 0; i < nodes.size(); ++i)
        law.setPoint(i, nodes[i].T, nodes[i].v);
    message = "ok n=" + std::to_string(nodes.size()) + " T=[" + std::to_string(nodes.front().T) + "," +
              std::to_string(nodes.back().T) + "]";
    return true;
}

struct ThesisMaterialCurveLoadReport
{
    bool ok = false;
    std::string k_path;
    std::string cp_path;
    std::string mu_path;
    std::string k_msg;
    std::string cp_msg;
    std::string mu_msg;
    std::string endpoint_policy =
        "endpoint_hold: evaluate clamps outside knot range; mu prepends T=303K with first-node value (curve starts ~502K)";
    std::string interpolation =
        "k,cp: linear in T/value; mu: linear in T/log(value) via OphelieTemperatureLaw::logarithmic_value_";
    Real k_1473 = 0;
    Real cp_1473 = 0;
    Real mu_1473 = 0;
    Real k_303 = 0;
    Real cp_303 = 0;
    Real mu_303 = 0;
    /** Doc156 CEP1473K anchors (applied after raw thesis load when requested). */
    bool cep_anchored = false;
    Real k_1473_raw = 0;
    Real cp_1473_raw = 0;
    Real mu_1473_raw = 0;
    Real s_k = 1;
    Real s_cp = 1;
    Real s_mu = 1;
    Real k_800 = 0;
    Real cp_800 = 0;
    Real mu_800 = 0;
    Real k_1200 = 0;
    Real cp_1200 = 0;
    Real mu_1200 = 0;
    Real k_1600 = 0;
    Real cp_1600 = 0;
    Real mu_1600 = 0;
};

/** Doc139 B: thesis Fig III.8=k, III.9=cp, III.10=μ digitized NUM (correct figure IDs). */
inline ThesisMaterialCurveLoadReport installThesisIII8910ThermalRheologyFromDigitizedCsv(
    OphelieFrenchGlassMaterialLaws &laws, const std::string &curves_dir,
    bool apply_cep1473_anchors = false, Real mu_anchor = Real(4), Real cp_anchor = Real(1150),
    Real k_anchor = Real(4))
{
    ThesisMaterialCurveLoadReport r;
    r.k_path = curves_dir + "/thesis_III8_k_T_digitized.csv";
    r.cp_path = curves_dir + "/thesis_III9_cp_T_digitized.csv";
    r.mu_path = curves_dir + "/thesis_III10_mu_T_digitized.csv";
    const bool k_ok = loadOphelieTemperatureLawFromCsv(laws.k, r.k_path, false, r.k_msg, true, Real(303));
    const bool cp_ok = loadOphelieTemperatureLawFromCsv(laws.cp, r.cp_path, false, r.cp_msg, true, Real(303));
    const bool mu_ok = loadOphelieTemperatureLawFromCsv(laws.mu, r.mu_path, true, r.mu_msg, true, Real(303));
    r.ok = k_ok && cp_ok && mu_ok;
    if (r.ok)
    {
        laws.fig2_thermal_rheology_digitized_ = false;
        laws.paper_digitized_ = true;
        laws.k_source_ = "thesis Fig.III.8 NUM digitized (french_reference_data_2026-10-15)";
        laws.cp_source_ = "thesis Fig.III.9 NUM digitized (french_reference_data_2026-10-15)";
        laws.mu_source_ = "thesis Fig.III.10 NUM digitized (french_reference_data_2026-10-15); log-interp; endpoint hold";
        r.k_1473_raw = laws.k.evaluate(Real(1473));
        r.cp_1473_raw = laws.cp.evaluate(Real(1473));
        r.mu_1473_raw = laws.mu.evaluate(Real(1473));
        r.k_1473 = r.k_1473_raw;
        r.cp_1473 = r.cp_1473_raw;
        r.mu_1473 = r.mu_1473_raw;
        if (apply_cep1473_anchors)
        {
            if (!(r.k_1473_raw > TinyReal && r.cp_1473_raw > TinyReal && r.mu_1473_raw > TinyReal))
            {
                r.ok = false;
                r.k_msg = "cep_anchor_failed:nonpositive_thesis_1473";
                return r;
            }
            r.s_k = k_anchor / r.k_1473_raw;
            r.s_cp = cp_anchor / r.cp_1473_raw;
            r.s_mu = mu_anchor / r.mu_1473_raw;
            laws.k.scaleValues(r.s_k);
            laws.cp.scaleValues(r.s_cp);
            laws.mu.scaleValues(r.s_mu);
            r.cep_anchored = true;
            laws.k_source_ =
                "thesis Fig.III.8 NUM digitized + Doc156 CEP1473K anchor scale s_k";
            laws.cp_source_ =
                "thesis Fig.III.9 NUM digitized + Doc156 CEP1473K anchor scale s_cp";
            laws.mu_source_ =
                "thesis Fig.III.10 NUM digitized + Doc156 CEP1473K anchor scale s_mu; log-interp";
            r.k_1473 = laws.k.evaluate(Real(1473));
            r.cp_1473 = laws.cp.evaluate(Real(1473));
            r.mu_1473 = laws.mu.evaluate(Real(1473));
        }
        r.k_303 = laws.k.evaluate(Real(303));
        r.cp_303 = laws.cp.evaluate(Real(303));
        r.mu_303 = laws.mu.evaluate(Real(303));
        r.k_800 = laws.k.evaluate(Real(800));
        r.cp_800 = laws.cp.evaluate(Real(800));
        r.mu_800 = laws.mu.evaluate(Real(800));
        r.k_1200 = laws.k.evaluate(Real(1200));
        r.cp_1200 = laws.cp.evaluate(Real(1200));
        r.mu_1200 = laws.mu.evaluate(Real(1200));
        r.k_1600 = laws.k.evaluate(Real(1600));
        r.cp_1600 = laws.cp.evaluate(Real(1600));
        r.mu_1600 = laws.mu.evaluate(Real(1600));
    }
    return r;
}

/**
 * Particle μ field for ViscousForceCK: harmonic mean of OphelieDynamicViscosity (not body-mean).
 * Bind after the field is registered. Inherits Viscosity so getMaterialProperty<Viscosity> still works for diagnostics.
 */
class ParticleViscosity : public Viscosity
{
    DiscreteVariable<Real> *dv_mu_ = nullptr;

  public:
    explicit ParticleViscosity(Real mu0) : Viscosity(mu0) {};
    virtual ~ParticleViscosity() {};

    void bindParticleViscosityField(BaseParticles &particles,
                                    const std::string &name = kOphelieDynamicViscosityField)
    {
        dv_mu_ = particles.getVariableByName<Real>(name);
    }

    DiscreteVariable<Real> *viscosityField() const { return dv_mu_; }

    typedef PairGeomAverageVariable<Real> InterParticleViscosity;
    typedef ParameterVariable<Real> OneSideViscosity;

    template <typename ExecutionPolicy>
    PairGeomAverageVariable<Real> getInterParticleViscosity(const ExecutionPolicy &ex_policy,
                                                            ParticleViscosity &other_viscosity)
    {
        return PairGeomAverageVariable<Real>(dv_mu_->DelegatedData(ex_policy),
                                             other_viscosity.dv_mu_->DelegatedData(ex_policy));
    }

    template <typename ExecutionPolicy>
    ParameterVariable<Real> getOneSideViscosity(const ExecutionPolicy &ex_policy)
    {
        return ParameterVariable<Real>(dv_mu_->DelegatedData(ex_policy));
    }
};

/** Host check: PairGeomAverageVariable of particle μ equals 2 μi μj/(μi+μj). */
inline bool hostVerifyParticleViscosityHarmonicMean(BaseParticles &particles, Real &max_rel_err, size_t max_pairs = 64)
{
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    const size_t n = particles.TotalRealParticles();
    PairGeomAverageVariable<Real> harm(mu, mu);
    max_rel_err = 0;
    size_t checked = 0;
    for (size_t i = 0; i < n && checked < max_pairs; i += std::max<size_t>(1, n / max_pairs))
    {
        const size_t j = (i + n / 3) % n;
        const Real expect = Real(2) * mu[i] * mu[j] / (mu[i] + mu[j] + TinyReal);
        const Real got = harm(i, j);
        const Real err = std::abs(got - expect) / (std::abs(expect) + TinyReal);
        max_rel_err = std::max(max_rel_err, err);
        ++checked;
    }
    return max_rel_err < Real(1.0e-12) || (sizeof(Real) == 4 && max_rel_err < Real(1.0e-5));
}

struct LocalViscosityForceEvidence
{
    bool harmonic_mean_ok = false;
    Real harmonic_max_rel_err = 0;
    Real mu_mean = 0;
    Real mu_max = 0;
    Real mu_min = 0;
    Real viscous_force_rms = 0;
    Real viscous_force_rms_high_mu = 0; /**< particles with mu > 0.5*(mean+max) */
    Real viscous_force_rms_low_mu = 0;
    size_t n_high_mu = 0;
    size_t n_low_mu = 0;
    bool local_mu_spatially_active = false; /**< mu_max/mu_mean > 1.05 */
    bool force_contrast_supports_local_mu = false;
};

inline OphelieMaterialLawCheckResult verifyFrenchMaterialLawsAtReferenceTemperatures(
    const OphelieFrenchGlassMaterialLaws &laws, Real rel_tol = Real(0.05))
{
    OphelieMaterialLawCheckResult result;
    result.mu_1473 = laws.mu.evaluate(Real(1473));
    result.cp_1473 = laws.cp.evaluate(Real(1473));
    result.k_1473 = laws.k.evaluate(Real(1473));
    result.sigma_1473 = laws.sigma.evaluate(Real(1473));
    const Real mu_lo = laws.mu.evaluate(Real(900));
    const Real mu_hi = laws.mu.evaluate(Real(1700));
    const bool table1_ok = std::abs(result.mu_1473 - Real(4)) <= rel_tol * Real(4) &&
                           std::abs(result.cp_1473 - Real(1150)) <= rel_tol * Real(1150) &&
                           std::abs(result.k_1473 - Real(4)) <= rel_tol * Real(4);
    const bool monotonic_mu = mu_lo + TinyReal >= result.mu_1473 && result.mu_1473 + TinyReal >= mu_hi;
    const bool finite_ok = std::isfinite(result.mu_1473) && std::isfinite(result.cp_1473) &&
                           std::isfinite(result.k_1473) && std::isfinite(result.sigma_1473) &&
                           result.mu_1473 > TinyReal && result.cp_1473 > TinyReal && result.k_1473 > TinyReal &&
                           result.sigma_1473 > TinyReal;
    result.passed = table1_ok && monotonic_mu && finite_ok;
    result.message = result.passed ? "material_laws_ok" : "material_laws_failed_table1_or_mu_trend";
    return result;
}

inline void registerOphelieFrenchThermalRheologyFields(BaseParticles &particles, Real initial_mu, Real initial_cp,
                                                       Real initial_k)
{
    particles.registerStateVariable<Real>(kOphelieDynamicViscosityField, initial_mu);
    particles.registerStateVariable<Real>(kOphelieHeatCapacityField, initial_cp);
    const size_t n = particles.TotalRealParticles();
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    Real *cp = particles.getVariableDataByName<Real>(kOphelieHeatCapacityField);
    // Conductivity is registered by registerOphelieThermalDiffusionAuxFields (same string name).
    Real *k = particles.getVariableDataByName<Real>("OphelieThermalConductivity");
    for (size_t i = 0; i < n; ++i)
    {
        mu[i] = initial_mu;
        cp[i] = initial_cp;
        k[i] = initial_k;
    }
}

inline void hostScalarFieldStats(const Real *field, const Real *vol, size_t n, Real &vmin, Real &vmax, Real &vmean)
{
    vmin = field[0];
    vmax = field[0];
    Real vol_sum = 0.0;
    Real weighted = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        vmin = std::min(vmin, field[i]);
        vmax = std::max(vmax, field[i]);
        vol_sum += vol[i];
        weighted += field[i] * vol[i];
    }
    vmean = weighted / (vol_sum + TinyReal);
}

inline LocalViscosityForceEvidence hostCollectLocalViscosityForceEvidence(BaseParticles &particles)
{
    LocalViscosityForceEvidence e;
    e.harmonic_mean_ok = hostVerifyParticleViscosityHarmonicMean(particles, e.harmonic_max_rel_err);
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    const Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *fv = particles.getVariableDataByName<Vecd>("ViscousForce");
    const size_t n = particles.TotalRealParticles();
    hostScalarFieldStats(mu, vol, n, e.mu_min, e.mu_max, e.mu_mean);
    e.local_mu_spatially_active = e.mu_mean > TinyReal && (e.mu_max / e.mu_mean) > Real(1.05);
    const Real thr = Real(0.5) * (e.mu_mean + e.mu_max);
    Real s2 = 0, s2h = 0, s2l = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real f2 = fv[i].squaredNorm();
        s2 += f2;
        if (mu[i] >= thr)
        {
            s2h += f2;
            ++e.n_high_mu;
        }
        else
        {
            s2l += f2;
            ++e.n_low_mu;
        }
    }
    e.viscous_force_rms = std::sqrt(s2 / Real(n));
    e.viscous_force_rms_high_mu = e.n_high_mu ? std::sqrt(s2h / Real(e.n_high_mu)) : 0;
    e.viscous_force_rms_low_mu = e.n_low_mu ? std::sqrt(s2l / Real(e.n_low_mu)) : 0;
    e.force_contrast_supports_local_mu =
        e.local_mu_spatially_active && e.n_high_mu > 0 && e.n_low_mu > 0 &&
        e.viscous_force_rms_high_mu + TinyReal >= Real(0.5) * e.viscous_force_rms_low_mu;
    return e;
}

/**
 * Evaluate μ/cp/k from T onto particle fields. When `constant_properties` is true, all particles
 * receive the law at `t_reference` (regression mode). Returns volume-mean μ for Viscosity refresh.
 */
inline Real hostUpdateFrenchThermalRheologyFromTemperature(BaseParticles &particles,
                                                           const OphelieFrenchGlassMaterialLaws &laws,
                                                           bool constant_properties, Real t_reference,
                                                           OphelieThermalRheologyStats &stats,
                                                           const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    Real *cp = particles.getVariableDataByName<Real>(kOphelieHeatCapacityField);
    Real *k = particles.getVariableDataByName<Real>("OphelieThermalConductivity");
    for (size_t i = 0; i < n; ++i)
    {
        const Real t = constant_properties ? t_reference : temperature[i];
        mu[i] = std::max(laws.mu.evaluate(t), TinyReal);
        cp[i] = std::max(laws.cp.evaluate(t), TinyReal);
        k[i] = std::max(laws.k.evaluate(t), TinyReal);
    }
    hostScalarFieldStats(mu, vol, n, stats.mu_raw_min, stats.mu_raw_max, stats.mu_raw_mean);
    stats.mu_min = stats.mu_raw_min;
    stats.mu_max = stats.mu_raw_max;
    stats.mu_mean = stats.mu_raw_mean;
    stats.mu_stats_stage = "after_evaluate_raw";
    hostScalarFieldStats(cp, vol, n, stats.cp_min, stats.cp_max, stats.cp_mean);
    hostScalarFieldStats(k, vol, n, stats.k_min, stats.k_max, stats.k_mean);
    syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
    syncVariableToDevice<Real>(particles, kOphelieHeatCapacityField);
    syncVariableToDevice<Real>(particles, "OphelieThermalConductivity");
    return stats.mu_raw_mean;
}

/**
 * Doc180: clamp particle μ to mu_cap in-place, then fill mu_used_* from the field operators read.
 * Does not recompute mu_raw_* (caller must have evaluated raw first).
 */
inline void hostClampParticleViscosityCap(BaseParticles &particles, Real mu_cap, OphelieThermalRheologyStats &stats)
{
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    size_t n_cap = 0;
    Real vol_cap = 0;
    Real vol_sum = 0;
    const Real cap = std::max(mu_cap, TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        vol_sum += vol[i];
        if (mu[i] > cap)
        {
            mu[i] = cap;
            ++n_cap;
            vol_cap += vol[i];
        }
    }
    hostScalarFieldStats(mu, vol, n, stats.mu_used_min, stats.mu_used_max, stats.mu_used_mean);
    stats.mu_min = stats.mu_used_min;
    stats.mu_max = stats.mu_used_max;
    stats.mu_mean = stats.mu_used_mean;
    stats.mu_cap = cap;
    stats.n_mu_capped = n_cap;
    stats.vol_frac_mu_capped = vol_cap / (vol_sum + TinyReal);
    stats.mu_stats_stage = "after_cap_used";
    syncVariableToDevice<Real>(particles, kOphelieDynamicViscosityField);
}

/**
 * Exact segment integral of a piecewise-linear (or clamp-constant) cp law.
 * Matches OphelieTemperatureLaw::evaluate, including out-of-range clamping.
 * Logarithmic laws use a local Simpson rule on each knot-split interval.
 */
inline Real ophelieIntegrateHeatCapacity(const OphelieTemperatureLaw &cp, Real t_from, Real t_to)
{
    if (cp.point_count_ == 0 || !(std::isfinite(t_from) && std::isfinite(t_to)))
    {
        return Real(0);
    }
    if (std::abs(t_to - t_from) <= TinyReal)
    {
        return Real(0);
    }
    const bool negate = t_to < t_from;
    const Real a = negate ? t_to : t_from;
    const Real b = negate ? t_from : t_to;

    Real splits[OphelieTemperatureLaw::kMaxPoints + 2];
    size_t n_splits = 0;
    splits[n_splits++] = a;
    for (size_t i = 0; i < cp.point_count_; ++i)
    {
        const Real tk = cp.temperature_k_[i];
        if (tk > a + TinyReal && tk < b - TinyReal)
        {
            splits[n_splits++] = tk;
        }
    }
    splits[n_splits++] = b;
    std::sort(splits, splits + n_splits);
    size_t unique = 0;
    for (size_t i = 0; i < n_splits; ++i)
    {
        if (unique == 0 || splits[i] > splits[unique - 1] + TinyReal)
        {
            splits[unique++] = splits[i];
        }
    }

    Real sum = Real(0);
    for (size_t i = 0; i + 1 < unique; ++i)
    {
        const Real x0 = splits[i];
        const Real x1 = splits[i + 1];
        const Real dx = x1 - x0;
        if (dx <= TinyReal)
        {
            continue;
        }
        if (!cp.logarithmic_value_)
        {
            // Linear (or constant-clamp) on each knot-bounded interval → trapezoid exact.
            sum += Real(0.5) * (cp.evaluate(x0) + cp.evaluate(x1)) * dx;
        }
        else
        {
            const Real xm = Real(0.5) * (x0 + x1);
            sum += dx / Real(6) * (cp.evaluate(x0) + Real(4) * cp.evaluate(xm) + cp.evaluate(x1));
        }
    }
    return negate ? -sum : sum;
}

/** Specific enthalpy h(T) - h(T_ref) = ∫_{T_ref}^T cp(θ) dθ [J/(kg·K)·K = J/kg]. */
inline Real ophelieSpecificEnthalpy(const OphelieTemperatureLaw &cp, Real temperature_k, Real t_ref)
{
    return ophelieIntegrateHeatCapacity(cp, t_ref, temperature_k);
}

/**
 * Constant-cp sensible energy relative to t_ref (doc 42 ledger):
 *   H = Σ ρ cp (T − t_ref) V
 * Accumulates in double. Do not use float Σ ρ cp T V: at T~1473 K that product
 * is ~1e8 J and WCSPH volume noise swamps ΔH. Matches hostThermalEnthalpyEnergy
 * for a one-point (constant) cp law.
 */
inline Real hostThermalEnergyConstantCp(BaseParticles &particles, Real rho, Real cp, Real t_ref,
                                        const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *t = particles.getVariableDataByName<Real>(temperature_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const double rho_d = static_cast<double>(rho);
    const double cp_d = static_cast<double>(cp);
    const double t_ref_d = static_cast<double>(t_ref);
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e += rho_d * cp_d * (static_cast<double>(t[i]) - t_ref_d) * static_cast<double>(vol[i]);
    }
    return static_cast<Real>(e);
}

/** Doc 45/47/66/68 F07: mass/volume/density conventions used by the energy path. */
struct OphelieThermalMassVolumeAudit
{
    size_t n = 0;
    double sum_mass = 0.0;
    double sum_vol = 0.0;
    double sum_rho0_vol = 0.0;
    double sum_rho_vol = 0.0;
    double rho_min = 0.0;
    double rho_max = 0.0;
    double rho_mean = 0.0;
    /** Particlewise Mass − ρV (same host sync). */
    double max_abs_mass_minus_rho_vol = 0.0;
    double rms_mass_minus_rho_vol = 0.0;
    double rel_l2_mass_vs_rho_vol = 0.0;
    size_t n_mass_rho_vol_checked = 0;
    bool has_mass = false;
    bool has_density = false;
};

inline OphelieThermalMassVolumeAudit hostThermalMassVolumeAudit(BaseParticles &particles, Real rho0)
{
    OphelieThermalMassVolumeAudit a;
    a.n = particles.TotalRealParticles();
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    DiscreteVariable<Real> *mass_var = particles.getVariableByName<Real>("Mass");
    a.has_mass = mass_var != nullptr;
    const Real *mass = nullptr;
    if (a.has_mass)
    {
        syncVariableToHost<Real>(particles, "Mass");
        mass = particles.getVariableDataByName<Real>("Mass");
    }
    DiscreteVariable<Real> *rho_var = particles.getVariableByName<Real>("Density");
    a.has_density = rho_var != nullptr;
    const Real *rho = nullptr;
    if (a.has_density)
    {
        syncVariableToHost<Real>(particles, "Density");
        rho = particles.getVariableDataByName<Real>("Density");
    }
    const double rho0_d = static_cast<double>(rho0);
    double rho_sum = 0.0;
    double sum_sq_mrv = 0.0;
    double sum_sq_mass = 0.0;
    for (size_t i = 0; i < a.n; ++i)
    {
        const double v = static_cast<double>(vol[i]);
        a.sum_vol += v;
        a.sum_rho0_vol += rho0_d * v;
        if (mass != nullptr)
        {
            a.sum_mass += static_cast<double>(mass[i]);
        }
        if (rho != nullptr)
        {
            const double ri = static_cast<double>(rho[i]);
            a.sum_rho_vol += ri * v;
            rho_sum += ri;
            if (i == 0)
            {
                a.rho_min = ri;
                a.rho_max = ri;
            }
            else
            {
                a.rho_min = std::min(a.rho_min, ri);
                a.rho_max = std::max(a.rho_max, ri);
            }
            if (mass != nullptr)
            {
                const double mi = static_cast<double>(mass[i]);
                const double diff = mi - ri * v;
                a.max_abs_mass_minus_rho_vol = std::max(a.max_abs_mass_minus_rho_vol, std::abs(diff));
                sum_sq_mrv += diff * diff;
                sum_sq_mass += mi * mi;
                ++a.n_mass_rho_vol_checked;
            }
        }
    }
    if (a.n > 0 && a.has_density)
    {
        a.rho_mean = rho_sum / static_cast<double>(a.n);
    }
    if (a.n_mass_rho_vol_checked > 0)
    {
        a.rms_mass_minus_rho_vol = std::sqrt(sum_sq_mrv / static_cast<double>(a.n_mass_rho_vol_checked));
        a.rel_l2_mass_vs_rho_vol =
            std::sqrt(sum_sq_mrv) / (std::sqrt(sum_sq_mass) + std::numeric_limits<double>::epsilon());
    }
    return a;
}

/**
 * Doc 66/68 diagnostic only: H_m = Σ m cp (T − t_ref). Production ledger still uses
 * hostThermalEnergyConstantCp(rho0,V). Do not swap production H without matching Q/BC.
 */
inline Real hostThermalEnergyConstantCpFromMass(BaseParticles &particles, Real cp, Real t_ref,
                                                const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "Mass");
    const Real *t = particles.getVariableDataByName<Real>(temperature_field);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t n = particles.TotalRealParticles();
    const double cp_d = static_cast<double>(cp);
    const double t_ref_d = static_cast<double>(t_ref);
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e += static_cast<double>(mass[i]) * cp_d * (static_cast<double>(t[i]) - t_ref_d);
    }
    return static_cast<Real>(e);
}

/** Doc 69 P1: same-state F07 snapshot row (CSV). */
inline void writeFrenchF07StateAuditCsvHeader(const std::string &path)
{
    std::ofstream out(path);
    out << std::setprecision(17)
        << "stage,advection_steps,acoustic_steps,t,n_particles,sum_mass,sum_vol,sum_rho0_vol,sum_rho_vol,"
           "rho_min,rho_max,rho_mean,max_abs_mass_minus_rho_vol,rms_mass_minus_rho_vol,rel_l2_mass_vs_rho_vol,"
           "H_rho0V,H_mass,E_Q_applied,E_boundary_out,schema\n";
}

inline void appendFrenchF07StateAuditCsv(const std::string &path, const char *stage, size_t advection_steps,
                                         size_t acoustic_steps, Real t, const OphelieThermalMassVolumeAudit &mv,
                                         Real h_rho0v, Real h_mass, Real e_q, Real e_b)
{
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17) << stage << "," << advection_steps << "," << acoustic_steps << "," << t << ","
        << mv.n << "," << mv.sum_mass << "," << mv.sum_vol << "," << mv.sum_rho0_vol << "," << mv.sum_rho_vol << ","
        << mv.rho_min << "," << mv.rho_max << "," << mv.rho_mean << "," << mv.max_abs_mass_minus_rho_vol << ","
        << mv.rms_mass_minus_rho_vol << "," << mv.rel_l2_mass_vs_rho_vol << "," << h_rho0v << "," << h_mass << ","
        << e_q << "," << e_b << ",1\n";
}

/**
 * Melt thermal energy consistent with T += Q dt/(ρ cp): H = Σ ρ V_i h(T_i),
 * h = ∫_{t_ref}^{T} cp(θ) dθ. Uses double accumulation.
 */
inline Real hostThermalEnthalpyEnergy(BaseParticles &particles, Real rho, const OphelieTemperatureLaw &cp_law,
                                      Real t_ref, const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *t = particles.getVariableDataByName<Real>(temperature_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e += static_cast<double>(rho) * static_cast<double>(vol[i]) *
             static_cast<double>(ophelieSpecificEnthalpy(cp_law, t[i], t_ref));
    }
    return static_cast<Real>(e);
}

/**
 * Host regression for P0-05: constant-cp adiabatic identity, piecewise Fig.2 integral,
 * and demonstration that Σ ρ cp(T) T V is not the enthalpy change.
 */
inline bool verifyOphelieEnthalpyEnergyAccounting(std::string &message)
{
    OphelieTemperatureLaw cp_const;
    cp_const.setPoint(0, Real(1473), Real(1150));
    const Real t_ref = Real(1473);
    const Real dT = Real(100);
    const Real h_const = ophelieSpecificEnthalpy(cp_const, t_ref + dT, t_ref);
    if (!(std::abs(h_const - Real(1150) * dT) <= Real(1.0e-9) * Real(1150) * dT))
    {
        message = "constant_cp_enthalpy_mismatch";
        return false;
    }

    const OphelieTemperatureLaw cp_var = makeJacoutotCes2008Fig2HeatCapacityTemperatureLaw();
    const Real t_a = Real(800);
    const Real t_b = Real(1600);
    const Real h_ab = ophelieSpecificEnthalpy(cp_var, t_b, t_a);
    // Independent trapezoid on a fine uniform grid (non-log law).
    constexpr int n_samples = 800;
    Real h_num = Real(0);
    for (int i = 0; i < n_samples; ++i)
    {
        const Real x0 = t_a + (t_b - t_a) * (Real(i) / Real(n_samples));
        const Real x1 = t_a + (t_b - t_a) * (Real(i + 1) / Real(n_samples));
        h_num += Real(0.5) * (cp_var.evaluate(x0) + cp_var.evaluate(x1)) * (x1 - x0);
    }
    if (!(std::abs(h_ab - h_num) <= Real(1.0e-6) * std::max(std::abs(h_num), Real(1))))
    {
        message = "piecewise_cp_integral_mismatch";
        return false;
    }

    // Explicit T += Q dt/(ρ cp) deposits a left Riemann sum; enthalpy is the exact
    // segment integral. Defect must shrink under refinement (doc 23 P0-05).
    const Real rho = Real(2750);
    const Real vol = Real(1.0e-3);
    const Real t0 = Real(1473);
    const Real t1 = Real(1523);
    const Real h_exact = rho * vol * ophelieSpecificEnthalpy(cp_var, t1, t0);
    auto adiabatic_defect = [&](int n_steps) -> Real {
        Real t = t0;
        const Real dt_step = (t1 - t0) / Real(n_steps);
        double deposited = 0.0;
        for (int s = 0; s < n_steps; ++s)
        {
            const Real cp_now = cp_var.evaluate(t);
            const Real q = rho * cp_now * dt_step; // synthetic dt = 1
            deposited += static_cast<double>(q) * static_cast<double>(vol);
            t += dt_step;
        }
        return std::abs(static_cast<Real>(deposited) - h_exact) / (std::abs(h_exact) + TinyReal);
    };
    const Real defect_coarse = adiabatic_defect(10);
    const Real defect_fine = adiabatic_defect(400);
    if (!(defect_fine + TinyReal < defect_coarse))
    {
        message = "adiabatic_defect_not_refined";
        return false;
    }
    if (!(defect_fine < Real(5.0e-4)))
    {
        message = "adiabatic_defect_too_large";
        return false;
    }

    const Real wrong_delta = cp_var.evaluate(t_b) * t_b - cp_var.evaluate(t_a) * t_a;
    if (!(std::abs(wrong_delta - h_ab) > Real(1)))
    {
        message = "cpT_bookkeeping_unexpectedly_matches_enthalpy";
        return false;
    }

    message = "enthalpy_ok";
    return true;
}

/**
 * Doc 23 P1-01: Viscosity::getInterParticleViscosity / getOneSideViscosity copy μ by value.
 * setReferenceViscosity alone does not update an already-built PairGeomAverageFixed / ParameterFixed.
 */
inline bool verifyFrenchViscosityFunctorCachesByValue(std::string &message)
{
    Viscosity viscosity(Real(2));
    execution::ParallelPolicy policy;
    auto pair_mu = viscosity.getInterParticleViscosity(policy, viscosity);
    auto one_side = viscosity.getOneSideViscosity(policy);
    const Real before_pair = pair_mu(0, 0);
    const Real before_one = one_side(0);
    viscosity.setReferenceViscosity(Real(8));
    const Real after_pair = pair_mu(0, 0);
    const Real after_one = one_side(0);
    if (!(std::abs(before_pair - Real(2)) <= Real(1.0e-12) && std::abs(before_one - Real(2)) <= Real(1.0e-12)))
    {
        message = "viscosity_functor_initial_value_wrong";
        return false;
    }
    if (!(std::abs(after_pair - before_pair) <= Real(1.0e-12) && std::abs(after_one - before_one) <= Real(1.0e-12)))
    {
        message = "viscosity_functor_unexpectedly_live";
        return false;
    }
    auto pair_mu_new = viscosity.getInterParticleViscosity(policy, viscosity);
    if (!(std::abs(pair_mu_new(0, 0) - Real(8)) <= Real(1.0e-12)))
    {
        message = "viscosity_functor_rebuild_failed";
        return false;
    }
    message = "viscosity_ck_cache_by_value_ok";
    return true;
}

/** Joule heating using per-particle cp(T): T += Q dt / (ρ cp_i). */
class ApplyJouleHeatWithParticleCpCK : public LocalDynamics
{
  public:
    ApplyJouleHeatWithParticleCpCK(SPHBody &sph_body, Real rho, Real t_initial,
                                   const std::string &q_field = "JouleHeat")
        : LocalDynamics(sph_body), inv_rho_(Real(1) / (rho + TinyReal)), t_initial_(t_initial),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_cp_(particles_->template getVariableByName<Real>(kOphelieHeatCapacityField)),
          dv_temperature_(particles_->template getVariableByName<Real>("Temperature")),
          dv_delta_t_(particles_->template getVariableByName<Real>("OphelieThermalDeltaT"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : inv_rho_(encloser.inv_rho_), t_initial_(encloser.t_initial_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)), cp_(encloser.dv_cp_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            // Doc 45 F01: accumulate heating on ΔT even for particle-cp path.
            const Real heating = q_[index_i] * dt * inv_rho_ / (cp_[index_i] + TinyReal);
            delta_t_[index_i] += heating;
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real inv_rho_;
        Real t_initial_;
        Real *q_;
        Real *cp_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real inv_rho_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_cp_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

/**
 * Legacy Σ ρ cp_i(T_i) T_i V_i — inconsistent with T += Q dt/(ρ cp) when cp varies.
 * Prefer hostThermalEnthalpyEnergy with the same cp law and t_ref.
 */
inline Real hostThermalEnergyWithParticleCp(BaseParticles &particles, Real rho,
                                            const std::string &temperature_field = "Temperature")
{
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, kOphelieHeatCapacityField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *t = particles.getVariableDataByName<Real>(temperature_field);
    const Real *cp = particles.getVariableDataByName<Real>(kOphelieHeatCapacityField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e += static_cast<double>(rho) * static_cast<double>(cp[i]) * static_cast<double>(t[i]) *
             static_cast<double>(vol[i]);
    }
    return static_cast<Real>(e);
}

/**
 * Provisional Arrhenius-like sigma(T) anchored at Jacoutot Table-1:
 *   sigma(1473 K)=16 S/m, Ea/R = 15000 K.
 *
 * This is NOT thesis Fig. III.12 / IV.10. It exists so Stage 3.2 can exercise
 * T→σ→EM→Q→T machinery before switching to the thesis closed-form laws.
 */
inline OphelieTemperatureLaw makeProvisionalJacoutotLikeSigmaTemperatureLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real t_ref = Real(1473);
    constexpr Real sigma_ref = Real(16);
    constexpr Real ea_over_r = Real(15000);
    const Real temperatures[5] = {Real(1273), Real(1373), Real(1473), Real(1573), Real(1673)};
    for (size_t i = 0; i < 5; ++i)
    {
        const Real t = temperatures[i];
        const Real value = sigma_ref * std::exp(ea_over_r * (Real(1) / t_ref - Real(1) / t));
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * Jacoutot thesis Fig. III.12 / law (III.6) [C_2], natural / static bath track:
 *   log10(sigma) = 3.7921 - 3179.8 / T
 * Sampled onto OphelieTemperatureLaw with log-space interpolation.
 * At 1473 K this gives ~43 S/m (Table-1's 16 is a representative point, not this law).
 *
 * GPT 2026-08-11: raw III.12 is sensitivity-only for French Natural CEP2008 reproduction.
 */
inline OphelieTemperatureLaw makeJacoutotThesisSigmaTemperatureLawIII12()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real a = Real(3.7921);
    constexpr Real b = Real(3179.8);
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1229), Real(1300),
                                   Real(1400), Real(1473), Real(1500), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        const Real value = std::pow(Real(10), a - b / t);
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * French Natural (CEP 2008) production σ(T):
 *   journal Table-1 anchor σ(1473K)=16, thesis III.12 temperature shape, clipped to [1e-6, 30].
 *   sigma_nat(T) = 16 * 10^(3179.8 * (1/1473 - 1/T))
 */
inline OphelieTemperatureLaw makeFrenchNaturalJournalAnchoredSigmaTemperatureLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real t_ref = Real(1473);
    constexpr Real sigma_ref = Real(16);
    constexpr Real b = Real(3179.8);
    constexpr Real sigma_min = Real(1.0e-6);
    constexpr Real sigma_max = Real(30);
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1550), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        const Real raw = sigma_ref * std::pow(Real(10), b * (Real(1) / t_ref - Real(1) / t));
        const Real value = std::min(sigma_max, std::max(sigma_min, raw));
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * Jacoutot thesis Fig. IV.10, stirred Uox2 + RuO2 glass [C_8]:
 *   log10(sigma) = 4.05726 - 3923.73 / T
 */
inline OphelieTemperatureLaw makeJacoutotThesisSigmaTemperatureLawIV10()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real a = Real(4.05726);
    constexpr Real b = Real(3923.73);
    const Real temperatures[11] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 11; ++i)
    {
        const Real t = temperatures[i];
        const Real value = std::pow(Real(10), a - b / t);
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

inline OphelieFrenchGlassMaterialLaws makeProvisionalJacoutotLikeMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeProvisionalJacoutotLikeSigmaTemperatureLaw();
    laws.paper_digitized_ = false;
    return laws;
}

inline OphelieFrenchGlassMaterialLaws makeJacoutotThesisIII12MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeJacoutotThesisSigmaTemperatureLawIII12();
    laws.paper_digitized_ = true;
    return laws;
}

/** CEP2008 Natural production material: journal-anchored σ(T) + CES Fig.2 μ/cp/k digitization. */
inline OphelieFrenchGlassMaterialLaws makeFrenchNaturalCep2008MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeFrenchNaturalJournalAnchoredSigmaTemperatureLaw();
    laws.sigma_source_ = "CEP2008 journal-anchored σ(T): 16*10^(3179.8*(1/1473-1/T)), clip [1e-6,30]";
    installJacoutotCes2008Fig2ThermalRheology(laws);
    laws.paper_digitized_ = true;
    return laws;
}

inline OphelieFrenchGlassMaterialLaws makeJacoutotThesisIV10MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeJacoutotThesisSigmaTemperatureLawIV10();
    laws.sigma_source_ = "thesis Fig.IV.10: log10(sigma)=4.05726-3923.73/T";
    installJacoutotCes2008Fig2ThermalRheology(laws);
    laws.paper_digitized_ = true;
    return laws;
}

/** Constant Table-1 thermal/rheology regression companion (σ unchanged from caller). */
inline void forceConstantTable1ThermalRheology(OphelieFrenchGlassMaterialLaws &laws)
{
    laws.mu.clear();
    laws.cp.clear();
    laws.k.clear();
    laws.mu.setPoint(0, Real(1473), Real(4));
    laws.cp.setPoint(0, Real(1473), Real(1150));
    laws.k.setPoint(0, Real(1473), Real(4));
    laws.mu.logarithmic_value_ = true;
    laws.fig2_thermal_rheology_digitized_ = false;
    laws.mu_source_ = "Jacoutot Table-1 @1473 K constant (regression)";
    laws.cp_source_ = "Jacoutot Table-1 @1473 K constant (regression)";
    laws.k_source_ = "Jacoutot Table-1 @1473 K constant (regression)";
}

inline constexpr const char *kOphelieSigmaRawTemperatureField = "OphelieSigmaRawTemperature";

inline void registerOphelieFrenchTemperatureMaterialFields(BaseParticles &particles, Real initial_sigma)
{
    particles.registerStateVariable<Real>(kOphelieSigmaRawTemperatureField, initial_sigma);
}

inline void hostOphelieSigmaFieldStats(BaseParticles &particles, const std::string &sigma_field, size_t n, Real &sigma_min,
                                       Real &sigma_max, Real &sigma_mean)
{
    syncVariableToHost<Real>(particles, sigma_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *sigma = particles.getVariableDataByName<Real>(sigma_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    sigma_min = sigma[0];
    sigma_max = sigma[0];
    Real vol_sum = 0.0;
    Real weighted = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sigma_min = std::min(sigma_min, sigma[i]);
        sigma_max = std::max(sigma_max, sigma[i]);
        vol_sum += vol[i];
        weighted += sigma[i] * vol[i];
    }
    sigma_mean = weighted / (vol_sum + TinyReal);
}

inline Real hostOphelieRelativeScalarFieldChange(BaseParticles &particles, const std::string &field_name,
                                                 const StdVec<Real> &previous, size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    const Real *field = particles.getVariableDataByName<Real>(field_name);
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real delta = field[i] - previous[i];
        numerator += delta * delta;
        denominator += field[i] * field[i];
    }
    return std::sqrt(numerator / (denominator + TinyReal));
}

inline void hostStoreScalarField(BaseParticles &particles, const std::string &field_name, StdVec<Real> &storage,
                                 size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    const Real *field = particles.getVariableDataByName<Real>(field_name);
    storage.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        storage[i] = field[i];
    }
}

/**
 * Applies sigma(T) with pointwise under-relaxation:
 * sigma^(n+1)=(1-alpha)sigma^n+alpha*sigma_raw(T).
 */
class UpdateOphelieGlassSigmaFromTemperatureCK : public LocalDynamics
{
  public:
    UpdateOphelieGlassSigmaFromTemperatureCK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
                                             const OphelieTemperatureLaw &sigma_law, Real relaxation_alpha,
                                             const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), sigma_law_(sigma_law),
          relaxation_alpha_(std::max(Real(0), std::min(relaxation_alpha, Real(1)))),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
          dv_sigma_raw_(particles_->template getVariableByName<Real>(kOphelieSigmaRawTemperatureField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_law_(encloser.sigma_law_), relaxation_alpha_(encloser.relaxation_alpha_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              sigma_raw_(encloser.dv_sigma_raw_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            const Real raw = sigma_law_.evaluate(temperature_[index_i]);
            sigma_raw_[index_i] = std::max(raw, TinyReal);
            sigma_[index_i] = (Real(1) - relaxation_alpha_) * sigma_[index_i] +
                              relaxation_alpha_ * sigma_raw_[index_i];
        }

      protected:
        OphelieTemperatureLaw sigma_law_;
        Real relaxation_alpha_;
        Real *temperature_;
        Real *sigma_;
        Real *sigma_raw_;
    };

  protected:
    OphelieTemperatureLaw sigma_law_;
    Real relaxation_alpha_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_sigma_raw_;
};

/** Diagnostic seed: only Temperature is set; DeltaT stays 0 so Joule closure remains heating-only. */
class SeedOphelieRadialTemperatureCK : public LocalDynamics
{
  public:
    SeedOphelieRadialTemperatureCK(SPHBody &sph_body, const Vecd &center, Real glass_radius, Real t0, Real delta_t,
                                   const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), center_(center), glass_radius_(std::max(glass_radius, TinyReal)), t0_(t0),
          delta_t_(delta_t), dv_pos_(particles_->template getVariableByName<Vecd>("Position")),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : center_(encloser.center_), glass_radius_(encloser.glass_radius_), t0_(encloser.t0_),
              delta_t_(encloser.delta_t_), pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            const Real dx = pos_[index_i][0] - center_[0];
            const Real dy = pos_[index_i][1] - center_[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real frac = std::min(r / glass_radius_, Real(1));
            temperature_[index_i] = t0_ + delta_t_ * frac;
        }

      protected:
        Vecd center_;
        Real glass_radius_;
        Real t0_;
        Real delta_t_;
        Vecd *pos_;
        Real *temperature_;
    };

  protected:
    Vecd center_;
    Real glass_radius_;
    Real t0_;
    Real delta_t_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_temperature_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H
