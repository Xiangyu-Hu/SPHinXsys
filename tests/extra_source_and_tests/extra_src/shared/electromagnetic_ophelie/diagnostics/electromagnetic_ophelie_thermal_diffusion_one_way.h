#ifndef ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H
#define ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_thermal_vtp.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Thermal BC mode for French one-way diffusion:
 * - DirichletRegression: cold-wall T=T0 (Stage 3.3 regression only)
 * - FrenchNaturalProduction: Robin side/bottom + free-surface convection/radiation (Stage 4.0)
 */
enum class OphelieThermalBoundaryMode
{
    None = 0,
    DirichletRegression = 1,
    FrenchNaturalProduction = 2
};

/** Isotropic conduction + optional thermal BCs (no EM feedback). */
struct OphelieThermalDiffusionOneWayOptions
{
    bool enable_diffusion = false;
    /** Legacy alias: Dirichlet regression shell (mutually exclusive with french-natural). */
    bool enable_cold_wall_dirichlet = false;
    bool enable_french_natural_bc = false;
    /** Shell width = factor * dp (match phi boundary convention). */
    Real boundary_width_factor = 1.5;
    Real pair_weight_regularization = Real(0.01);
    /** Natural production BC (journal sensitivity range; h_side default selected, not unique). */
    Real h_side = Real(300);     // W/(m^2 K)
    Real h_bottom = Real(35);    // W/(m^2 K)
    Real h_free = Real(20);      // W/(m^2 K)
    Real emissivity = Real(0.8); // first-version free-surface emissivity
    Real t_cool = Real(300);     // crucible coolant / wall ambient [K]
    Real t_ambient = Real(300);  // free-surface convection ambient [K]
    Real t_rad_ambient = Real(300); // radiative ambient [K]
    Real stefan_boltzmann = Real(5.670374419e-8);
    /**
     * Level-1 companion: free-surface crust conduction in series with surface conv+rad.
     * q = (T_melt - T_s) / (delta/k) where T_s solves surface balance at the outer face.
     * OPEN_PARAMETER thickness; not a resolved crust geometry.
     */
    bool enable_free_crust = false;
    Real free_crust_delta_m = Real(0.005);
    Real free_crust_k = Real(1.0);
    /** Filled by setup helpers: boundary_width_factor * dp. */
    Real boundary_shell_thickness = Real(0);
    /**
     * Doc 40 R03 / Doc 42: after multi-face tagging, scale each face weight so
     * Σ w_f V_i/δ matches the geometric face area of the tagging cylinder.
     * Same scales enter the T update and the loss audit. Not a silent h/skull tune;
     * not `--balance-heat-loss`. French literature uses continuum geometry; this only
     * corrects the SPH shell measure.
     */
    bool renormalize_boundary_area_to_geometry = false;
    Real area_scale_side = Real(1);
    Real area_scale_bottom = Real(1);
    Real area_scale_free = Real(1);
    /**
     * Doc207 B thermal wall model (OPHELIE-local; production default remains legacy-shell).
     * - legacy-shell: q = h(T_particle - T∞) on shell (no Tw reconstruction)
     * - wall-temperature: series gap q = k(Ti-Tw)/d = h(Tw-T∞)+… with d = geometric wall gap
     */
    std::string thermal_wall_model = "legacy-shell";
};

/** Doc207 B: true iff options request the wall-temperature series-gap candidate. */
inline bool ophelieThermalWallModelIsWallTemperature(const OphelieThermalDiffusionOneWayOptions &options)
{
    return options.thermal_wall_model == "wall-temperature" || options.thermal_wall_model == "wall_temperature" ||
           options.thermal_wall_model == "series-gap";
}

inline OphelieThermalBoundaryMode ophelieThermalBoundaryModeFromOptions(
    const OphelieThermalDiffusionOneWayOptions &options)
{
    if (options.enable_french_natural_bc)
    {
        return OphelieThermalBoundaryMode::FrenchNaturalProduction;
    }
    if (options.enable_cold_wall_dirichlet)
    {
        return OphelieThermalBoundaryMode::DirichletRegression;
    }
    return OphelieThermalBoundaryMode::None;
}

inline const char *ophelieThermalBoundaryModeName(OphelieThermalBoundaryMode mode)
{
    switch (mode)
    {
    case OphelieThermalBoundaryMode::DirichletRegression:
        return "dirichlet_regression";
    case OphelieThermalBoundaryMode::FrenchNaturalProduction:
        return "french_natural_production";
    default:
        return "none";
    }
}

/**
 * Free-surface heat flux with optional crust conduction (device-safe; no std / lambda).
 * Without crust: q = h(T-Ta) + εσ(T^4 - Tar^4).
 * With crust: find outer T_s from (T - T_s)/R = h(T_s-Ta)+εσ(T_s^4-Tar^4), q=(T-T_s)/R.
 * Convective / radiative parts are evaluated at T_s (or T if no crust).
 */
inline void frenchFreeSurfaceHeatFlux(Real t_melt, Real h_free, Real emissivity, Real t_ambient, Real t_rad_ambient,
                                      Real stefan_boltzmann, Real enable_crust_flag, Real crust_delta, Real crust_k,
                                      Real &q_total, Real &q_conv, Real &q_rad)
{
    const Real tar2 = t_rad_ambient * t_rad_ambient;
    if (enable_crust_flag < Real(0.5) || crust_delta <= TinyReal)
    {
        const Real t2 = t_melt * t_melt;
        q_conv = h_free * (t_melt - t_ambient);
        q_rad = emissivity * stefan_boltzmann * (t2 * t2 - tar2 * tar2);
        q_total = q_conv + q_rad;
        return;
    }

    const Real k_safe = crust_k > TinyReal ? crust_k : TinyReal;
    const Real r_cond = crust_delta / k_safe;
    Real t_s = t_melt;
    for (int it = 0; it < 12; ++it)
    {
        const Real t2 = t_s * t_s;
        const Real qc = h_free * (t_s - t_ambient);
        const Real qr = emissivity * stefan_boltzmann * (t2 * t2 - tar2 * tar2);
        const Real f = (t_melt - t_s) / r_cond - (qc + qr);
        const Real df = Real(-1) / r_cond - h_free - Real(4) * emissivity * stefan_boltzmann * t_s * t_s * t_s;
        const Real denom = df >= Real(0) ? (df + TinyReal) : (df - TinyReal);
        const Real step = f / denom;
        t_s -= step;
        if (t_s < Real(200))
        {
            t_s = Real(200);
        }
        if (t_s > t_melt)
        {
            t_s = t_melt;
        }
        const Real abs_step = step >= Real(0) ? step : -step;
        if (abs_step < Real(1.0e-3))
        {
            break;
        }
    }
    {
        const Real t2 = t_s * t_s;
        q_conv = h_free * (t_s - t_ambient);
        q_rad = emissivity * stefan_boltzmann * (t2 * t2 - tar2 * tar2);
    }
    q_total = (t_melt - t_s) / r_cond;
}

/** Combined wall flux [W/m^2] using non-exclusive face weights. Shared by kernel and host audit. */
inline Real frenchCombinedBoundaryHeatFluxWm2(Real t, Real w_side, Real w_bottom, Real w_free, Real h_side,
                                             Real h_bottom, Real t_cool, Real h_free, Real emissivity, Real t_ambient,
                                             Real t_rad_ambient, Real stefan_boltzmann, Real enable_crust,
                                             Real crust_delta, Real crust_k, Real &q_side, Real &q_bottom,
                                             Real &q_free_total, Real &q_free_conv, Real &q_free_rad)
{
    q_side = h_side * (t - t_cool);
    q_bottom = h_bottom * (t - t_cool);
    q_free_total = Real(0);
    q_free_conv = Real(0);
    q_free_rad = Real(0);
    if (w_free > TinyReal)
    {
        frenchFreeSurfaceHeatFlux(t, h_free, emissivity, t_ambient, t_rad_ambient, stefan_boltzmann, enable_crust,
                                  crust_delta, crust_k, q_free_total, q_free_conv, q_free_rad);
    }
    return w_side * q_side + w_bottom * q_bottom + w_free * q_free_total;
}

struct OphelieThermalDiffusionOneWayStepResult : public OphelieJouleHeatOneWayStepResult
{
};

inline void registerOphelieThermalDiffusionAuxFields(BaseParticles &particles, Real thermal_conductivity)
{
    particles.registerStateVariable<Real>(kOphelieThermalLaplaceTField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalConductivityField, thermal_conductivity);
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryMaskField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryFaceField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryWeightSideField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryWeightBottomField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryWeightFreeField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalWallGapField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalWallTemperatureField, Real(0));
    const size_t n = particles.TotalRealParticles();
    Real *conductivity = particles.getVariableDataByName<Real>(kOphelieThermalConductivityField);
    for (size_t i = 0; i < n; ++i)
    {
        conductivity[i] = thermal_conductivity;
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalLaplaceTField);
    syncVariableToDevice<Real>(particles, kOphelieThermalConductivityField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryFaceField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToDevice<Real>(particles, kOphelieThermalWallGapField);
    syncVariableToDevice<Real>(particles, kOphelieThermalWallTemperatureField);
}

/** Host: mark boundary shell particles for Dirichlet regression T = T_wall. */
inline size_t setupOphelieThermalDirichletBoundaryMask(BaseParticles &particles,
                                                       OphelieThermalDiffusionOneWayOptions &options,
                                                       const OpheliePhiBoundaryGeometryContext &geom, Real dp)
{
    if (!options.enable_cold_wall_dirichlet)
    {
        return 0;
    }
    options.boundary_shell_thickness = options.boundary_width_factor * dp;
    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = options.boundary_shell_thickness;
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *mask = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryMaskField);
    Real *face = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryFaceField);
    size_t n_boundary = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dist = boundaryDistanceFromContext(pos[i], geom);
        if (dist <= boundary_width)
        {
            mask[i] = 1.0;
            face[i] = 1.0; // undifferentiated wall for regression
            ++n_boundary;
        }
        else
        {
            mask[i] = 0.0;
            face[i] = 0.0;
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryFaceField);
    return n_boundary;
}

/**
 * Host: classify French cylinder shell into side / bottom / free-top.
 * Exclusive `face` remains the closest-face diagnostic.
 * Weights are non-exclusive (doc 40 R03): a corner particle in two shells contributes to both fluxes.
 */
inline void renormalizeOphelieThermalFrenchBoundaryAreasToGeometry(BaseParticles &particles,
                                                                  OphelieThermalDiffusionOneWayOptions &options,
                                                                  const OphelieFrenchReducedCaseParams &french);

inline size_t setupOphelieThermalFrenchNaturalBoundaryFaces(BaseParticles &particles,
                                                            OphelieThermalDiffusionOneWayOptions &options,
                                                            const OphelieFrenchReducedCaseParams &french)
{
    if (!options.enable_french_natural_bc)
    {
        return 0;
    }
    options.boundary_shell_thickness = options.boundary_width_factor * french.dp;
    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = options.boundary_shell_thickness;
    const Vecd &center = french.glass_center;
    const Real z_lo = center[2] - french.glass_half_height;
    const Real z_hi = center[2] + french.glass_half_height;
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *mask = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryMaskField);
    Real *face = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryFaceField);
    Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    size_t n_boundary = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real dist_side = french.glass_radius - r;
        const Real dist_bottom = pos[i][2] - z_lo;
        const Real dist_top = z_hi - pos[i][2];
        w_side[i] = dist_side <= boundary_width ? Real(1) : Real(0);
        w_bottom[i] = dist_bottom <= boundary_width ? Real(1) : Real(0);
        w_free[i] = dist_top <= boundary_width ? Real(1) : Real(0);
        const Real wsum = w_side[i] + w_bottom[i] + w_free[i];
        if (wsum <= TinyReal)
        {
            mask[i] = 0.0;
            face[i] = 0.0;
            continue;
        }
        mask[i] = 1.0;
        if (dist_side <= dist_bottom && dist_side <= dist_top)
        {
            face[i] = 1.0;
        }
        else if (dist_bottom <= dist_top)
        {
            face[i] = 2.0;
        }
        else
        {
            face[i] = 3.0;
        }
        ++n_boundary;
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryFaceField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    if (options.renormalize_boundary_area_to_geometry)
    {
        renormalizeOphelieThermalFrenchBoundaryAreasToGeometry(particles, options, french);
    }
    return n_boundary;
}

/**
 * Scale non-exclusive face weights so discrete A_f = Σ w_f V_i/δ equals geometric A_f
 * of the tagging cylinder (doc 40 R03). Re-run after every retag. Logs the three scales.
 */
inline void renormalizeOphelieThermalFrenchBoundaryAreasToGeometry(BaseParticles &particles,
                                                                  OphelieThermalDiffusionOneWayOptions &options,
                                                                  const OphelieFrenchReducedCaseParams &french)
{
    options.area_scale_side = Real(1);
    options.area_scale_bottom = Real(1);
    options.area_scale_free = Real(1);
    if (!options.enable_french_natural_bc || options.boundary_shell_thickness <= TinyReal)
    {
        return;
    }
    const Real height = Real(2) * french.glass_half_height;
    const Real a_geom_side = Real(2) * Real(M_PI) * french.glass_radius * height;
    const Real a_geom_bottom = Real(M_PI) * french.glass_radius * french.glass_radius;
    const Real a_geom_free = a_geom_bottom;
    const Real inv_delta = Real(1) / (options.boundary_shell_thickness + TinyReal);

    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();

    double a_side = 0.0;
    double a_bottom = 0.0;
    double a_free = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double area = static_cast<double>(vol[i]) * static_cast<double>(inv_delta);
        a_side += static_cast<double>(w_side[i]) * area;
        a_bottom += static_cast<double>(w_bottom[i]) * area;
        a_free += static_cast<double>(w_free[i]) * area;
    }

    const double s_side = a_side > 1.0e-30 ? static_cast<double>(a_geom_side) / a_side : 1.0;
    const double s_bottom = a_bottom > 1.0e-30 ? static_cast<double>(a_geom_bottom) / a_bottom : 1.0;
    const double s_free = a_free > 1.0e-30 ? static_cast<double>(a_geom_free) / a_free : 1.0;
    options.area_scale_side = static_cast<Real>(s_side);
    options.area_scale_bottom = static_cast<Real>(s_bottom);
    options.area_scale_free = static_cast<Real>(s_free);

    for (size_t i = 0; i < n; ++i)
    {
        w_side[i] = static_cast<Real>(static_cast<double>(w_side[i]) * s_side);
        w_bottom[i] = static_cast<Real>(static_cast<double>(w_bottom[i]) * s_bottom);
        w_free[i] = static_cast<Real>(static_cast<double>(w_free[i]) * s_free);
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightFreeField);

    std::cout << "[ophelie][boundary-area-normalize] tagging_cylinder R=" << french.glass_radius
              << " H=" << height << " shell=" << options.boundary_shell_thickness
              << " A_disc_raw side/bottom/free=" << a_side << "/" << a_bottom << "/" << a_free
              << " A_geom=" << a_geom_side << "/" << a_geom_bottom << "/" << a_geom_free
              << " scale side/bottom/free=" << s_side << "/" << s_bottom << "/" << s_free
              << " (explicit SPH→geometry measure; h and skull unchanged)" << std::endl;
}

class ApplyOphelieJouleHeatDiffusionCombinedStepCK : public LocalDynamics
{
  public:
    /**
     * Laplace field must already include conductivity (OpheliePairwiseLaplaceCK harmonic-k).
     * Therefore the diffusion factor is only 1/(ρ cp), not k/(ρ cp). Passing `k` is kept for
     * call-site compatibility and is intentionally unused.
     */
    ApplyOphelieJouleHeatDiffusionCombinedStepCK(SPHBody &sph_body, const std::string &q_field,
                                                 const std::string &delta_t_field,
                                                 const std::string &temperature_field,
                                                 const std::string &laplace_t_field, Real t_initial, Real dt,
                                                 Real rho, Real cp, Real /*k*/, bool apply_diffusion)
        : LocalDynamics(sph_body), t_initial_(t_initial), joule_factor_(dt / (rho * cp + TinyReal)),
          diffusion_factor_(apply_diffusion ? -dt / (rho * cp + TinyReal) : Real(0)),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_laplace_t_(particles_->template getVariableByName<Real>(laplace_t_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), joule_factor_(encloser.joule_factor_),
              diffusion_factor_(encloser.diffusion_factor_), q_(encloser.dv_q_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              laplace_t_(encloser.dv_laplace_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            delta_t_[index_i] += joule_factor_ * q_[index_i] + diffusion_factor_ * laplace_t_[index_i];
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real joule_factor_;
        Real diffusion_factor_;
        Real *q_;
        Real *delta_t_;
        Real *temperature_;
        Real *laplace_t_;
    };

  protected:
    Real t_initial_;
    Real joule_factor_;
    Real diffusion_factor_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_laplace_t_;
};

class ApplyOphelieThermalDirichletWallCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalDirichletWallCK(SPHBody &sph_body, const std::string &delta_t_field,
                                       const std::string &temperature_field, Real t_wall)
        : LocalDynamics(sph_body), t_wall_(t_wall),
          dv_boundary_mask_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryMaskField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_wall_(encloser.t_wall_), boundary_mask_(encloser.dv_boundary_mask_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (boundary_mask_[index_i] < Real(0.5))
            {
                return;
            }
            delta_t_[index_i] = Real(0);
            temperature_[index_i] = t_wall_;
        }

      protected:
        Real t_wall_;
        Real *boundary_mask_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_wall_;
    DiscreteVariable<Real> *dv_boundary_mask_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

/**
 * Stage 4.0: Robin side/bottom + free-surface convection/radiation.
 * Surface flux q [W/m^2] is distributed over the boundary shell volume via q / δ.
 */
class ApplyOphelieThermalFrenchNaturalBcCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalFrenchNaturalBcCK(SPHBody &sph_body, const std::string &delta_t_field,
                                         const std::string &temperature_field, Real t_initial, Real rho, Real cp,
                                         Real shell_thickness, const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamics(sph_body), t_initial_(t_initial),
          inv_rho_cp_delta_(Real(1) / ((rho * cp + TinyReal) * (shell_thickness + TinyReal))), h_side_(bc.h_side),
          h_bottom_(bc.h_bottom), h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool),
          t_ambient_(bc.t_ambient), t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          enable_free_crust_(bc.enable_free_crust ? Real(1) : Real(0)), free_crust_delta_(bc.free_crust_delta_m),
          free_crust_k_(bc.free_crust_k),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_rho_cp_delta_(encloser.inv_rho_cp_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_), enable_free_crust_(encloser.enable_free_crust_),
              free_crust_delta_(encloser.free_crust_delta_), free_crust_k_(encloser.free_crust_k_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            // Continuous area weights (doc 45): inactive face is zero weight, not <0.5.
            if (ws + wb + wf <= TinyReal)
            {
                return;
            }
            const Real t = t_initial_ + delta_t_[index_i];
            Real q_side = 0, q_bottom = 0, q_free = 0, q_conv = 0, q_rad = 0;
            const Real q_out = frenchCombinedBoundaryHeatFluxWm2(
                t, ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_, t_ambient_, t_rad_ambient_,
                stefan_boltzmann_, enable_free_crust_, free_crust_delta_, free_crust_k_, q_side, q_bottom, q_free,
                q_conv, q_rad);
            (void)q_conv;
            (void)q_rad;
            delta_t_[index_i] -= dt * inv_rho_cp_delta_ * q_out;
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real inv_rho_cp_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real enable_free_crust_;
        Real free_crust_delta_;
        Real free_crust_k_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    Real inv_rho_cp_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    Real enable_free_crust_;
    Real free_crust_delta_;
    Real free_crust_k_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

/**
 * Doc207 B: series-gap Robin flux [W/m^2].
 * Continuum gap from particle representative Ti at distance d to the physical face:
 *   q = k(Ti - Tw)/d = h(Tw - T∞) + ε σ (Tw^4 - Trad^4)
 * Convection-only closed form; with radiation, bounded bisection on Tw.
 * d→0 recovers q = h(Ti-T∞)+rad(Ti) (legacy particle-as-wall limit).
 * Does NOT model interior SPH Laplace resistance — only the unresolved particle→face gap.
 */
inline Real ophelieSeriesGapFaceFluxWm2(Real Ti, Real d, Real k_cond, Real h, Real t_inf, Real emissivity,
                                        Real t_rad_ambient, Real stefan_boltzmann, Real &Tw_out)
{
    if (h <= TinyReal && emissivity <= TinyReal)
    {
        Tw_out = Ti;
        return Real(0);
    }
    if (d <= TinyReal || k_cond <= TinyReal)
    {
        Tw_out = Ti;
        Real q = h * (Ti - t_inf);
        if (emissivity > TinyReal)
        {
            const Real t2 = Ti * Ti;
            const Real tr2 = t_rad_ambient * t_rad_ambient;
            q += emissivity * stefan_boltzmann * (t2 * t2 - tr2 * tr2);
        }
        return q;
    }
    if (emissivity <= TinyReal)
    {
        const Real Bi = h * d / (k_cond + TinyReal);
        Tw_out = (Ti + Bi * t_inf) / (Real(1) + Bi);
        return h * (Tw_out - t_inf);
    }
    // Nonlinear (no lambda — SYCL device): f(Tw)=k(Ti-Tw)/d - h(Tw-t_inf) - εσ(Tw^4-Trad^4)
    Real lo = Ti < t_inf ? Ti : t_inf;
    lo = lo < t_rad_ambient ? lo : t_rad_ambient;
    Real hi = Ti > t_inf ? Ti : t_inf;
    hi = hi > t_rad_ambient ? hi : t_rad_ambient;
    lo = lo > Real(1) ? lo - Real(50) : Real(1);
    hi = hi + Real(50);
    const Real tr2 = t_rad_ambient * t_rad_ambient;
    const Real rad_c = emissivity * stefan_boltzmann;
    const Real inv_d = Real(1) / d;
    Real t2lo = lo * lo;
    Real flo = k_cond * (Ti - lo) * inv_d - h * (lo - t_inf) - rad_c * (t2lo * t2lo - tr2 * tr2);
    Real t2hi = hi * hi;
    Real fhi = k_cond * (Ti - hi) * inv_d - h * (hi - t_inf) - rad_c * (t2hi * t2hi - tr2 * tr2);
    if (flo * fhi > Real(0))
    {
        const Real Bi = h * d / (k_cond + TinyReal);
        Tw_out = (Ti + Bi * t_inf) / (Real(1) + Bi);
    }
    else
    {
        for (int it = 0; it < 40; ++it)
        {
            const Real mid = Real(0.5) * (lo + hi);
            const Real t2m = mid * mid;
            const Real fm = k_cond * (Ti - mid) * inv_d - h * (mid - t_inf) - rad_c * (t2m * t2m - tr2 * tr2);
            if (flo * fm <= Real(0))
            {
                hi = mid;
                fhi = fm;
            }
            else
            {
                lo = mid;
                flo = fm;
            }
        }
        Tw_out = Real(0.5) * (lo + hi);
    }
    const Real t2 = Tw_out * Tw_out;
    return h * (Tw_out - t_inf) + rad_c * (t2 * t2 - tr2 * tr2);
}

/**
 * Combined face flux for wall-temperature model (same weight assembly as frenchCombinedBoundaryHeatFluxWm2).
 * Side/bottom: series-gap convection to t_cool. Free: series-gap conv+rad (no crust in Doc207 B candidate).
 */
inline Real ophelieWallTemperatureCombinedFluxWm2(Real Ti, Real d, Real k_cond, Real w_side, Real w_bottom, Real w_free,
                                                 Real h_side, Real h_bottom, Real t_cool, Real h_free, Real emissivity,
                                                 Real t_ambient, Real t_rad_ambient, Real stefan_boltzmann,
                                                 Real &q_side, Real &q_bottom, Real &q_free, Real &Tw_rep)
{
    q_side = q_bottom = q_free = Real(0);
    Tw_rep = Ti;
    Real Tw_s = Ti, Tw_b = Ti, Tw_f = Ti;
    if (w_side > TinyReal)
        q_side = ophelieSeriesGapFaceFluxWm2(Ti, d, k_cond, h_side, t_cool, Real(0), t_rad_ambient, stefan_boltzmann,
                                             Tw_s);
    if (w_bottom > TinyReal)
        q_bottom = ophelieSeriesGapFaceFluxWm2(Ti, d, k_cond, h_bottom, t_cool, Real(0), t_rad_ambient, stefan_boltzmann,
                                               Tw_b);
    if (w_free > TinyReal)
        q_free = ophelieSeriesGapFaceFluxWm2(Ti, d, k_cond, h_free, t_ambient, emissivity, t_rad_ambient,
                                             stefan_boltzmann, Tw_f);
    const Real wsum = w_side + w_bottom + w_free + TinyReal;
    Tw_rep = (w_side * Tw_s + w_bottom * Tw_b + w_free * Tw_f) / wsum;
    return w_side * q_side + w_bottom * q_bottom + w_free * q_free;
}

/**
 * Doc207 B candidate: Apply Robin using reconstructed Tw from geometric wall gap.
 * Same shell area distribution ΔT -= dt * q_out/(ρ cp δ) as legacy FrenchNaturalBcCK — replaces flux only.
 */
class ApplyOphelieThermalWallTemperatureBcCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalWallTemperatureBcCK(SPHBody &sph_body, const std::string &delta_t_field,
                                           const std::string &temperature_field, Real t_initial, Real rho, Real cp,
                                           Real shell_thickness, const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamics(sph_body), t_initial_(t_initial),
          inv_rho_cp_delta_(Real(1) / ((rho * cp + TinyReal) * (shell_thickness + TinyReal))), h_side_(bc.h_side),
          h_bottom_(bc.h_bottom), h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool),
          t_ambient_(bc.t_ambient), t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_gap_(particles_->template getVariableByName<Real>(kOphelieThermalWallGapField)),
          dv_k_(particles_->template getVariableByName<Real>(kOphelieThermalConductivityField)),
          dv_tw_(particles_->template getVariableByName<Real>(kOphelieThermalWallTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_rho_cp_delta_(encloser.inv_rho_cp_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)),
              gap_(encloser.dv_gap_->DelegatedData(ex_policy)), k_(encloser.dv_k_->DelegatedData(ex_policy)),
              tw_(encloser.dv_tw_->DelegatedData(ex_policy)), delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            if (ws + wb + wf <= TinyReal)
            {
                return;
            }
            const Real Ti = t_initial_ + delta_t_[index_i];
            Real q_side = 0, q_bottom = 0, q_free = 0, Tw = Ti;
            const Real q_out = ophelieWallTemperatureCombinedFluxWm2(
                Ti, gap_[index_i], k_[index_i], ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_,
                t_ambient_, t_rad_ambient_, stefan_boltzmann_, q_side, q_bottom, q_free, Tw);
            (void)q_side;
            (void)q_bottom;
            (void)q_free;
            tw_[index_i] = Tw;
            delta_t_[index_i] -= dt * inv_rho_cp_delta_ * q_out;
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real inv_rho_cp_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *gap_;
        Real *k_;
        Real *tw_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    Real inv_rho_cp_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_gap_;
    DiscreteVariable<Real> *dv_k_;
    DiscreteVariable<Real> *dv_tw_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

/** Ledger twin of ApplyOphelieThermalWallTemperatureBcCK (call before BC update). */
class OphelieWallTemperatureBoundaryLossPowerReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieWallTemperatureBoundaryLossPowerReduceCK(SPHBody &sph_body, Real t_initial, Real shell_thickness,
                                                    const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), t_initial_(t_initial),
          inv_delta_(Real(1) / (shell_thickness + TinyReal)), h_side_(bc.h_side), h_bottom_(bc.h_bottom),
          h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool), t_ambient_(bc.t_ambient),
          t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_gap_(particles_->template getVariableByName<Real>(kOphelieThermalWallGapField)),
          dv_k_(particles_->template getVariableByName<Real>(kOphelieThermalConductivityField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieWallTemperatureBoundaryLossPower";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_delta_(encloser.inv_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)),
              gap_(encloser.dv_gap_->DelegatedData(ex_policy)), k_(encloser.dv_k_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            if (ws + wb + wf <= TinyReal)
                return Real(0);
            const Real Ti = t_initial_ + delta_t_[index_i];
            const Real area = vol_[index_i] * inv_delta_;
            Real q_side = 0, q_bottom = 0, q_free = 0, Tw = Ti;
            const Real q = ophelieWallTemperatureCombinedFluxWm2(
                Ti, gap_[index_i], k_[index_i], ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_,
                t_ambient_, t_rad_ambient_, stefan_boltzmann_, q_side, q_bottom, q_free, Tw);
            (void)Tw;
            return q * area;
        }

      protected:
        Real t_initial_;
        Real inv_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *gap_;
        Real *k_;
        Real *delta_t_;
        Real *vol_;
    };

  protected:
    Real t_initial_;
    Real inv_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_gap_;
    DiscreteVariable<Real> *dv_k_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_vol_;
};

/** Doc207 B + variable-cp: same series-gap flux as WallTemperatureBcCK; ΔT uses per-particle cp. */
class ApplyOphelieThermalWallTemperatureBcParticleCpCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalWallTemperatureBcParticleCpCK(SPHBody &sph_body, const std::string &delta_t_field,
                                                     const std::string &temperature_field, Real t_initial, Real rho,
                                                     Real shell_thickness, const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamics(sph_body), t_initial_(t_initial),
          inv_rho_delta_(Real(1) / ((rho + TinyReal) * (shell_thickness + TinyReal))), h_side_(bc.h_side),
          h_bottom_(bc.h_bottom), h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool),
          t_ambient_(bc.t_ambient), t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_gap_(particles_->template getVariableByName<Real>(kOphelieThermalWallGapField)),
          dv_k_(particles_->template getVariableByName<Real>(kOphelieThermalConductivityField)),
          dv_tw_(particles_->template getVariableByName<Real>(kOphelieThermalWallTemperatureField)),
          dv_cp_(particles_->template getVariableByName<Real>("OphelieHeatCapacity")),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_rho_delta_(encloser.inv_rho_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)),
              gap_(encloser.dv_gap_->DelegatedData(ex_policy)), k_(encloser.dv_k_->DelegatedData(ex_policy)),
              tw_(encloser.dv_tw_->DelegatedData(ex_policy)), cp_(encloser.dv_cp_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            if (ws + wb + wf <= TinyReal)
                return;
            const Real Ti = t_initial_ + delta_t_[index_i];
            Real q_side = 0, q_bottom = 0, q_free = 0, Tw = Ti;
            const Real q_out = ophelieWallTemperatureCombinedFluxWm2(
                Ti, gap_[index_i], k_[index_i], ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_,
                t_ambient_, t_rad_ambient_, stefan_boltzmann_, q_side, q_bottom, q_free, Tw);
            (void)q_side;
            (void)q_bottom;
            (void)q_free;
            tw_[index_i] = Tw;
            delta_t_[index_i] -= dt * inv_rho_delta_ * q_out / (cp_[index_i] + TinyReal);
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real inv_rho_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *gap_;
        Real *k_;
        Real *tw_;
        Real *cp_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    Real inv_rho_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_gap_;
    DiscreteVariable<Real> *dv_k_;
    DiscreteVariable<Real> *dv_tw_;
    DiscreteVariable<Real> *dv_cp_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

/** Host face-power audit for wall-temperature (series-gap); free rad folded into free_rad_w approx via residual. */
inline void hostOphelieThermalWallTemperatureHeatLossPowers(BaseParticles &particles,
                                                            const OphelieThermalDiffusionOneWayOptions &options,
                                                            Real shell_thickness, Real &side_w, Real &bottom_w,
                                                            Real &free_conv_w, Real &free_rad_w, Real &total_w)
{
    side_w = bottom_w = free_conv_w = free_rad_w = total_w = Real(0);
    if (!options.enable_french_natural_bc)
        return;
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToHost<Real>(particles, kOphelieThermalWallGapField);
    syncVariableToHost<Real>(particles, kOphelieThermalConductivityField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    const Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    const Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    const Real *gap = particles.getVariableDataByName<Real>(kOphelieThermalWallGapField);
    const Real *kcond = particles.getVariableDataByName<Real>(kOphelieThermalConductivityField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real inv_delta = Real(1) / (shell_thickness + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        if (w_side[i] + w_bottom[i] + w_free[i] <= TinyReal)
            continue;
        const Real area = vol[i] * inv_delta;
        Real q_side = 0, q_bottom = 0, q_free = 0, Tw = temperature[i];
        (void)ophelieWallTemperatureCombinedFluxWm2(
            temperature[i], gap[i], kcond[i], w_side[i], w_bottom[i], w_free[i], options.h_side, options.h_bottom,
            options.t_cool, options.h_free, options.emissivity, options.t_ambient, options.t_rad_ambient,
            options.stefan_boltzmann, q_side, q_bottom, q_free, Tw);
        // Split free into conv vs rad at reconstructed Tw for audit (same Tw used in combined).
        Real q_conv = options.h_free * (Tw - options.t_ambient);
        const Real t2 = Tw * Tw;
        const Real tr2 = options.t_rad_ambient * options.t_rad_ambient;
        Real q_rad = options.emissivity * options.stefan_boltzmann * (t2 * t2 - tr2 * tr2);
        if (w_free[i] <= TinyReal)
        {
            q_conv = 0;
            q_rad = 0;
        }
        side_w += w_side[i] * q_side * area;
        bottom_w += w_bottom[i] * q_bottom * area;
        free_conv_w += w_free[i] * q_conv * area;
        free_rad_w += w_free[i] * q_rad * area;
        total_w += (w_side[i] * q_side + w_bottom[i] * q_bottom + w_free[i] * q_free) * area;
    }
}

/** Host audit of instantaneous Natural BC heat-loss powers at the current temperature field. */
inline void hostOphelieThermalFrenchNaturalHeatLossPowers(BaseParticles &particles,
                                                          const OphelieThermalDiffusionOneWayOptions &options,
                                                          Real shell_thickness, Real &side_w, Real &bottom_w,
                                                          Real &free_conv_w, Real &free_rad_w, Real &total_w)
{
    side_w = bottom_w = free_conv_w = free_rad_w = total_w = Real(0);
    if (!options.enable_french_natural_bc)
    {
        return;
    }
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    const Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    const Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real inv_delta = Real(1) / (shell_thickness + TinyReal);
    const Real enable_crust = options.enable_free_crust ? Real(1) : Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        // Doc 45: continuous area weights — inactive is zero weight, not the old binary <0.5 tag.
        if (w_side[i] + w_bottom[i] + w_free[i] <= TinyReal)
        {
            continue;
        }
        const Real area = vol[i] * inv_delta;
        Real q_side = 0, q_bottom = 0, q_free = 0, q_conv = 0, q_rad = 0;
        (void)frenchCombinedBoundaryHeatFluxWm2(temperature[i], w_side[i], w_bottom[i], w_free[i], options.h_side,
                                                options.h_bottom, options.t_cool, options.h_free, options.emissivity,
                                                options.t_ambient, options.t_rad_ambient, options.stefan_boltzmann,
                                                enable_crust, options.free_crust_delta_m, options.free_crust_k, q_side,
                                                q_bottom, q_free, q_conv, q_rad);
        side_w += w_side[i] * q_side * area;
        bottom_w += w_bottom[i] * q_bottom * area;
        free_conv_w += w_free[i] * q_conv * area;
        free_rad_w += w_free[i] * q_rad * area;
        total_w += (w_side[i] * q_side + w_bottom[i] * q_bottom + w_free[i] * q_free) * area;
    }
}

/**
 * Doc208 A: single host dispatch for face-power audit.
 * Matches the applied BC path (legacy-shell vs wall-temperature). Read-only: no T advance, no ledger mutate.
 * @return formula tag written into monitors ("wall-temperature" | "legacy-shell").
 */
inline const char *hostOphelieThermalAppliedBcHeatLossPowers(BaseParticles &particles,
                                                             const OphelieThermalDiffusionOneWayOptions &options,
                                                             Real shell_thickness, bool use_wall_temperature,
                                                             Real &side_w, Real &bottom_w, Real &free_conv_w,
                                                             Real &free_rad_w, Real &total_w)
{
    if (use_wall_temperature)
    {
        hostOphelieThermalWallTemperatureHeatLossPowers(particles, options, shell_thickness, side_w, bottom_w,
                                                        free_conv_w, free_rad_w, total_w);
        return "wall-temperature";
    }
    hostOphelieThermalFrenchNaturalHeatLossPowers(particles, options, shell_thickness, side_w, bottom_w, free_conv_w,
                                                  free_rad_w, total_w);
    return "legacy-shell";
}

/**
 * Device reduce of French natural BC heat-loss power [W] at the current thermal state.
 * Uses T = T0 + ΔT (doc 45 F01) and the same flux weights as ApplyOphelieThermalFrenchNaturalBcCK.
 * Call immediately before the BC update so the ledger integrates the flux the operator applies (F02).
 */
class OphelieFrenchBoundaryLossPowerReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieFrenchBoundaryLossPowerReduceCK(SPHBody &sph_body, Real t_initial, Real shell_thickness,
                                           const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), t_initial_(t_initial),
          inv_delta_(Real(1) / (shell_thickness + TinyReal)), h_side_(bc.h_side), h_bottom_(bc.h_bottom),
          h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool), t_ambient_(bc.t_ambient),
          t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          enable_free_crust_(bc.enable_free_crust ? Real(1) : Real(0)), free_crust_delta_(bc.free_crust_delta_m),
          free_crust_k_(bc.free_crust_k),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieFrenchBoundaryLossPower";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_delta_(encloser.inv_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_), enable_free_crust_(encloser.enable_free_crust_),
              free_crust_delta_(encloser.free_crust_delta_), free_crust_k_(encloser.free_crust_k_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            if (ws + wb + wf <= TinyReal)
            {
                return Real(0);
            }
            const Real t = t_initial_ + delta_t_[index_i];
            const Real area = vol_[index_i] * inv_delta_;
            Real q_side = 0, q_bottom = 0, q_free = 0, q_conv = 0, q_rad = 0;
            (void)frenchCombinedBoundaryHeatFluxWm2(t, ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_,
                                                    t_ambient_, t_rad_ambient_, stefan_boltzmann_, enable_free_crust_,
                                                    free_crust_delta_, free_crust_k_, q_side, q_bottom, q_free, q_conv,
                                                    q_rad);
            return (ws * q_side + wb * q_bottom + wf * q_free) * area;
        }

      protected:
        Real t_initial_;
        Real inv_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real enable_free_crust_;
        Real free_crust_delta_;
        Real free_crust_k_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *delta_t_;
        Real *vol_;
    };

  protected:
    Real t_initial_;
    Real inv_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    Real enable_free_crust_;
    Real free_crust_delta_;
    Real free_crust_k_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_vol_;
};

/** Discrete vs geometric boundary-area audit (doc 23 P0-06). Area proxy: A_i = V_i / shell_thickness. */
struct OphelieFrenchBoundaryAreaAudit
{
    Real shell_thickness = 0.0;
    Real geometric_area_side = 0.0;
    Real geometric_area_bottom = 0.0;
    Real geometric_area_free = 0.0;
    Real discrete_area_side = 0.0;
    Real discrete_area_bottom = 0.0;
    Real discrete_area_free = 0.0;
    size_t n_side = 0;
    size_t n_bottom = 0;
    size_t n_free = 0;
    size_t n_corner = 0;
    size_t n_multi_face = 0;
    Real mean_t_side = 0.0;
    Real mean_t_bottom = 0.0;
    Real mean_t_free = 0.0;
    Real side_w = 0.0;
    Real bottom_w = 0.0;
    Real free_conv_w = 0.0;
    Real free_rad_w = 0.0;
    Real total_w = 0.0;
    Real area_side_ratio = 0.0;
    Real area_bottom_ratio = 0.0;
    Real area_free_ratio = 0.0;
};

inline OphelieFrenchBoundaryAreaAudit hostOphelieThermalFrenchNaturalBoundaryAreaAudit(
    BaseParticles &particles, const OphelieThermalDiffusionOneWayOptions &options,
    const OphelieFrenchReducedCaseParams &french, Real shell_thickness)
{
    OphelieFrenchBoundaryAreaAudit a;
    a.shell_thickness = shell_thickness;
    const Real height = Real(2) * french.glass_half_height;
    a.geometric_area_side = Real(2) * Real(M_PI) * french.glass_radius * height;
    a.geometric_area_bottom = Real(M_PI) * french.glass_radius * french.glass_radius;
    a.geometric_area_free = a.geometric_area_bottom;

    if (!options.enable_french_natural_bc || shell_thickness <= TinyReal)
    {
        return a;
    }

    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    const Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    const Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real inv_delta = Real(1) / (shell_thickness + TinyReal);
    const Real enable_crust = options.enable_free_crust ? Real(1) : Real(0);

    Real t_side_sum = 0, t_bottom_sum = 0, t_free_sum = 0;
    Real a_side_sum = 0, a_bottom_sum = 0, a_free_sum = 0;

    for (size_t i = 0; i < n; ++i)
    {
        const Real ws = w_side[i];
        const Real wb = w_bottom[i];
        const Real wf = w_free[i];
        if (ws + wb + wf < Real(0.5))
        {
            continue;
        }
        const int nfaces = (ws > Real(0.5) ? 1 : 0) + (wb > Real(0.5) ? 1 : 0) + (wf > Real(0.5) ? 1 : 0);
        if (nfaces >= 2)
        {
            ++a.n_multi_face;
        }
        if (ws > Real(0.5) && wb > Real(0.5))
        {
            ++a.n_corner;
        }
        const Real area = vol[i] * inv_delta;
        const Real t = temperature[i];
        Real q_side = 0, q_bottom = 0, q_free = 0, q_conv = 0, q_rad = 0;
        (void)frenchCombinedBoundaryHeatFluxWm2(t, ws, wb, wf, options.h_side, options.h_bottom, options.t_cool,
                                                options.h_free, options.emissivity, options.t_ambient,
                                                options.t_rad_ambient, options.stefan_boltzmann, enable_crust,
                                                options.free_crust_delta_m, options.free_crust_k, q_side, q_bottom,
                                                q_free, q_conv, q_rad);
        if (ws > Real(0.5))
        {
            ++a.n_side;
            a.discrete_area_side += ws * area;
            a_side_sum += ws * area;
            t_side_sum += t * ws * area;
            a.side_w += ws * q_side * area;
        }
        if (wb > Real(0.5))
        {
            ++a.n_bottom;
            a.discrete_area_bottom += wb * area;
            a_bottom_sum += wb * area;
            t_bottom_sum += t * wb * area;
            a.bottom_w += wb * q_bottom * area;
        }
        if (wf > Real(0.5))
        {
            ++a.n_free;
            a.discrete_area_free += wf * area;
            a_free_sum += wf * area;
            t_free_sum += t * wf * area;
            a.free_conv_w += wf * q_conv * area;
            a.free_rad_w += wf * q_rad * area;
        }
        a.total_w += (ws * q_side + wb * q_bottom + wf * q_free) * area;
    }

    a.mean_t_side = a_side_sum > TinyReal ? t_side_sum / a_side_sum : Real(0);
    a.mean_t_bottom = a_bottom_sum > TinyReal ? t_bottom_sum / a_bottom_sum : Real(0);
    a.mean_t_free = a_free_sum > TinyReal ? t_free_sum / a_free_sum : Real(0);
    a.area_side_ratio = a.discrete_area_side / (a.geometric_area_side + TinyReal);
    a.area_bottom_ratio = a.discrete_area_bottom / (a.geometric_area_bottom + TinyReal);
    a.area_free_ratio = a.discrete_area_free / (a.geometric_area_free + TinyReal);
    return a;
}

inline void logOphelieFrenchBoundaryAreaAudit(const OphelieFrenchBoundaryAreaAudit &a, const char *tag)
{
    std::cout << "[ophelie][boundary-area-audit] " << tag << " shell_thickness=" << a.shell_thickness << " m\n"
              << "  geometric A side/bottom/free=" << a.geometric_area_side << "/" << a.geometric_area_bottom << "/"
              << a.geometric_area_free << " m^2\n"
              << "  discrete  A side/bottom/free=" << a.discrete_area_side << "/" << a.discrete_area_bottom << "/"
              << a.discrete_area_free << " m^2\n"
              << "  A_disc/A_geom side/bottom/free=" << a.area_side_ratio << "/" << a.area_bottom_ratio << "/"
              << a.area_free_ratio << "\n"
              << "  n_particles side/bottom/free=" << a.n_side << "/" << a.n_bottom << "/" << a.n_free
              << " n_corner=" << a.n_corner << " n_multi_face=" << a.n_multi_face << " tagging=multi_face\n"
              << "  mean_T side/bottom/free=" << a.mean_t_side << "/" << a.mean_t_bottom << "/" << a.mean_t_free
              << " K\n"
              << "  loss side/bottom/free_conv+rad=" << a.side_w << "/" << a.bottom_w << "/"
              << (a.free_conv_w + a.free_rad_w) << " W  total=" << a.total_w << " W" << std::endl;
}

/**
 * Host Stage-1 identities (doc 40 / 42). No GPU.
 * 1) constant-cp Q update: ΔH = Σ Q V dt
 * 2) Robin cooling: ΔH = −P_loss dt with the same q·(V/δ) weights as the kernel
 * 3) synthetic cylinder: multi-face recovers exclusive-tag corner loss
 *
 * Accumulators are double. Production Real is float: at T0=1473 K, a French-scale
 * explicit ΔT can be smaller than ulp(T), so the same algebra in Real is not a
 * valid GPU identity. Host gate is 1e-8 relative (round-trip through T, not 1e-12).
 */
inline bool verifyFrenchThermalStage1HostAlgebra(std::string &message)
{
    const double rho = 2800.0;
    const double cp = 1200.0;
    const double dt = 1.0e-4;
    const double q = 1.0e6;
    const double vol = 1.0e-6;
    const double t0 = 1473.0;
    const double h0 = rho * cp * t0 * vol;
    const double t1 = t0 + q * dt / (rho * cp);
    const double h1 = rho * cp * t1 * vol;
    const double e_q = q * vol * dt;
    const double q_rel = std::abs(h1 - h0 - e_q) / (std::abs(e_q) + 1.0e-30);
    const float t0f = 1473.0f;
    const float dtf = static_cast<float>(q * dt / (rho * cp));
    const int float_ulp_swallows_dt = (t0f + dtf == t0f) ? 1 : 0;
    if (q_rel > 1.0e-8)
    {
        message = "q_only_identity_failed";
        return false;
    }

    const double delta = 0.0075;
    const double h_side = 300.0;
    const double t_cool = 300.0;
    const double q_wm2 = h_side * (t0 - t_cool);
    const double p_loss = q_wm2 * (vol / delta);
    const double t_cool_end = t0 - dt * q_wm2 / (rho * cp * delta);
    const double h_end = rho * cp * t_cool_end * vol;
    const double cool_rel = std::abs((h_end - h0) + p_loss * dt) / (std::abs(p_loss * dt) + 1.0e-30);
    if (cool_rel > 1.0e-8)
    {
        message = "cooling_identity_failed";
        return false;
    }

    // Doc 45 F01/F05: ΔT-authoritative diffusion with nonzero Laplace must retain sub-ULP increments.
    {
        const float t0f = 1473.0f;
        const float rho_f = 2800.0f;
        const float cp_f = 1150.0f;
        const float dt_f = 3.8e-4f;
        const float laplace_f = -1.0e5f; // W/m·K units absorbed into pairwise Laplace·k
        const float diff_factor = -1.0f / (rho_f * cp_f);
        const float d_step = dt_f * diff_factor * laplace_f; // positive heating-like if laplace negative
        float t_bug = t0f;
        float delta_fix = 0.0f;
        for (int s = 0; s < 10000; ++s)
        {
            t_bug = t_bug + d_step; // buggy absolute path
            delta_fix = delta_fix + d_step;
        }
        const float t_fix = t0f + delta_fix;
        if (!(t_bug == t0f && std::abs(delta_fix) > 1.0e-4f && std::abs(t_fix - t0f - delta_fix) < 1.0e-3f))
        {
            message = "delta_t_diffusion_ulp_contract_failed";
            return false;
        }
        // Inactive clamp must preserve delta (doc 45 F01).
        float delta_clamp = delta_fix;
        float t_clamp = t0f + delta_clamp;
        const float t_min = 300.0f;
        const float t_max = 2000.0f;
        if (t_clamp >= t_min && t_clamp <= t_max)
        {
            // no change to delta
        }
        else
        {
            t_clamp = std::min(t_max, std::max(t_min, t_clamp));
            delta_clamp = t_clamp - t0f;
        }
        if (std::abs(delta_clamp - delta_fix) > 1.0e-6f)
        {
            message = "inactive_clamp_mutated_delta";
            return false;
        }
    }

    // Doc 45 F05: diffusion rate scales with conductivity already inside Laplace.
    // With L = k * L0, ΔT ∝ -dt/(ρ cp) * L ⇒ ratio of increments = k2/k1.
    {
        const double rho = 2800.0;
        const double cp = 1150.0;
        const double dt = 1.0e-3;
        const double L0 = 1.0e4; // Laplace without k
        const double d1 = -dt / (rho * cp) * (1.0 * L0);
        const double d4 = -dt / (rho * cp) * (4.0 * L0);
        if (std::abs(d4 / d1 - 4.0) > 1.0e-12)
        {
            message = "diffusion_k_scaling_failed";
            return false;
        }
    }

    // Doc 47 §6 / Doc 54: production stirring order Q → Laplace → inactive clamp (ΔT-authoritative).
    {
        const float t0f = 1473.0f;
        const float rho_f = 2800.0f;
        const float cp_f = 1150.0f;
        const float dt_f = 3.8e-4f;
        const float q_f = 5.0e4f;
        const float laplace_f = -1.0e5f;
        const float inv_rho_cp = 1.0f / (rho_f * cp_f);
        const float diff_factor = -1.0f / (rho_f * cp_f);
        const float t_min = 300.0f;
        const float t_max = 2000.0f;
        float delta = 0.0f;
        for (int s = 0; s < 5000; ++s)
        {
            delta += q_f * dt_f * inv_rho_cp;           // ApplyGridJouleHeatToTemperatureCK
            delta += dt_f * diff_factor * laplace_f;     // ApplyLaplaceDiffusionDtCK
            float t = t0f + delta;                      // inactive ClampTemperatureCK
            if (t < t_min || t > t_max)
            {
                t = std::min(t_max, std::max(t_min, t));
                delta = t - t0f;
            }
        }
        const float t_end = t0f + delta;
        const float expected_q = q_f * dt_f * inv_rho_cp * 5000.0f;
        const float expected_d = dt_f * diff_factor * laplace_f * 5000.0f;
        if (!(std::abs(delta) > 1.0e-3f && std::abs(t_end - t0f - delta) < 1.0e-3f &&
              std::abs(delta - (expected_q + expected_d)) < 1.0e-2f * std::abs(expected_q + expected_d)))
        {
            message = "production_q_diffusion_clamp_sequence_failed";
            return false;
        }
        // Active clamp must rewrite ΔT when T would exceed band.
        float delta_hi = t_max - t0f + 10.0f;
        float t_hi = t0f + delta_hi;
        if (t_hi > t_max)
        {
            t_hi = t_max;
            delta_hi = t_hi - t0f;
        }
        if (std::abs(delta_hi - (t_max - t0f)) > 1.0e-3f)
        {
            message = "active_clamp_delta_rewrite_failed";
            return false;
        }
    }

    const Real R = Real(0.25);
    const Real H = Real(0.21);
    const Real dp = Real(0.02);
    const Real shell = Real(1.5) * dp;
    const Real geom_side = Real(2) * Real(M_PI) * R * H;
    const Real geom_bottom = Real(M_PI) * R * R;
    Real a_side_ex = 0, a_bot_ex = 0, a_side_mf = 0, a_bot_mf = 0;
    size_t n_corner = 0;
    const int nr = static_cast<int>(std::ceil(R / dp));
    const int nz = static_cast<int>(std::ceil(H / dp));
    for (int iz = 0; iz < nz; ++iz)
    {
        const Real z = (static_cast<Real>(iz) + Real(0.5)) * dp;
        if (z > H)
        {
            continue;
        }
        for (int ir = 0; ir < nr; ++ir)
        {
            const Real r = (static_cast<Real>(ir) + Real(0.5)) * dp;
            if (r > R)
            {
                continue;
            }
            const int nth = std::max(1, static_cast<int>(std::round(Real(2) * Real(M_PI) * r / dp)));
            const Real dtheta = Real(2) * Real(M_PI) / static_cast<Real>(nth);
            const Real vol_cell = r * dtheta * dp * dp * static_cast<Real>(nth);
            const Real dist_side = R - r;
            const Real dist_bot = z;
            const Real dist_top = H - z;
            const Real area = vol_cell / (shell + TinyReal);
            const bool in_side = dist_side <= shell;
            const bool in_bot = dist_bot <= shell;
            const bool in_top = dist_top <= shell;
            if (!in_side && !in_bot && !in_top)
            {
                continue;
            }
            if (in_side && in_bot)
            {
                ++n_corner;
            }
            if (dist_side <= dist_bot && dist_side <= dist_top)
            {
                a_side_ex += area;
            }
            else if (dist_bot <= dist_top)
            {
                a_bot_ex += area;
            }
            if (in_side)
            {
                a_side_mf += area;
            }
            if (in_bot)
            {
                a_bot_mf += area;
            }
            (void)in_top;
        }
    }
    const Real r_side_ex = a_side_ex / (geom_side + TinyReal);
    const Real r_bot_ex = a_bot_ex / (geom_bottom + TinyReal);
    const Real r_side_mf = a_side_mf / (geom_side + TinyReal);
    const Real r_bot_mf = a_bot_mf / (geom_bottom + TinyReal);
    if (n_corner == 0)
    {
        message = "multi_face_no_corners";
        return false;
    }
    if (!(r_side_mf + Real(1.0e-6) >= r_side_ex && r_bot_mf + Real(1.0e-6) >= r_bot_ex))
    {
        message = "multi_face_did_not_recover_corner_area";
        return false;
    }

    // Nonuniform-T flux identity (doc 40 R03): same multi-face weights + normalize as production,
    // T = T0 + az*z + ar*r. Discrete Σ w h (T−Tcool) V/δ must match the geometric face integrals.
    {
        const double R_d = static_cast<double>(R);
        const double H_d = static_cast<double>(H);
        const double shell_d = static_cast<double>(shell);
        const double t0_nu = 1473.0;
        const double t_cool_nu = 300.0;
        const double h_side_nu = 300.0;
        const double h_bot_nu = 35.0;
        const double az = 80.0; // K/m
        const double ar = 40.0; // K/m
        const double a_geom_side = 2.0 * M_PI * R_d * H_d;
        const double a_geom_bot = M_PI * R_d * R_d;
        // Analytic: side at r=R, ∫(T0+az z+ar R−Tcool) R dθ dz
        const double mean_dT_side = (t0_nu - t_cool_nu) + 0.5 * az * H_d + ar * R_d;
        const double p_side_geom = h_side_nu * mean_dT_side * a_geom_side;
        // Bottom at z=0: ∫(T0+ar r−Tcool) 2π r dr → mean ΔT = (T0−Tcool) + (2/3) ar R
        const double mean_dT_bot = (t0_nu - t_cool_nu) + (2.0 / 3.0) * ar * R_d;
        const double p_bot_geom = h_bot_nu * mean_dT_bot * a_geom_bot;

        double a_side_raw = 0.0;
        double a_bot_raw = 0.0;
        double p_side_raw = 0.0;
        double p_bot_raw = 0.0;
        for (int iz = 0; iz < nz; ++iz)
        {
            const double z = (static_cast<double>(iz) + 0.5) * static_cast<double>(dp);
            if (z > H_d)
            {
                continue;
            }
            for (int ir = 0; ir < nr; ++ir)
            {
                const double r = (static_cast<double>(ir) + 0.5) * static_cast<double>(dp);
                if (r > R_d)
                {
                    continue;
                }
                const int nth = std::max(1, static_cast<int>(std::round(2.0 * M_PI * r / static_cast<double>(dp))));
                const double dtheta = 2.0 * M_PI / static_cast<double>(nth);
                const double vol_cell = r * dtheta * static_cast<double>(dp) * static_cast<double>(dp) *
                                       static_cast<double>(nth);
                const double dist_side = R_d - r;
                const double dist_bot = z;
                const double dist_top = H_d - z;
                const bool in_side = dist_side <= shell_d;
                const bool in_bot = dist_bot <= shell_d;
                const bool in_top = dist_top <= shell_d;
                if (!in_side && !in_bot && !in_top)
                {
                    continue;
                }
                const double area = vol_cell / (shell_d + 1.0e-30);
                const double t = t0_nu + az * z + ar * r;
                const double w_s = in_side ? 1.0 : 0.0;
                const double w_b = in_bot ? 1.0 : 0.0;
                if (w_s > 0.0)
                {
                    a_side_raw += area;
                    p_side_raw += w_s * h_side_nu * (t - t_cool_nu) * area;
                }
                if (w_b > 0.0)
                {
                    a_bot_raw += area;
                    p_bot_raw += w_b * h_bot_nu * (t - t_cool_nu) * area;
                }
            }
        }
        const double s_side = a_side_raw > 1.0e-30 ? a_geom_side / a_side_raw : 1.0;
        const double s_bot = a_bot_raw > 1.0e-30 ? a_geom_bot / a_bot_raw : 1.0;
        const double p_side_norm = p_side_raw * s_side;
        const double p_bot_norm = p_bot_raw * s_bot;
        const double rel_side = std::abs(p_side_norm - p_side_geom) / (std::abs(p_side_geom) + 1.0e-30);
        const double rel_bot = std::abs(p_bot_norm - p_bot_geom) / (std::abs(p_bot_geom) + 1.0e-30);
        // Coarse lattice (dp=0.02): allow a few percent; production uses finer packing + normalize.
        if (rel_side > 0.05 || rel_bot > 0.05)
        {
            message = "nonuniform_flux_normalize_failed";
            return false;
        }
        std::ostringstream oss;
        oss << "stage1_host_ok q_rel=" << q_rel << " cool_rel=" << cool_rel
            << " float_ulp_swallows_dT=" << float_ulp_swallows_dt << " n_corner=" << n_corner
            << " A_side_ex/mf=" << r_side_ex << "/" << r_side_mf << " A_bot_ex/mf=" << r_bot_ex << "/" << r_bot_mf
            << " nonuniform_P_rel side/bot=" << rel_side << "/" << rel_bot
            << " (synthetic cylinder; not a 2% production packing gate)";
        message = oss.str();
        return true;
    }
}

inline bool verifyFrenchBoundaryFluxMatchesUpdateWeights(BaseParticles &particles,
                                                         const OphelieThermalDiffusionOneWayOptions &options,
                                                         Real shell_thickness, std::string &message)
{
    if (!options.enable_french_natural_bc)
    {
        message = "french_bc_disabled";
        return true;
    }
    Real side_w = 0, bottom_w = 0, free_c = 0, free_r = 0, total = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(particles, options, shell_thickness, side_w, bottom_w, free_c, free_r,
                                                  total);
    const OphelieFrenchBoundaryAreaAudit audit =
        hostOphelieThermalFrenchNaturalBoundaryAreaAudit(particles, options,
                                                         /*french unused for flux compare*/ OphelieFrenchReducedCaseParams{},
                                                         shell_thickness);
    // Audit.total_w uses the same combined flux; compare to host heat-loss total.
    const Real rel = std::abs(audit.total_w - total) / (std::abs(total) + TinyReal);
    if (rel > Real(1.0e-5))
    {
        message = "boundary_audit_flux_mismatch";
        return false;
    }
    std::ostringstream oss;
    oss << "boundary_flux_weights_ok rel=" << rel << " P_loss=" << total << " n_corner=" << audit.n_corner
        << " A_side_ratio=" << audit.area_side_ratio << " A_bot_ratio=" << audit.area_bottom_ratio;
    message = oss.str();
    return true;
}

/**
 * Production packing check (doc 40 R03): impose T = T0 + az*(z−zlo) + ar*r, compare
 * discrete side/bottom loss (same weights as the update) to geometric face integrals.
 * Restores the previous temperature field. Expects area-normalize on for a tight gate.
 */
inline bool verifyFrenchNonuniformBoundaryFluxOnParticles(BaseParticles &particles,
                                                          OphelieThermalDiffusionOneWayOptions &options,
                                                          const OphelieFrenchReducedCaseParams &french, Real t0,
                                                          std::string &message)
{
    if (!options.enable_french_natural_bc || options.boundary_shell_thickness <= TinyReal)
    {
        message = "nonuniform_flux_skip_no_bc";
        return true;
    }
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Vecd>(particles, "Position");
    Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> t_backup(n);
    for (size_t i = 0; i < n; ++i)
    {
        t_backup[i] = temperature[i];
    }

    const Real z_lo = french.glass_center[2] - french.glass_half_height;
    const Real H = Real(2) * french.glass_half_height;
    const Real R = french.glass_radius;
    const Real az = Real(80);
    const Real ar = Real(40);
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - french.glass_center[0];
        const Real dy = pos[i][1] - french.glass_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real z = pos[i][2] - z_lo;
        temperature[i] = t0 + az * z + ar * r;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);

    Real side_w = 0, bottom_w = 0, free_c = 0, free_r = 0, total = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(particles, options, options.boundary_shell_thickness, side_w,
                                                  bottom_w, free_c, free_r, total);

    const double R_d = static_cast<double>(R);
    const double H_d = static_cast<double>(H);
    const double t0_d = static_cast<double>(t0);
    const double t_cool = static_cast<double>(options.t_cool);
    const double h_side = static_cast<double>(options.h_side);
    const double h_bot = static_cast<double>(options.h_bottom);
    const double az_d = static_cast<double>(az);
    const double ar_d = static_cast<double>(ar);
    const double a_side = 2.0 * M_PI * R_d * H_d;
    const double a_bot = M_PI * R_d * R_d;
    const double p_side_geom = h_side * ((t0_d - t_cool) + 0.5 * az_d * H_d + ar_d * R_d) * a_side;
    const double p_bot_geom = h_bot * ((t0_d - t_cool) + (2.0 / 3.0) * ar_d * R_d) * a_bot;
    const double rel_side = std::abs(static_cast<double>(side_w) - p_side_geom) / (std::abs(p_side_geom) + 1.0e-30);
    const double rel_bot = std::abs(static_cast<double>(bottom_w) - p_bot_geom) / (std::abs(p_bot_geom) + 1.0e-30);

    for (size_t i = 0; i < n; ++i)
    {
        temperature[i] = t_backup[i];
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);

    // With normalize: expect ~few % on production packing. Without: area bias dominates.
    const double tol = options.renormalize_boundary_area_to_geometry ? 0.05 : 0.25;
    if (rel_side > tol || rel_bot > tol)
    {
        std::ostringstream oss;
        oss << "nonuniform_flux_particles_failed rel_side=" << rel_side << " rel_bot=" << rel_bot
            << " tol=" << tol << " normalize=" << (options.renormalize_boundary_area_to_geometry ? 1 : 0);
        message = oss.str();
        return false;
    }
    std::ostringstream oss;
    oss << "nonuniform_flux_particles_ok rel_side=" << rel_side << " rel_bot=" << rel_bot
        << " P_side=" << side_w << "/" << p_side_geom << " P_bot=" << bottom_w << "/" << p_bot_geom
        << " normalize=" << (options.renormalize_boundary_area_to_geometry ? 1 : 0);
    message = oss.str();
    return true;
}

class OphelieThermalBoundaryComplianceReduceCK
    : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    OphelieThermalBoundaryComplianceReduceCK(SPHBody &sph_body, const std::string &temperature_field, Real t_wall,
                                             Real tolerance)
        : BaseDynamicsType(sph_body), t_wall_(t_wall), tolerance_(tolerance),
          dv_boundary_mask_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryMaskField)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
        quantity_name_ = "OphelieThermalBoundaryCompliance";
    }

    class FinishDynamics
    {
      public:
        using OutputType = Real;
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser)
        {
            (void)encloser;
        }
        Real Result(const ReduceReturnType &reduced_value)
        {
            return reduced_value.first / (reduced_value.second + TinyReal);
        }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_wall_(encloser.t_wall_), tolerance_(encloser.tolerance_),
              boundary_mask_(encloser.dv_boundary_mask_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (boundary_mask_[index_i] < Real(0.5))
            {
                return ReduceReturnType(0.0, 0.0);
            }
            const Real compliant = std::abs(temperature_[index_i] - t_wall_) <= tolerance_ ? 1.0 : 0.0;
            return ReduceReturnType(compliant, 1.0);
        }

      protected:
        Real t_wall_;
        Real tolerance_;
        Real *boundary_mask_;
        Real *temperature_;
    };

  protected:
    Real t_wall_;
    Real tolerance_;
    DiscreteVariable<Real> *dv_boundary_mask_;
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieThermalMaxTemperatureReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    OphelieThermalMaxTemperatureReduceCK(SPHBody &sph_body, const std::string &temperature_field)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
        quantity_name_ = "OphelieThermalMaxTemperature";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return temperature_[index_i];
        }

      protected:
        Real *temperature_;
    };

  protected:
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieThermalMinTemperatureReduceCK : public BaseLocalDynamicsReduce<ReduceMin, SPHBody>
{
  public:
    OphelieThermalMinTemperatureReduceCK(SPHBody &sph_body, const std::string &temperature_field)
        : BaseLocalDynamicsReduce<ReduceMin, SPHBody>(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
        quantity_name_ = "OphelieThermalMinTemperature";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return temperature_[index_i];
        }

      protected:
        Real *temperature_;
    };

  protected:
    DiscreteVariable<Real> *dv_temperature_;
};

/**
 * Explicit Euler diffusion using a precomputed pairwise Laplace that already includes k:
 *   L ≈ −∇·(k ∇T)  (pair form Σ C_ij (T_i − T_j) with C_ij ∝ harmonic-mean k)
 *   T += dt * (−1/(ρ cp)) * L
 * Do not multiply by k again here (Round 3A energy-audit fix).
 */
class ApplyLaplaceDiffusionDtCK : public LocalDynamics
{
  public:
    ApplyLaplaceDiffusionDtCK(SPHBody &sph_body, Real rho, Real cp, Real t_initial)
        : LocalDynamics(sph_body), diff_factor_(Real(-1) / (rho * cp + TinyReal)), t_initial_(t_initial),
          dv_laplace_(particles_->template getVariableByName<Real>(kOphelieThermalLaplaceTField)),
          dv_temperature_(particles_->template getVariableByName<Real>(kOphelieTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField))
    {
    }

    /** Deprecated overload: `k` is ignored because conductivity lives in the Laplace field. */
    ApplyLaplaceDiffusionDtCK(SPHBody &sph_body, Real rho, Real cp, Real /*k*/, Real t_initial)
        : ApplyLaplaceDiffusionDtCK(sph_body, rho, cp, t_initial)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : diff_factor_(encloser.diff_factor_), t_initial_(encloser.t_initial_),
              laplace_(encloser.dv_laplace_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            // Doc 45 F01: accumulate on ΔT (authoritative), then derive absolute T.
            // Adding a sub-ULP increment to T≈1473 K in float32 can round to zero.
            delta_t_[index_i] += dt * diff_factor_ * laplace_[index_i];
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real diff_factor_;
        Real t_initial_;
        Real *laplace_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real diff_factor_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_laplace_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

/**
 * Doc139: diffusion update with per-particle cp:
 *   ΔT += dt * (−1/(ρ cp_i)) * L ,   T = T0 + ΔT
 * Requires OphelieHeatCapacity field (see french_material_laws.h).
 */
class ApplyLaplaceDiffusionDtParticleCpCK : public LocalDynamics
{
  public:
    ApplyLaplaceDiffusionDtParticleCpCK(SPHBody &sph_body, Real rho, Real t_initial)
        : LocalDynamics(sph_body), inv_rho_(Real(1) / (rho + TinyReal)), t_initial_(t_initial),
          dv_laplace_(particles_->template getVariableByName<Real>(kOphelieThermalLaplaceTField)),
          dv_cp_(particles_->template getVariableByName<Real>("OphelieHeatCapacity")),
          dv_temperature_(particles_->template getVariableByName<Real>(kOphelieTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : inv_rho_(encloser.inv_rho_), t_initial_(encloser.t_initial_),
              laplace_(encloser.dv_laplace_->DelegatedData(ex_policy)),
              cp_(encloser.dv_cp_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            delta_t_[index_i] += dt * (Real(-1) * inv_rho_ / (cp_[index_i] + TinyReal)) * laplace_[index_i];
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real inv_rho_;
        Real t_initial_;
        Real *laplace_;
        Real *cp_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real inv_rho_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_laplace_;
    DiscreteVariable<Real> *dv_cp_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

/**
 * Doc139: French natural BC with per-particle cp:
 *   ΔT -= dt * q_out / (ρ cp_i δ)
 */
class ApplyOphelieThermalFrenchNaturalBcParticleCpCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalFrenchNaturalBcParticleCpCK(SPHBody &sph_body, const std::string &delta_t_field,
                                                   const std::string &temperature_field, Real t_initial, Real rho,
                                                   Real shell_thickness, const OphelieThermalDiffusionOneWayOptions &bc)
        : LocalDynamics(sph_body), t_initial_(t_initial), inv_rho_delta_(Real(1) / ((rho + TinyReal) * (shell_thickness + TinyReal))),
          h_side_(bc.h_side), h_bottom_(bc.h_bottom), h_free_(bc.h_free), emissivity_(bc.emissivity), t_cool_(bc.t_cool),
          t_ambient_(bc.t_ambient), t_rad_ambient_(bc.t_rad_ambient), stefan_boltzmann_(bc.stefan_boltzmann),
          enable_free_crust_(bc.enable_free_crust ? Real(1) : Real(0)), free_crust_delta_(bc.free_crust_delta_m),
          free_crust_k_(bc.free_crust_k),
          dv_w_side_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightSideField)),
          dv_w_bottom_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightBottomField)),
          dv_w_free_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryWeightFreeField)),
          dv_cp_(particles_->template getVariableByName<Real>("OphelieHeatCapacity")),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), inv_rho_delta_(encloser.inv_rho_delta_), h_side_(encloser.h_side_),
              h_bottom_(encloser.h_bottom_), h_free_(encloser.h_free_), emissivity_(encloser.emissivity_),
              t_cool_(encloser.t_cool_), t_ambient_(encloser.t_ambient_), t_rad_ambient_(encloser.t_rad_ambient_),
              stefan_boltzmann_(encloser.stefan_boltzmann_), enable_free_crust_(encloser.enable_free_crust_),
              free_crust_delta_(encloser.free_crust_delta_), free_crust_k_(encloser.free_crust_k_),
              w_side_(encloser.dv_w_side_->DelegatedData(ex_policy)),
              w_bottom_(encloser.dv_w_bottom_->DelegatedData(ex_policy)),
              w_free_(encloser.dv_w_free_->DelegatedData(ex_policy)), cp_(encloser.dv_cp_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real ws = w_side_[index_i];
            const Real wb = w_bottom_[index_i];
            const Real wf = w_free_[index_i];
            if (ws + wb + wf <= TinyReal)
            {
                return;
            }
            const Real t = t_initial_ + delta_t_[index_i];
            Real q_side = 0, q_bottom = 0, q_free = 0, q_conv = 0, q_rad = 0;
            const Real q_out = frenchCombinedBoundaryHeatFluxWm2(
                t, ws, wb, wf, h_side_, h_bottom_, t_cool_, h_free_, emissivity_, t_ambient_, t_rad_ambient_,
                stefan_boltzmann_, enable_free_crust_, free_crust_delta_, free_crust_k_, q_side, q_bottom, q_free,
                q_conv, q_rad);
            (void)q_conv;
            (void)q_rad;
            delta_t_[index_i] -= dt * inv_rho_delta_ * q_out / (cp_[index_i] + TinyReal);
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real inv_rho_delta_;
        Real h_side_;
        Real h_bottom_;
        Real h_free_;
        Real emissivity_;
        Real t_cool_;
        Real t_ambient_;
        Real t_rad_ambient_;
        Real stefan_boltzmann_;
        Real enable_free_crust_;
        Real free_crust_delta_;
        Real free_crust_k_;
        Real *w_side_;
        Real *w_bottom_;
        Real *w_free_;
        Real *cp_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    Real inv_rho_delta_;
    Real h_side_;
    Real h_bottom_;
    Real h_free_;
    Real emissivity_;
    Real t_cool_;
    Real t_ambient_;
    Real t_rad_ambient_;
    Real stefan_boltzmann_;
    Real enable_free_crust_;
    Real free_crust_delta_;
    Real free_crust_k_;
    DiscreteVariable<Real> *dv_w_side_;
    DiscreteVariable<Real> *dv_w_bottom_;
    DiscreteVariable<Real> *dv_w_free_;
    DiscreteVariable<Real> *dv_cp_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

template <class ExecutionPolicy>
inline void execOphelieJouleHeatDiffusionOneWayStep(SolidBody &glass_body, Inner<> &glass_inner,
                                                    const std::string &q_field, const std::string &temperature_field,
                                                    Real t_initial, Real dt, const OphelieJouleHeatOneWayMaterialProps &material,
                                                    const OphelieThermalDiffusionOneWayOptions &options)
{
    if (options.enable_diffusion)
    {
        UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
        UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace_temperature(
            glass_inner, temperature_field, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
            options.pair_weight_regularization);
        laplace_temperature.exec();
    }
    else
    {
        StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_laplace(glass_body, kOphelieThermalLaplaceTField);
        zero_laplace.exec();
    }

    StateDynamics<ExecutionPolicy, ApplyOphelieJouleHeatDiffusionCombinedStepCK> combined_step(
        glass_body, q_field, kOphelieThermalDeltaTField, temperature_field, kOphelieThermalLaplaceTField, t_initial, dt,
        material.rho, material.cp, material.k, options.enable_diffusion);
    combined_step.exec();

    if (options.enable_cold_wall_dirichlet)
    {
        StateDynamics<ExecutionPolicy, ApplyOphelieThermalDirichletWallCK> dirichlet_wall(
            glass_body, kOphelieThermalDeltaTField, temperature_field, t_initial);
        dirichlet_wall.exec();
    }
    else if (options.enable_french_natural_bc)
    {
        const Real shell_thickness =
            options.boundary_shell_thickness > TinyReal
                ? options.boundary_shell_thickness
                : options.boundary_width_factor; // fallback; prefer explicit setup
        StateDynamics<ExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> natural_bc(
            glass_body, kOphelieThermalDeltaTField, temperature_field, t_initial, material.rho, material.cp,
            shell_thickness, options);
        natural_bc.exec(dt);
    }
}

template <class ExecutionPolicy>
inline OphelieThermalDiffusionOneWayStepResult execOphelieThermalDiffusionDiagnostics(
    SolidBody &glass_body, const std::string &q_field, Real effective_dt, Real rho, Real cp, size_t n,
    Real t_wall, Real q_threshold, const OphelieThermalDiffusionOneWayOptions &options)
{
    OphelieThermalDiffusionOneWayStepResult result;
    const OphelieJouleHeatOneWayStepResult base =
        execOphelieJouleHeatOneWayClosureDiagnostics<ExecutionPolicy>(
            glass_body, q_field, kOphelieThermalDeltaTField, effective_dt, rho, cp, n, q_threshold);
    result.n_steps = base.n_steps;
    result.max_per_particle_rel_err = base.max_per_particle_rel_err;
    result.mean_delta_t = base.mean_delta_t;
    result.max_delta_t = base.max_delta_t;
    result.total_joule_energy_j = base.total_joule_energy_j;
    result.total_thermal_energy_j = base.total_thermal_energy_j;
    result.energy_balance_rel_err = base.energy_balance_rel_err;
    result.vol_weighted_delta_t = base.vol_weighted_delta_t;
    result.vol_weighted_expected_delta_t = base.vol_weighted_expected_delta_t;
    result.closure_mismatch_vol_fraction = base.closure_mismatch_vol_fraction;
    result.closure_inline_energy_rel_err = base.closure_inline_energy_rel_err;
    if (options.enable_cold_wall_dirichlet)
    {
        ReduceDynamicsCK<ExecutionPolicy, OphelieThermalBoundaryComplianceReduceCK> boundary_compliance(
            glass_body, kOphelieTemperatureField, t_wall, std::max(Real(1.0e-6), Real(1.0e-4) * std::abs(t_wall)));
        result.boundary_dirichlet_compliance = boundary_compliance.exec();
    }
    else
    {
        result.boundary_dirichlet_compliance = 1.0;
    }
    if (options.enable_french_natural_bc)
    {
        const Real shell =
            options.boundary_shell_thickness > TinyReal ? options.boundary_shell_thickness : options.boundary_width_factor;
        hostOphelieThermalFrenchNaturalHeatLossPowers(glass_body.getBaseParticles(), options, shell,
                                                      result.wall_heat_loss_side_w, result.wall_heat_loss_bottom_w,
                                                      result.free_surface_conv_loss_w, result.free_surface_rad_loss_w,
                                                      result.total_heat_loss_w);
    }
    ReduceDynamicsCK<ExecutionPolicy, OphelieThermalMaxTemperatureReduceCK> max_temperature(
        glass_body, kOphelieTemperatureField);
    result.max_temperature = max_temperature.exec();
    return result;
}

inline void copyOphelieJouleHeatOneWayStepResult(OphelieJouleHeatOneWayStepResult &dst,
                                                 const OphelieJouleHeatOneWayStepResult &src)
{
    dst.n_steps = src.n_steps;
    dst.max_per_particle_rel_err = src.max_per_particle_rel_err;
    dst.mean_delta_t = src.mean_delta_t;
    dst.max_delta_t = src.max_delta_t;
    dst.total_joule_energy_j = src.total_joule_energy_j;
    dst.total_thermal_energy_j = src.total_thermal_energy_j;
    dst.energy_balance_rel_err = src.energy_balance_rel_err;
    dst.vol_weighted_delta_t = src.vol_weighted_delta_t;
    dst.vol_weighted_expected_delta_t = src.vol_weighted_expected_delta_t;
    dst.closure_mismatch_vol_fraction = src.closure_mismatch_vol_fraction;
    dst.closure_inline_energy_rel_err = src.closure_inline_energy_rel_err;
    dst.boundary_dirichlet_compliance = src.boundary_dirichlet_compliance;
    dst.max_temperature = src.max_temperature;
    dst.wall_heat_loss_side_w = src.wall_heat_loss_side_w;
    dst.wall_heat_loss_bottom_w = src.wall_heat_loss_bottom_w;
    dst.free_surface_conv_loss_w = src.free_surface_conv_loss_w;
    dst.free_surface_rad_loss_w = src.free_surface_rad_loss_w;
    dst.total_heat_loss_w = src.total_heat_loss_w;
}

template <class ExecutionPolicy>
inline OphelieThermalDiffusionOneWayStepResult applyOphelieJouleHeatDiffusionOneWaySteps(
    SolidBody &glass_body, Inner<> &glass_inner, BaseParticles &particles, const std::string &q_field,
    const std::string &temperature_field, Real dt, const OphelieJouleHeatOneWayMaterialProps &material, size_t n,
    size_t n_steps, OphelieThermalDiffusionOneWayOptions &options, Real q_threshold = 1.0e-6,
    const OphelieThermalVtpRecordingOptions *recording = nullptr)
{
    OphelieThermalDiffusionOneWayStepResult result;
    result.n_steps = n_steps;
    const Real effective_dt = dt * static_cast<Real>(n_steps);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, q_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalDeltaTField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, temperature_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalLaplaceTField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalConductivityField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryMaskField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryFaceField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryWeightSideField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, "VolumetricMeasure");
    syncOphelieJouleHeatOneWayThermalFieldsToDevice(particles, kOphelieThermalDeltaTField, temperature_field);

    for (size_t step = 0; step < n_steps; ++step)
    {
        execOphelieJouleHeatDiffusionOneWayStep<ExecutionPolicy>(glass_body, glass_inner, q_field, temperature_field,
                                                                 material.t_initial, dt, material, options);
        writeOphelieThermalVtpIfDue(recording, step, n_steps);
    }

    result = execOphelieThermalDiffusionDiagnostics<ExecutionPolicy>(
        glass_body, q_field, effective_dt, material.rho, material.cp, n, material.t_initial, q_threshold, options);
    return result;
}

class ResetOphelieThermalStateToInitialCK : public LocalDynamics
{
  public:
    ResetOphelieThermalStateToInitialCK(SPHBody &sph_body, Real t_initial,
                                        const std::string &delta_t_field = kOphelieThermalDeltaTField,
                                        const std::string &temperature_field = kOphelieTemperatureField)
        : LocalDynamics(sph_body), t_initial_(t_initial),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            delta_t_[index_i] = Real(0);
            temperature_[index_i] = t_initial_;
        }

      protected:
        Real t_initial_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

/**
 * Stage 3.3 verification helper:
 * keep frozen EM/Q, reset T→T0 (clear σ-seed heating history), then run diffusion + cold wall.
 * Used after σ(T) coupling so energy_cap / Dirichlet gates stay auditable.
 */
template <class ExecutionPolicy>
inline OphelieThermalDiffusionOneWayStepResult applyOphelieFrozenQDiffusionFromUniformT0(
    SolidBody &glass_body, Inner<> &glass_inner, BaseParticles &particles, const std::string &q_field,
    const OphelieFrenchReducedCaseParams &french, Real thermal_dt, size_t thermal_steps,
    const OphelieJouleHeatOneWayMaterialProps &material, OphelieThermalDiffusionOneWayOptions &options,
    const OphelieThermalVtpRecordingOptions *recording = nullptr)
{
    registerOphelieJouleHeatTemperatureField(particles, material.t_initial);
    registerOphelieThermalDiffusionAuxFields(particles, material.k);
    StateDynamics<ExecutionPolicy, ResetOphelieThermalStateToInitialCK> reset_thermal(glass_body, material.t_initial);
    reset_thermal.exec();

    if (options.enable_cold_wall_dirichlet)
    {
        OpheliePhiBoundaryGeometryContext geom;
        geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
        geom.french = french;
        (void)setupOphelieThermalDirichletBoundaryMask(particles, options, geom, french.dp);
    }
    else if (options.enable_french_natural_bc)
    {
        (void)setupOphelieThermalFrenchNaturalBoundaryFaces(particles, options, french);
    }

    const size_t n = particles.TotalRealParticles();
    return applyOphelieJouleHeatDiffusionOneWaySteps<ExecutionPolicy>(
        glass_body, glass_inner, particles, q_field, kOphelieTemperatureField, thermal_dt, material, n, thermal_steps,
        options, 1.0e-6, recording);
}

/** French EM + Joule heat with optional isotropic diffusion and thermal BC. */
template <class ExecutionPolicy>
inline OphelieFrenchEmJouleHeatOneWayResult runFrenchReducedEmThenJouleHeatDiffusionOneWay(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french, Real thermal_dt, size_t thermal_steps,
    const OphelieJouleHeatOneWayMaterialProps &material, OphelieThermalDiffusionOneWayOptions &thermal_options,
    const OphelieThermalVtpRecordingOptions *recording = nullptr)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    OphelieFrenchEmJouleHeatOneWayResult result;
    runFrenchReducedEmOrSelfInductionForThermalHandoff<ExecutionPolicy>(glass_body, glass_inner, names, params, french,
                                                                        result);

    const size_t n = particles.TotalRealParticles();
    registerOphelieJouleHeatTemperatureField(particles, material.t_initial);
    registerOphelieThermalDiffusionAuxFields(particles, material.k);

    if (thermal_options.enable_cold_wall_dirichlet)
    {
        OpheliePhiBoundaryGeometryContext geom;
        geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
        geom.french = french;
        (void)setupOphelieThermalDirichletBoundaryMask(particles, thermal_options, geom, french.dp);
    }
    else if (thermal_options.enable_french_natural_bc)
    {
        (void)setupOphelieThermalFrenchNaturalBoundaryFaces(particles, thermal_options, french);
    }

    const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
    const OphelieThermalDiffusionOneWayStepResult thermal_diffusion =
        applyOphelieJouleHeatDiffusionOneWaySteps<ExecutionPolicy>(
            glass_body, glass_inner, particles, q_field, kOphelieTemperatureField, thermal_dt, material, n,
            thermal_steps, thermal_options, 1.0e-6, recording);
    copyOphelieJouleHeatOneWayStepResult(result.thermal, thermal_diffusion);
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H
