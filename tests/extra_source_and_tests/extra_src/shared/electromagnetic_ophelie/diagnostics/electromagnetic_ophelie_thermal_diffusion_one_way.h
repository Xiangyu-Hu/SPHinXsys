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
    /** Filled by setup helpers: boundary_width_factor * dp. */
    Real boundary_shell_thickness = Real(0);
};

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

struct OphelieThermalDiffusionOneWayStepResult : public OphelieJouleHeatOneWayStepResult
{
};

inline void registerOphelieThermalDiffusionAuxFields(BaseParticles &particles, Real thermal_conductivity)
{
    particles.registerStateVariable<Real>(kOphelieThermalLaplaceTField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalConductivityField, thermal_conductivity);
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryMaskField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryFaceField, Real(0));
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
 * Host: classify French cylinder shell into side / bottom / free-top for Natural production BC.
 * Face codes: 0 interior, 1 side, 2 bottom, 3 free-top.
 */
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
    size_t n_boundary = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real dist_side = french.glass_radius - r;
        const Real dist_bottom = pos[i][2] - z_lo;
        const Real dist_top = z_hi - pos[i][2];
        const Real dist = std::min(dist_side, std::min(dist_bottom, dist_top));
        if (dist > boundary_width)
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
    return n_boundary;
}

class ApplyOphelieJouleHeatDiffusionCombinedStepCK : public LocalDynamics
{
  public:
    ApplyOphelieJouleHeatDiffusionCombinedStepCK(SPHBody &sph_body, const std::string &q_field,
                                                 const std::string &delta_t_field,
                                                 const std::string &temperature_field,
                                                 const std::string &laplace_t_field, Real t_initial, Real dt,
                                                 Real rho, Real cp, Real k, bool apply_diffusion)
        : LocalDynamics(sph_body), t_initial_(t_initial), joule_factor_(dt / (rho * cp + TinyReal)),
          diffusion_factor_(apply_diffusion ? -dt * k / (rho * cp + TinyReal) : Real(0)),
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
          dv_face_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryFaceField)),
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
              stefan_boltzmann_(encloser.stefan_boltzmann_), face_(encloser.dv_face_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real face = face_[index_i];
            if (face < Real(0.5))
            {
                return;
            }
            const Real t = temperature_[index_i];
            Real q_out = Real(0);
            if (face < Real(1.5))
            {
                q_out = h_side_ * (t - t_cool_);
            }
            else if (face < Real(2.5))
            {
                q_out = h_bottom_ * (t - t_cool_);
            }
            else
            {
                const Real t2 = t * t;
                const Real tar2 = t_rad_ambient_ * t_rad_ambient_;
                q_out = h_free_ * (t - t_ambient_) +
                        emissivity_ * stefan_boltzmann_ * (t2 * t2 - tar2 * tar2);
            }
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
        Real *face_;
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
    DiscreteVariable<Real> *dv_face_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

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
    syncVariableToHost<Real>(particles, kOphelieThermalBoundaryFaceField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *face = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryFaceField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real inv_delta = Real(1) / (shell_thickness + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real f = face[i];
        if (f < Real(0.5))
        {
            continue;
        }
        const Real area = vol[i] * inv_delta;
        const Real t = temperature[i];
        if (f < Real(1.5))
        {
            const Real q = options.h_side * (t - options.t_cool);
            side_w += q * area;
        }
        else if (f < Real(2.5))
        {
            const Real q = options.h_bottom * (t - options.t_cool);
            bottom_w += q * area;
        }
        else
        {
            const Real t2 = t * t;
            const Real tar2 = options.t_rad_ambient * options.t_rad_ambient;
            const Real q_conv = options.h_free * (t - options.t_ambient);
            const Real q_rad =
                options.emissivity * options.stefan_boltzmann * (t2 * t2 - tar2 * tar2);
            free_conv_w += q_conv * area;
            free_rad_w += q_rad * area;
        }
    }
    total_w = side_w + bottom_w + free_conv_w + free_rad_w;
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
