#ifndef APHI_SOURCE_DRIVEN_EM_SOLVE_HELPERS_H
#define APHI_SOURCE_DRIVEN_EM_SOLVE_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_cold_crucible_demo_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include "io_vtk.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Source-driven EM solve (GMRES), not impressed-A. Uses TEAM7-like single-body region tags. */
struct AphiSourceDrivenEmSolveSpec
{
    Real body_length = benchmark::AphiTeam7CanonicalCaseSpec::body_length;
    Real body_height = benchmark::AphiTeam7CanonicalCaseSpec::body_height;
    Real body_width = benchmark::AphiTeam7CanonicalCaseSpec::body_width;

    Real dp = benchmark::AphiTeam7CanonicalCaseSpec::candidate_dp;
    Real omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    Real phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    Real tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    UnsignedInt restart_dimension = benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension;
    UnsignedInt max_outer_iterations = benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations;

    Real min_solution_norm = 1.0e-8;
    Real min_conductor_E = 1.0e-10;
    Real min_conductor_J = 1.0e-12;
    Real min_conductor_joule_integral = 1.0e-14;
    /** Air sigma is tiny but nonzero; require air Joule ≪ conductor Joule. */
    Real max_air_to_conductor_joule_ratio = 0.05;
    /** Scales outer SPH padding only; TEAM7 interior box size is unchanged. */
    Real boundary_width_scale = 1.0;
    /** Passive sigma=0 support shell outside physical box (Strategy C diagnostic only). Width = scale * dp. */
    Real passive_air_shell_width_scale = 0.0;
    /** @deprecated legacy alias for passive_air_shell_width_scale */
    Real dummy_shell_width_scale = 0.0;
    bool write_vtp = true;
    bool write_probe_csv = false;
    UnsignedInt probe_centerline_samples = 32;
    std::string probe_output_dir = "output";
};

struct AphiSourceDrivenEmSolveMetrics
{
    bool converged = false;
    UnsignedInt num_iterations = 0;
    Real final_residual = 0.0;
    Real max_abs_A = 0.0;
    Real max_abs_phi = 0.0;
    Real max_abs_B = 0.0;
    Real max_abs_H = 0.0;
    Real max_abs_E = 0.0;
    Real max_abs_J = 0.0;
    Real max_abs_J_conductor = 0.0;
    Real max_Joule = 0.0;
    Real max_Joule_conductor = 0.0;
    Real air_Joule_integral = 0.0;
    Real conductor_Joule_integral = 0.0;
    Real source_Joule_integral = 0.0;
    Real physical_box_Joule_integral = 0.0;
    Real total_Joule_integral = 0.0;
    Real source_rhs_l2 = 0.0;
    Real source_rhs_l2_source_region = 0.0;
    size_t particle_count_physical_box = 0;
    size_t particle_count_shell = 0;
    size_t particle_count_conductor = 0;
    size_t particle_count_source = 0;
    size_t particle_count_air = 0;
    Real shell_sigma_max = 0.0;
    Real shell_joule_integral = 0.0;
    bool finite_field_check = false;
    size_t particles = 0;
    std::string vtp_output_dir = "output";
    AphiProbeCsvWriteResult probe_csv{};
    bool probe_csv_requested = false;
};

struct AphiSourceDrivenEmSolveFieldNames
{
    std::string b_real = "BReal";
    std::string b_imag = "BImag";
    std::string h_real = "HReal";
    std::string h_imag = "HImag";
    std::string a_magnitude = "AMagnitude";
    std::string b_magnitude = "BMagnitude";
    std::string e_magnitude = "EMagnitude";
    std::string j_magnitude = "JMagnitude";
    std::string h_magnitude = "HMagnitude";
    std::string material_region_id = "MaterialRegionId";
};

inline Real hostParticleRegionBlockNorm(BaseParticles &particles, const AphiBlockNames &block_names,
                                        const Vecd *positions, size_t total_real_particles,
                                        const std::function<bool(const Vecd &)> &in_region)
{
    syncAphiBlockToHost(particles, block_names);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_names.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_names.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        const Real block_squared = a_real[i].squaredNorm() + a_imag[i].squaredNorm() + phi_real[i] * phi_real[i] +
                                   phi_imag[i] * phi_imag[i];
        sum_squared += vol[i] * block_squared;
    }
    return std::sqrt(sum_squared);
}

inline Real effectivePassiveAirShellWidthScale(const AphiSourceDrivenEmSolveSpec &spec)
{
    if (spec.passive_air_shell_width_scale > 0.0)
    {
        return spec.passive_air_shell_width_scale;
    }
    return spec.dummy_shell_width_scale;
}

/** Prescribed-current source region: non-conductive, RHS-only (sigma=0 in coil box). */
class AssignZeroSigmaInTeam7CoilRegionCK : public LocalDynamics
{
  public:
    AssignZeroSigmaInTeam7CoilRegionCK(SPHBody &sph_body, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                       const AphiMaterialNames &material_names)
        : LocalDynamics(sph_body), layout_(layout),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : layout_(encloser.layout_), position_(encloser.dv_position_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (team7ParticleInRegion(position_[index_i], layout_, AphiBenchmarkMaterialRegion::Coil))
            {
                sigma_[index_i] = 0.0;
            }
        }

      protected:
        benchmark::AphiTeam7LikeUnitBoxLayout layout_;
        Vecd *position_;
        Real *sigma_;
    };

  protected:
    benchmark::AphiTeam7LikeUnitBoxLayout layout_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
};

class AssignZeroSigmaOutsidePhysicalBoxCK : public LocalDynamics
{
  public:
    AssignZeroSigmaOutsidePhysicalBoxCK(SPHBody &sph_body, Real body_length, Real body_height, Real body_width,
                                        const AphiMaterialNames &material_names)
        : LocalDynamics(sph_body), body_length_(body_length), body_height_(body_height), body_width_(body_width),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : body_length_(encloser.body_length_), body_height_(encloser.body_height_), body_width_(encloser.body_width_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (!isInsideTeam7PhysicalBox(position_[index_i], body_length_, body_height_, body_width_))
            {
                sigma_[index_i] = 0.0;
            }
        }

      protected:
        Real body_length_;
        Real body_height_;
        Real body_width_;
        Vecd *position_;
        Real *sigma_;
    };

  protected:
    Real body_length_;
    Real body_height_;
    Real body_width_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
};

class AssignTeam7MaterialRegionIdCK : public LocalDynamics
{
  public:
    AssignTeam7MaterialRegionIdCK(SPHBody &sph_body, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                  const std::string &region_id_name)
        : LocalDynamics(sph_body), layout_(layout),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_region_id_(particles_->template registerStateVariable<Real>(region_id_name, Real(0)))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              region_id_(encloser.dv_region_id_->DelegatedData(ex_policy)), layout_(encloser.layout_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            if (team7ParticleInRegion(position, layout_, AphiBenchmarkMaterialRegion::Conductor))
            {
                region_id_[index_i] = 1.0;
            }
            else if (team7ParticleInRegion(position, layout_, AphiBenchmarkMaterialRegion::Coil))
            {
                region_id_[index_i] = 2.0;
            }
            else
            {
                region_id_[index_i] = 0.0;
            }
        }

      protected:
        Vecd *position_;
        Real *region_id_;
        benchmark::AphiTeam7LikeUnitBoxLayout layout_;
    };

  protected:
    benchmark::AphiTeam7LikeUnitBoxLayout layout_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_region_id_;
};

inline void registerSourceDrivenEmSolveObservableFields(SPHBody &body, const AphiSourceDrivenEmSolveFieldNames &fields)
{
    BaseParticles &particles = body.getBaseParticles();
    particles.registerStateVariable<Vecd>(fields.b_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(fields.b_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(fields.h_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(fields.h_imag, ZeroData<Vecd>::value);
}

inline void writeSourceDrivenEmSolveVtp(SPHSystem &sph_system, SPHBody &body, const AphiVariableNames &names,
                                        const AphiJouleHeatingFieldNames &joule_fields,
                                        const AphiSourceDrivenEmSolveFieldNames &obs_fields,
                                        const AphiMaterialNames &material_names)
{
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(body, names.solution.a_real);
    write_states.addToWrite<Vecd>(body, names.solution.a_imag);
    write_states.addToWrite<Real>(body, names.solution.phi_real);
    write_states.addToWrite<Real>(body, names.solution.phi_imag);
    write_states.addToWrite<Vecd>(body, obs_fields.b_real);
    write_states.addToWrite<Vecd>(body, obs_fields.b_imag);
    write_states.addToWrite<Vecd>(body, obs_fields.h_real);
    write_states.addToWrite<Vecd>(body, obs_fields.h_imag);
    write_states.addToWrite<Vecd>(body, joule_fields.electric_field_a_real);
    write_states.addToWrite<Vecd>(body, joule_fields.electric_field_a_imag);
    write_states.addToWrite<Vecd>(body, joule_fields.current_density_real);
    write_states.addToWrite<Vecd>(body, joule_fields.current_density_imag);
    write_states.addToWrite<Real>(body, joule_fields.joule_heat_source);
    write_states.addToWrite<Real>(body, material_names.sigma);
    write_states.addToWrite<Real>(body, material_names.nu);
    write_states.addToWrite<Real>(body, obs_fields.material_region_id);
    write_states.addToWrite<Real>(body, obs_fields.a_magnitude);
    write_states.addToWrite<Real>(body, obs_fields.b_magnitude);
    write_states.addToWrite<Real>(body, obs_fields.e_magnitude);
    write_states.addToWrite<Real>(body, obs_fields.j_magnitude);
    write_states.addToWrite<Real>(body, obs_fields.h_magnitude);
    write_states.writeToFile(0);
}

inline bool hostVecdFieldFinite(BaseParticles &particles, const std::string &name, size_t total_real_particles)
{
    syncVariableToHost<Vecd>(particles, name);
    const Vecd *values = particles.getVariableDataByName<Vecd>(name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!values[i].allFinite())
        {
            return false;
        }
    }
    return true;
}

inline bool hostScalarFieldFinite(BaseParticles &particles, const std::string &name, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, name);
    const Real *values = particles.getVariableDataByName<Real>(name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!std::isfinite(values[i]))
        {
            return false;
        }
    }
    return true;
}

inline bool hostSourceDrivenFieldsFinite(BaseParticles &particles, const AphiVariableNames &names,
                                         const AphiJouleHeatingFieldNames &joule_fields,
                                         const AphiSourceDrivenEmSolveFieldNames &obs_fields, size_t total_real_particles)
{
    return hostVecdFieldFinite(particles, names.solution.a_real, total_real_particles) &&
           hostVecdFieldFinite(particles, names.solution.a_imag, total_real_particles) &&
           hostScalarFieldFinite(particles, names.solution.phi_real, total_real_particles) &&
           hostScalarFieldFinite(particles, names.solution.phi_imag, total_real_particles) &&
           hostVecdFieldFinite(particles, obs_fields.b_real, total_real_particles) &&
           hostVecdFieldFinite(particles, joule_fields.electric_field_a_real, total_real_particles) &&
           hostVecdFieldFinite(particles, joule_fields.current_density_real, total_real_particles) &&
           hostScalarFieldFinite(particles, joule_fields.joule_heat_source, total_real_particles);
}

/** GMRES + E/J/Joule on an existing test body (shared by P1 and P2). */
inline AphiMatrixFreeSolverResult execSourceDrivenEmJoulePipelineOnBody(
    AphiLhsTestBody &test_body, const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const AphiSourceDrivenEmSolveSpec &spec,
    const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_fields,
    AphiSourceDrivenEmSolveFieldNames *obs_fields = nullptr, Real passive_air_shell_width = 0.0)
{
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    if (obs_fields != nullptr)
    {
        registerSourceDrivenEmSolveObservableFields(test_body.body, *obs_fields);
    }
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_material(
        test_body.body, layout, names.material);
    if (obs_fields != nullptr)
    {
        StateDynamics<MainExecutionPolicy, AssignTeam7MaterialRegionIdCK> assign_region_id(
            test_body.body, layout, obs_fields->material_region_id);
        assign_region_id.exec();
    }
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, spec.impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = spec.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = spec.phi_gauge_penalty;
    options.use_a_divergence_penalty = false;
    options.a_divergence_penalty = 0.0;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(spec.tolerance, spec.restart_dimension, spec.max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options, solver_options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    initialize_aphi_variables.exec();
    assign_material.exec();
    StateDynamics<MainExecutionPolicy, AssignZeroSigmaInTeam7CoilRegionCK> assign_source_only_coil(
        test_body.body, layout, names.material);
    assign_source_only_coil.exec();
    if (passive_air_shell_width > 0.0)
    {
        StateDynamics<MainExecutionPolicy, AssignZeroSigmaOutsidePhysicalBoxCK> assign_passive_air_shell(
            test_body.body, spec.body_length, spec.body_height, spec.body_width, names.material);
        assign_passive_air_shell.exec();
    }
    assign_coil_source.exec();
    test_body.updateRelations();
    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult em_result = solver.solve();

    if (obs_fields != nullptr)
    {
        execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, obs_fields->b_real, obs_fields->b_imag,
                                     AphiBCurlDiagnosticMode::BCorrectedGrad);
        StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByNuCK> compute_h_real(
            test_body.body, obs_fields->b_real, obs_fields->h_real, names.material);
        StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByNuCK> compute_h_imag(
            test_body.body, obs_fields->b_imag, obs_fields->h_imag, names.material);
        compute_h_real.exec();
        compute_h_imag.exec();
    }
    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    if (obs_fields != nullptr)
    {
        StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> a_magnitude(
            test_body.body, names.solution.a_real, names.solution.a_imag, obs_fields->a_magnitude);
        StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> b_magnitude(
            test_body.body, obs_fields->b_real, obs_fields->b_imag, obs_fields->b_magnitude);
        StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> h_magnitude(
            test_body.body, obs_fields->h_real, obs_fields->h_imag, obs_fields->h_magnitude);
        StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> e_magnitude(
            test_body.body, joule_fields.electric_field_a_real, joule_fields.electric_field_a_imag,
            obs_fields->e_magnitude);
        StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> j_magnitude(
            test_body.body, joule_fields.current_density_real, joule_fields.current_density_imag,
            obs_fields->j_magnitude);
        a_magnitude.exec();
        b_magnitude.exec();
        h_magnitude.exec();
        e_magnitude.exec();
        j_magnitude.exec();
    }
    test_body.updateRelations();
    return em_result;
}

inline AphiSourceDrivenEmSolveMetrics runSourceDrivenEmSolveWithLayout(
    int ac, char *av[], const AphiSourceDrivenEmSolveSpec &spec,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout)
{
    const Real boundary_width = 3.0 * spec.dp * spec.boundary_width_scale;
    const Real passive_air_shell_width =
        effectivePassiveAirShellWidthScale(spec) > 0.0 ? effectivePassiveAirShellWidthScale(spec) * spec.dp : 0.0;

    AphiJouleHeatingFieldNames joule_fields;
    AphiSourceDrivenEmSolveFieldNames obs_fields;
    AphiVariableNames names;

    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av,
                              passive_air_shell_width);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    (void)register_joule_fields;
    const AphiMatrixFreeSolverResult em_result = execSourceDrivenEmJoulePipelineOnBody(
        test_body, layout, spec, names, joule_fields, &obs_fields, passive_air_shell_width);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiSourceDrivenEmSolveMetrics metrics;
    metrics.particles = total_real_particles;
    metrics.num_iterations = em_result.outer_iteration_count;
    metrics.final_residual = em_result.final_true_relative_residual;
    metrics.converged = gmresConvergencePassed(em_result, spec.tolerance);

    const auto all_particles = [](const Vecd &) { return true; };
    const auto in_physical_conductor = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Conductor);
    };
    const auto in_physical_air = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Air);
    };
    const auto in_physical_source = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Coil);
    };
    metrics.max_abs_A = hostParticleRegionScalarMax(particles, obs_fields.a_magnitude, positions, total_real_particles,
                                                    all_particles);
    metrics.max_abs_B = hostParticleRegionScalarMax(particles, obs_fields.b_magnitude, positions, total_real_particles,
                                                    all_particles);
    metrics.max_abs_H = hostParticleRegionScalarMax(particles, obs_fields.h_magnitude, positions, total_real_particles,
                                                    all_particles);
    metrics.max_abs_E = hostParticleRegionScalarMax(particles, obs_fields.e_magnitude, positions, total_real_particles,
                                                    all_particles);
    metrics.max_abs_J = hostParticleRegionScalarMax(particles, obs_fields.j_magnitude, positions, total_real_particles,
                                                    all_particles);
    metrics.max_abs_J_conductor = hostParticleRegionScalarMax(particles, obs_fields.j_magnitude, positions,
                                                              total_real_particles, in_physical_conductor);
    metrics.max_Joule = hostParticleRegionScalarMax(particles, joule_fields.joule_heat_source, positions,
                                                    total_real_particles, all_particles);
    metrics.max_Joule_conductor = hostParticleRegionScalarMax(particles, joule_fields.joule_heat_source, positions,
                                                              total_real_particles, in_physical_conductor);

    syncVariableToHost<Real>(particles, names.solution.phi_real);
    syncVariableToHost<Real>(particles, names.solution.phi_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(names.solution.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(names.solution.phi_imag);
    metrics.max_abs_phi = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        metrics.max_abs_phi =
            std::max(metrics.max_abs_phi, std::sqrt(phi_real[i] * phi_real[i] + phi_imag[i] * phi_imag[i]));
    }

    metrics.conductor_Joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_conductor, joule_fields.joule_heat_source);
    metrics.source_Joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_source, joule_fields.joule_heat_source);
    metrics.air_Joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                         in_physical_air, joule_fields.joule_heat_source);
    const auto in_physical_box = [&](const Vecd &position) {
        return !isInsidePassiveAirShellParticle(position, spec.body_length, spec.body_height, spec.body_width);
    };
    metrics.physical_box_Joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_box, joule_fields.joule_heat_source);
    metrics.total_Joule_integral = hostVolWeightedSum(particles, joule_fields.joule_heat_source, total_real_particles,
                                                      [](const Vecd &) { return true; });
    metrics.source_rhs_l2 = hostBlockNorm(particles, names.rhs, total_real_particles);
    metrics.source_rhs_l2_source_region =
        hostParticleRegionBlockNorm(particles, names.rhs, positions, total_real_particles, in_physical_source);

    syncVariableToHost<Real>(particles, names.material.sigma);
    const Real *sigma = particles.getVariableDataByName<Real>(names.material.sigma);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (isInsidePassiveAirShellParticle(positions[i], spec.body_length, spec.body_height, spec.body_width))
        {
            metrics.particle_count_shell += 1;
            metrics.shell_sigma_max = std::max(metrics.shell_sigma_max, sigma[i]);
            continue;
        }
        metrics.particle_count_physical_box += 1;
        if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
        {
            metrics.particle_count_conductor += 1;
        }
        else if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Coil))
        {
            metrics.particle_count_source += 1;
        }
        else
        {
            metrics.particle_count_air += 1;
        }
    }
    const auto in_shell = [&](const Vecd &position) {
        return isInsidePassiveAirShellParticle(position, spec.body_length, spec.body_height, spec.body_width);
    };
    metrics.shell_joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_shell, joule_fields.joule_heat_source);

    metrics.finite_field_check =
        hostSourceDrivenFieldsFinite(particles, names, joule_fields, obs_fields, total_real_particles);

    if (spec.write_vtp)
    {
        writeSourceDrivenEmSolveVtp(test_body.sph_system, test_body.body, names, joule_fields, obs_fields, names.material);
    }

    if (spec.write_probe_csv)
    {
        AphiProbeObservedFieldNames probe_fields;
        probe_fields.b_real = obs_fields.b_real;
        probe_fields.b_imag = obs_fields.b_imag;
        probe_fields.b_magnitude = obs_fields.b_magnitude;
        probe_fields.e_magnitude = obs_fields.e_magnitude;
        probe_fields.j_magnitude = obs_fields.j_magnitude;
        probe_fields.material_region_id = obs_fields.material_region_id;
        metrics.probe_csv_requested = true;
        metrics.probe_csv = writeSourceDrivenProbeCsvArtifacts(
            particles, layout, spec.body_length, spec.body_height, spec.body_width, probe_fields, joule_fields,
            names.material, metrics.converged, metrics.num_iterations, metrics.final_residual, spec.omega, spec.dp,
            metrics.particles, metrics.conductor_Joule_integral, metrics.air_Joule_integral, metrics.total_Joule_integral,
            spec.probe_output_dir, spec.probe_centerline_samples);
    }

    return metrics;
}

inline AphiSourceDrivenEmSolveMetrics runSourceDrivenEmSolve(int ac, char *av[], const AphiSourceDrivenEmSolveSpec &spec)
{
    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    return runSourceDrivenEmSolveWithLayout(ac, av, spec, layout);
}

inline bool sourceDrivenEmSolvePassed(const AphiSourceDrivenEmSolveMetrics &metrics,
                                      const AphiSourceDrivenEmSolveSpec &spec)
{
    const Real air_to_conductor_joule_ratio =
        metrics.air_Joule_integral / (metrics.conductor_Joule_integral + TinyReal);
    return metrics.converged && metrics.finite_field_check && metrics.max_abs_A > spec.min_solution_norm &&
           metrics.max_abs_E > spec.min_conductor_E && metrics.max_abs_J_conductor > spec.min_conductor_J &&
           metrics.max_Joule_conductor > 0.0 && metrics.conductor_Joule_integral > spec.min_conductor_joule_integral &&
           air_to_conductor_joule_ratio <= spec.max_air_to_conductor_joule_ratio &&
           metrics.source_rhs_l2_source_region > 0.0 && metrics.particle_count_source > 0 &&
           metrics.particle_count_conductor > 0 && metrics.particle_count_air > 0;
}

inline void printSourceDrivenEmSolveMetrics(const char *test_name, const AphiSourceDrivenEmSolveMetrics &metrics,
                                             bool passed)
{
    std::cout << test_name << " passed=" << (passed ? 1 : 0) << " converged=" << (metrics.converged ? 1 : 0)
              << " num_iterations=" << metrics.num_iterations << " final_residual=" << metrics.final_residual
              << " max_abs_A=" << metrics.max_abs_A << " max_abs_phi=" << metrics.max_abs_phi
              << " max_abs_B=" << metrics.max_abs_B << " max_abs_H=" << metrics.max_abs_H
              << " max_abs_E=" << metrics.max_abs_E << " max_abs_J=" << metrics.max_abs_J
              << " max_abs_J_conductor=" << metrics.max_abs_J_conductor << " max_Joule=" << metrics.max_Joule
              << " max_Joule_conductor=" << metrics.max_Joule_conductor
              << " air_Joule_integral=" << metrics.air_Joule_integral
              << " conductor_Joule_integral=" << metrics.conductor_Joule_integral
              << " source_Joule_integral=" << metrics.source_Joule_integral
              << " physical_box_Joule_integral=" << metrics.physical_box_Joule_integral
              << " total_Joule_integral=" << metrics.total_Joule_integral
              << " source_rhs_l2_source_region=" << metrics.source_rhs_l2_source_region
              << " shell_Joule_integral=" << metrics.shell_joule_integral
              << " finite_field_check=" << (metrics.finite_field_check ? 1 : 0) << " particles=" << metrics.particles
              << " vtp_output_dir=" << metrics.vtp_output_dir << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_SOURCE_DRIVEN_EM_SOLVE_HELPERS_H
