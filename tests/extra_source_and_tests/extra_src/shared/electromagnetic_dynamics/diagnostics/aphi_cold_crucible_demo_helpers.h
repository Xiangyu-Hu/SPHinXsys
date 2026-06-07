#ifndef APHI_COLD_CRUCIBLE_DEMO_HELPERS_H
#define APHI_COLD_CRUCIBLE_DEMO_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_cold_crucible_case_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_h_curl_h_j_observable_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include "io_vtk.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Pipeline smoke-test fields (not production TEAM7 validation). */
struct AphiColdCrucibleDemoFieldNames
{
    std::string b_real = "BReal";
    std::string b_imag = "BImag";
    std::string b_magnitude = "BMagnitude";
    std::string h_real = "HReal";
    std::string h_imag = "HImag";
    std::string h_magnitude = "HMagnitude";
    std::string a_magnitude = "AMagnitude";
    std::string e_magnitude = "EMagnitude";
    std::string j_magnitude = "JMagnitude";
    std::string material_region_id = "MaterialRegionId";
};

struct AphiColdCrucibleDemoSpec
{
    static constexpr Real body_length = benchmark::AphiColdCruciblePhysicalDimensions::length;
    static constexpr Real body_height = benchmark::AphiColdCruciblePhysicalDimensions::height;
    static constexpr Real body_width = benchmark::AphiColdCruciblePhysicalDimensions::width;

    Real dp = 0.1;
    Real omega = 100.0;
    Real impressed_a0 = 0.05;
    /** Coil ring radius as fraction of body_length (centered in x-y). */
    Real impressed_r_coil = 0.50;
    Real impressed_w_r = 0.12;
    Real impressed_w_z = 0.20;
    Real sigma_air = 1.0e-8;
    Real sigma_melt = 1.0;
    Real sigma_wall = 0.1;
    Real nu_uniform = 1.0;
    Real max_air_j_norm = 1.0e-4;
    Real max_air_joule_integral = 1.0e-8;
    Real min_melt_j_norm = 1.0e-10;
    Real min_melt_joule_integral = 1.0e-14;
    bool write_vtp = true;
};

struct AphiColdCrucibleDemoMetrics
{
    Real max_a = 0.0;
    Real max_b = 0.0;
    Real max_e = 0.0;
    Real max_j = 0.0;
    Real max_joule = 0.0;
    Real melt_joule_integral = 0.0;
    Real air_joule_integral = 0.0;
    Real melt_j_max = 0.0;
    Real air_j_max = 0.0;
    size_t particles = 0;
    std::string vtp_output_dir = "output";
};

inline Vecd impressedCoilVectorPotentialReal(const Vecd &position, Real a0, Real r_coil, Real w_r, Real w_z,
                                             Real z_center)
{
    const Real x = position[0];
    const Real y = position[1];
    const Real z = position[2];
    const Real r = std::sqrt(x * x + y * y) + TinyReal;
    const Real dz = z - z_center;
    const Real exponent = -((r - r_coil) * (r - r_coil) / (w_r * w_r + TinyReal) + dz * dz / (w_z * w_z + TinyReal));
    const Real a_theta = a0 * std::exp(exponent);
    return Vecd(-a_theta * y / r, a_theta * x / r, 0.0);
}

inline benchmark::AphiColdCrucibleUnitBoxLayout buildColdCrucibleDemoLayout(const AphiColdCrucibleDemoSpec &spec)
{
    benchmark::AphiColdCrucibleUnitBoxLayout layout =
        benchmark::buildColdCrucibleLayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    layout.air.sigma = spec.sigma_air;
    layout.air.nu = spec.nu_uniform;
    layout.melt_material.sigma = spec.sigma_melt;
    layout.melt_material.nu = spec.nu_uniform;
    layout.crucible_wall_material.sigma = spec.sigma_wall;
    layout.crucible_wall_material.nu = spec.nu_uniform;
    layout.coil_material.sigma = spec.sigma_air;
    layout.coil_material.nu = spec.nu_uniform;
    return layout;
}

class AssignImpressedCoilVectorPotentialCK : public LocalDynamics
{
  public:
    AssignImpressedCoilVectorPotentialCK(SPHBody &sph_body, const AphiBlockNames &solution_block, Real a0, Real r_coil,
                                           Real w_r, Real w_z, Real z_center)
        : LocalDynamics(sph_body), a0_(a0), r_coil_(r_coil), w_r_(w_r), w_z_(w_z), z_center_(z_center),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(solution_block.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(solution_block.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(solution_block.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(solution_block.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : a0_(encloser.a0_), r_coil_(encloser.r_coil_), w_r_(encloser.w_r_), w_z_(encloser.w_z_),
              z_center_(encloser.z_center_), position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            a_real_[index_i] = impressedCoilVectorPotentialReal(position, a0_, r_coil_, w_r_, w_z_, z_center_);
            a_imag_[index_i] = Vecd::Zero();
            phi_real_[index_i] = 0.0;
            phi_imag_[index_i] = 0.0;
        }

      protected:
        Real a0_;
        Real r_coil_;
        Real w_r_;
        Real w_z_;
        Real z_center_;
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    Real a0_;
    Real r_coil_;
    Real w_r_;
    Real w_z_;
    Real z_center_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

class AssignColdCrucibleDemoRegionIdCK : public LocalDynamics
{
  public:
    AssignColdCrucibleDemoRegionIdCK(SPHBody &sph_body, const benchmark::AphiColdCrucibleUnitBoxLayout &layout,
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
            if (benchmark::insideBoxRegion(position, layout_.melt))
            {
                region_id_[index_i] = 2.0;
                return;
            }
            if (benchmark::insideBoxRegion(position, layout_.crucible_wall))
            {
                region_id_[index_i] = 1.0;
                return;
            }
            if (benchmark::insideAnyCoilRegion(position, layout_))
            {
                region_id_[index_i] = 3.0;
                return;
            }
            region_id_[index_i] = 0.0;
        }

      protected:
        Vecd *position_;
        Real *region_id_;
        benchmark::AphiColdCrucibleUnitBoxLayout layout_;
    };

  protected:
    benchmark::AphiColdCrucibleUnitBoxLayout layout_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_region_id_;
};

class AphiComplexVecdMagnitudeCK : public LocalDynamics
{
  public:
    AphiComplexVecdMagnitudeCK(SPHBody &sph_body, const std::string &real_name, const std::string &imag_name,
                               const std::string &magnitude_name)
        : LocalDynamics(sph_body),
          dv_real_(particles_->template getVariableByName<Vecd>(real_name)),
          dv_imag_(particles_->template getVariableByName<Vecd>(imag_name)),
          dv_magnitude_(particles_->template registerStateVariable<Real>(magnitude_name, Real(0)))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : real_(encloser.dv_real_->DelegatedData(ex_policy)),
              imag_(encloser.dv_imag_->DelegatedData(ex_policy)),
              magnitude_(encloser.dv_magnitude_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            magnitude_[index_i] = std::sqrt(real_[index_i].squaredNorm() + imag_[index_i].squaredNorm());
        }

      protected:
        Vecd *real_;
        Vecd *imag_;
        Real *magnitude_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_real_;
    DiscreteVariable<Vecd> *dv_imag_;
    DiscreteVariable<Real> *dv_magnitude_;
};

class AphiScaleVecdFieldByNuCK : public LocalDynamics
{
  public:
    AphiScaleVecdFieldByNuCK(SPHBody &sph_body, const std::string &input_name, const std::string &output_name,
                             const AphiMaterialNames &material_names)
        : LocalDynamics(sph_body), dv_input_(particles_->template getVariableByName<Vecd>(input_name)),
          dv_output_(particles_->template getVariableByName<Vecd>(output_name)),
          dv_nu_(particles_->template getVariableByName<Real>(material_names.nu))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : input_(encloser.dv_input_->DelegatedData(ex_policy)),
              output_(encloser.dv_output_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            output_[index_i] = nu_[index_i] * input_[index_i];
        }

      protected:
        Vecd *input_;
        Vecd *output_;
        Real *nu_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Vecd> *dv_output_;
    DiscreteVariable<Real> *dv_nu_;
};

inline void registerColdCrucibleDemoObservableFields(SPHBody &body, const AphiColdCrucibleDemoFieldNames &demo_fields)
{
    BaseParticles &particles = body.getBaseParticles();
    particles.registerStateVariable<Vecd>(demo_fields.b_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(demo_fields.b_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(demo_fields.h_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(demo_fields.h_imag, ZeroData<Vecd>::value);
}

inline void syncColdCrucibleDemoFieldsToHost(SPHBody &body, const AphiVariableNames &names,
                                             const AphiJouleHeatingFieldNames &joule_fields,
                                             const AphiColdCrucibleDemoFieldNames &demo_fields)
{
    BaseParticles &particles = body.getBaseParticles();
    syncAphiBlockToHost(particles, names.solution);
    syncAphiMaterialToHost(particles, names.material);
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_real);
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_imag);
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_real);
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_imag);
    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToHost<Vecd>(particles, demo_fields.b_real);
    syncVariableToHost<Vecd>(particles, demo_fields.b_imag);
    syncVariableToHost<Vecd>(particles, demo_fields.h_real);
    syncVariableToHost<Vecd>(particles, demo_fields.h_imag);
    syncVariableToHost<Real>(particles, demo_fields.a_magnitude);
    syncVariableToHost<Real>(particles, demo_fields.b_magnitude);
    syncVariableToHost<Real>(particles, demo_fields.e_magnitude);
    syncVariableToHost<Real>(particles, demo_fields.j_magnitude);
    syncVariableToHost<Real>(particles, demo_fields.h_magnitude);
    syncVariableToHost<Real>(particles, demo_fields.material_region_id);
}

inline void writeColdCrucibleDemoVtp(SPHSystem &sph_system, SPHBody &body, const AphiVariableNames &names,
                                     const AphiJouleHeatingFieldNames &joule_fields,
                                     const AphiColdCrucibleDemoFieldNames &demo_fields,
                                     const AphiMaterialNames &material_names)
{
    syncColdCrucibleDemoFieldsToHost(body, names, joule_fields, demo_fields);
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(body, names.solution.a_real);
    write_states.addToWrite<Vecd>(body, names.solution.a_imag);
    write_states.addToWrite<Real>(body, names.solution.phi_real);
    write_states.addToWrite<Real>(body, names.solution.phi_imag);
    write_states.addToWrite<Vecd>(body, demo_fields.b_real);
    write_states.addToWrite<Vecd>(body, demo_fields.b_imag);
    write_states.addToWrite<Vecd>(body, demo_fields.h_real);
    write_states.addToWrite<Vecd>(body, demo_fields.h_imag);
    write_states.addToWrite<Vecd>(body, joule_fields.electric_field_a_real);
    write_states.addToWrite<Vecd>(body, joule_fields.electric_field_a_imag);
    write_states.addToWrite<Vecd>(body, joule_fields.current_density_real);
    write_states.addToWrite<Vecd>(body, joule_fields.current_density_imag);
    write_states.addToWrite<Real>(body, joule_fields.joule_heat_source);
    write_states.addToWrite<Real>(body, material_names.sigma);
    write_states.addToWrite<Real>(body, material_names.nu);
    write_states.addToWrite<Real>(body, demo_fields.a_magnitude);
    write_states.addToWrite<Real>(body, demo_fields.b_magnitude);
    write_states.addToWrite<Real>(body, demo_fields.e_magnitude);
    write_states.addToWrite<Real>(body, demo_fields.j_magnitude);
    write_states.addToWrite<Real>(body, demo_fields.h_magnitude);
    write_states.addToWrite<Real>(body, demo_fields.material_region_id);
    write_states.writeToFile(0);
}

inline AphiColdCrucibleDemoMetrics hostColdCrucibleDemoGlobalMaxima(
    BaseParticles &particles, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_fields,
    const AphiColdCrucibleDemoFieldNames &demo_fields, const Vecd *positions, size_t total_real_particles)
{
    const auto all_particles = [](const Vecd &) { return true; };
    AphiColdCrucibleDemoMetrics metrics;
    metrics.max_a = hostParticleRegionVecdMax(particles, names.solution.a_real, positions, total_real_particles,
                                            all_particles);
    metrics.max_b = hostParticleRegionScalarMax(particles, demo_fields.b_magnitude, positions, total_real_particles,
                                                all_particles);
    metrics.max_e = hostParticleRegionScalarMax(particles, demo_fields.e_magnitude, positions, total_real_particles,
                                                all_particles);
    metrics.max_j = hostParticleRegionScalarMax(particles, demo_fields.j_magnitude, positions, total_real_particles,
                                                all_particles);
    metrics.max_joule = hostParticleRegionScalarMax(particles, joule_fields.joule_heat_source, positions,
                                                    total_real_particles, all_particles);
    return metrics;
}

inline bool coldCrucibleDemoMetricsPassed(const AphiColdCrucibleDemoMetrics &metrics,
                                          const AphiColdCrucibleDemoSpec &spec)
{
    return std::isfinite(metrics.max_a) && metrics.max_a > 0.0 && std::isfinite(metrics.max_b) && metrics.max_b > 0.0 &&
           std::isfinite(metrics.max_e) && metrics.max_e > 0.0 && std::isfinite(metrics.max_j) &&
           metrics.max_j > 0.0 && std::isfinite(metrics.max_joule) && metrics.max_joule > 0.0 &&
           metrics.melt_j_max > spec.min_melt_j_norm && metrics.melt_joule_integral > spec.min_melt_joule_integral &&
           metrics.air_joule_integral <= spec.max_air_joule_integral && metrics.air_j_max <= spec.max_air_j_norm;
}

inline AphiColdCrucibleDemoMetrics runColdCrucibleImpressedDemo(int ac, char *av[], const AphiColdCrucibleDemoSpec &spec)
{
    const Real boundary_width = 3.0 * spec.dp;
    const benchmark::AphiColdCrucibleUnitBoxLayout layout = buildColdCrucibleDemoLayout(spec);
    const Real z_center = 0.5 * spec.body_width;
    const Real r_coil = spec.impressed_r_coil * spec.body_length;
    const Real w_r = spec.impressed_w_r * spec.body_length;
    const Real w_z = spec.impressed_w_z * spec.body_width;

    AphiJouleHeatingFieldNames joule_fields;
    AphiColdCrucibleDemoFieldNames demo_fields;
    AphiVariableNames names;

    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    registerColdCrucibleDemoObservableFields(test_body.body, demo_fields);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignColdCrucibleRegionMaterialsCK> assign_material(
        test_body.body, layout, names.material);
    StateDynamics<MainExecutionPolicy, AssignColdCrucibleDemoRegionIdCK> assign_region_id(
        test_body.body, layout, demo_fields.material_region_id);
    StateDynamics<MainExecutionPolicy, AssignImpressedCoilVectorPotentialCK> assign_impressed_a(
        test_body.body, names.solution, spec.impressed_a0, r_coil, w_r, w_z, z_center);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> a_magnitude(
        test_body.body, names.solution.a_real, names.solution.a_imag, demo_fields.a_magnitude);
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> e_magnitude(
        test_body.body, joule_fields.electric_field_a_real, joule_fields.electric_field_a_imag, demo_fields.e_magnitude);
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> j_magnitude(
        test_body.body, joule_fields.current_density_real, joule_fields.current_density_imag, demo_fields.j_magnitude);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_region_id.exec();
    assign_impressed_a.exec();
    test_body.updateRelations();

    execBodyCurlFromVectorFieldDiagnostic(test_body.body, test_body.inner(), names.solution.a_real,
                                          names.solution.a_imag, demo_fields.b_real, demo_fields.b_imag,
                                          AphiBCurlDiagnosticMode::BCorrectedGrad);
    StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByNuCK> compute_h_real(
        test_body.body, demo_fields.b_real, demo_fields.h_real, names.material);
    StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByNuCK> compute_h_imag(
        test_body.body, demo_fields.b_imag, demo_fields.h_imag, names.material);
    compute_h_real.exec();
    compute_h_imag.exec();

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    a_magnitude.exec();
    e_magnitude.exec();
    j_magnitude.exec();
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> b_magnitude(
        test_body.body, demo_fields.b_real, demo_fields.b_imag, demo_fields.b_magnitude);
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> h_magnitude(
        test_body.body, demo_fields.h_real, demo_fields.h_imag, demo_fields.h_magnitude);
    b_magnitude.exec();
    h_magnitude.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiColdCrucibleDemoMetrics metrics =
        hostColdCrucibleDemoGlobalMaxima(particles, names, joule_fields, demo_fields, positions, total_real_particles);
    metrics.particles = total_real_particles;
    const auto in_melt = [&](const Vecd &position) {
        return coldCrucibleParticleInRegion(position, layout, AphiColdCrucibleMaterialRegion::Melt);
    };
    const auto in_air = [&](const Vecd &position) {
        return coldCrucibleParticleInRegion(position, layout, AphiColdCrucibleMaterialRegion::Air);
    };
    metrics.melt_j_max = hostParticleRegionScalarMax(particles, demo_fields.j_magnitude, positions, total_real_particles,
                                                     in_melt);
    metrics.air_j_max = hostParticleRegionScalarMax(particles, demo_fields.j_magnitude, positions, total_real_particles,
                                                    in_air);
    metrics.melt_joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                        in_melt, joule_fields.joule_heat_source);
    metrics.air_joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                         in_air, joule_fields.joule_heat_source);

    if (spec.write_vtp)
    {
        writeColdCrucibleDemoVtp(test_body.sph_system, test_body.body, names, joule_fields, demo_fields,
                                 names.material);
        metrics.vtp_output_dir = "output";
    }
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_COLD_CRUCIBLE_DEMO_HELPERS_H
