#ifndef APHI_TEAM7_CONTACT_TEST_HELPERS_H
#define APHI_TEAM7_CONTACT_TEST_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_assemble_lhs_debug_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline Vecd boxRegionCenter(const benchmark::AphiBoxRegion &region)
{
    return Vecd(0.5 * (region.xmin + region.xmax), 0.5 * (region.ymin + region.ymax),
                0.5 * (region.zmin + region.zmax));
}

inline Vecd boxRegionHalfSize(const benchmark::AphiBoxRegion &region)
{
    return Vecd(0.5 * (region.xmax - region.xmin), 0.5 * (region.ymax - region.ymin),
                0.5 * (region.zmax - region.zmin));
}

/** Non-overlapping TEAM7-like air slabs: left, coil-conductor gap, right. */
class AphiTeam7AirSlabsShape : public ComplexShape
{
  public:
    AphiTeam7AirSlabsShape(const std::string &shape_name, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                           Real body_length, Real body_height, Real body_width)
        : ComplexShape(shape_name)
    {
        const benchmark::AphiBoxRegion left_air{0.0, layout.coil.xmin, 0.0, body_height, 0.0, body_width};
        const benchmark::AphiBoxRegion gap_air{layout.coil.xmax, layout.conductor.xmin, 0.0, body_height, 0.0,
                                               body_width};
        const benchmark::AphiBoxRegion right_air{layout.conductor.xmax, body_length, 0.0, body_height, 0.0, body_width};
        for (const benchmark::AphiBoxRegion &region : {left_air, gap_air, right_air})
        {
            add<GeometricShapeBox>(Transform(boxRegionCenter(region)), boxRegionHalfSize(region));
        }
    }
};

class AphiTeam7RegionBoxShape : public ComplexShape
{
  public:
    AphiTeam7RegionBoxShape(const std::string &shape_name, const benchmark::AphiBoxRegion &region)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(boxRegionCenter(region)), boxRegionHalfSize(region));
    }
};

inline bool isAwayFromTeam7ContactInterfaces(const Vecd &position, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                             Real dp_0)
{
    const Real margin = dp_0 + TinyReal;
    const Real interfaces[] = {layout.coil.xmin, layout.coil.xmax, layout.conductor.xmin, layout.conductor.xmax};
    for (const Real x_interface : interfaces)
    {
        if (std::abs(position[0] - x_interface) <= margin)
        {
            return false;
        }
    }
    return true;
}

inline void collectTeam7CoreBlockByPosition(AphiBlockMapByPosition &block_map, BaseParticles &particles,
                                            const AphiBlockNames &block, Real body_length, Real body_height,
                                            Real body_width, Real core_shell,
                                            const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real dp_0)
{
    syncAphiBlockToHost(particles, block);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block.phi_imag);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!isAwayFromTeam7ContactInterfaces(positions[i], layout, dp_0))
        {
            continue;
        }

        AphiContactBlockByPosition value;
        value.a_real = a_real[i];
        value.a_imag = a_imag[i];
        value.phi_real = phi_real[i];
        value.phi_imag = phi_imag[i];
        block_map[makeContactPositionKey(positions[i])] = value;
    }
}

inline void collectTeam7CoreScalarByPosition(std::map<AphiContactPositionKey, Real> &value_map, BaseParticles &particles,
                                             const std::string &variable_name, Real body_length, Real body_height,
                                             Real body_width, Real core_shell,
                                             const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real dp_0)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!isAwayFromTeam7ContactInterfaces(positions[i], layout, dp_0))
        {
            continue;
        }
        value_map[makeContactPositionKey(positions[i])] = values[i];
    }
}

struct AphiTeam7ThreeBodyContactCase
{
    SPHSystem sph_system;
    benchmark::AphiTeam7LikeUnitBoxLayout layout;
    SolidBody air_body;
    SolidBody coil_body;
    SolidBody plate_body;
    UniquePtr<Inner<>> air_inner_ck;
    UniquePtr<Inner<>> coil_inner_ck;
    UniquePtr<Inner<>> plate_inner_ck;
    UniquePtr<Contact<>> air_contact_ck;
    UniquePtr<Contact<>> coil_to_air_ck;
    UniquePtr<Contact<>> plate_to_air_ck;

    AphiTeam7ThreeBodyContactCase(Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
                                  int ac, char *av[])
        : sph_system(BoundingBoxd(Vecd(-boundary_width, -boundary_width, -boundary_width),
                                  Vecd(body_length + boundary_width, body_height + boundary_width,
                                       body_width + boundary_width)),
                     dp_0),
          layout(benchmark::buildTeam7LayoutForBox(body_length, body_height, body_width)),
          air_body(sph_system,
                   makeShared<AphiTeam7AirSlabsShape>("AirBody", layout, body_length, body_height, body_width)),
          coil_body(sph_system, makeShared<AphiTeam7RegionBoxShape>("CoilBody", layout.coil)),
          plate_body(sph_system, makeShared<AphiTeam7RegionBoxShape>("PlateBody", layout.conductor))
    {
        if (ac > 0)
        {
            sph_system.handleCommandlineOptions(ac, av);
        }
        for (auto *body_ptr : {&air_body, &coil_body, &plate_body})
        {
            body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
            body_ptr->defineMaterial<Solid>();
            body_ptr->defineBodyLevelSetShape();
            body_ptr->generateParticles<BaseParticles, Lattice>();
        }
        air_inner_ck = makeUnique<Inner<>>(air_body);
        coil_inner_ck = makeUnique<Inner<>>(coil_body);
        plate_inner_ck = makeUnique<Inner<>>(plate_body);
        air_contact_ck = makeUnique<Contact<>>(air_body, StdVec<RealBody *>{&coil_body, &plate_body});
        coil_to_air_ck = makeUnique<Contact<>>(coil_body, StdVec<RealBody *>{&air_body});
        plate_to_air_ck = makeUnique<Contact<>>(plate_body, StdVec<RealBody *>{&air_body});
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
    }

    Inner<> &air_inner() { return *air_inner_ck; }
    Inner<> &coil_inner() { return *coil_inner_ck; }
    Inner<> &plate_inner() { return *plate_inner_ck; }
    Contact<> &air_contact() { return *air_contact_ck; }
    Contact<> &coil_to_air() { return *coil_to_air_ck; }
    Contact<> &plate_to_air() { return *plate_to_air_ck; }

    void updateRelations()
    {
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_air_cell_linked_list(air_body);
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_coil_cell_linked_list(coil_body);
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_plate_cell_linked_list(plate_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_air_inner(air_inner());
        UpdateRelation<MainExecutionPolicy, Inner<>> update_coil_inner(coil_inner());
        UpdateRelation<MainExecutionPolicy, Inner<>> update_plate_inner(plate_inner());
        UpdateRelation<MainExecutionPolicy, Contact<>> update_air_contact(air_contact());
        UpdateRelation<MainExecutionPolicy, Contact<>> update_coil_to_air(coil_to_air());
        UpdateRelation<MainExecutionPolicy, Contact<>> update_plate_to_air(plate_to_air());
        update_air_cell_linked_list.exec();
        update_coil_cell_linked_list.exec();
        update_plate_cell_linked_list.exec();
        update_air_inner.exec();
        update_coil_inner.exec();
        update_plate_inner.exec();
        update_air_contact.exec();
        update_coil_to_air.exec();
        update_plate_to_air.exec();
    }
};

inline void initializeTeam7ThreeBodyFields(AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
                                           bool source_only_coil = false)
{
    const Real coil_sigma =
        source_only_coil ? 0.0 : case_setup.layout.coil_material.sigma;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_air(
        case_setup.air_body, case_setup.layout.air.sigma, case_setup.layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_coil(
        case_setup.coil_body, coil_sigma, case_setup.layout.coil_material.nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_plate(
        case_setup.plate_body, case_setup.layout.conductor_material.sigma, case_setup.layout.conductor_material.nu,
        names);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_air_sigma(
        case_setup.air_body, case_setup.layout.air.sigma, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_coil_sigma(
        case_setup.coil_body, coil_sigma, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_plate_sigma(
        case_setup.plate_body, case_setup.layout.conductor_material.sigma, names.material);

    initialize_air.exec();
    initialize_coil.exec();
    initialize_plate.exec();
    assign_air_sigma.exec();
    assign_coil_sigma.exec();
    assign_plate_sigma.exec();
}

inline Real hostRegionLhsBlockMax(BaseParticles &particles, const AphiBlockNames &block, const Vecd *positions,
                                  size_t total_real_particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                  AphiBenchmarkMaterialRegion region, Real body_length, Real body_height, Real body_width,
                                  Real core_shell, Real dp_0)
{
    syncAphiBlockToHost(particles, block);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block.phi_imag);
    Real max_norm = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!team7ParticleInRegion(positions[i], layout, region))
        {
            continue;
        }
        max_norm = std::max(max_norm, a_real[i].norm());
        max_norm = std::max(max_norm, a_imag[i].norm());
        max_norm = std::max(max_norm, std::abs(phi_real[i]));
        max_norm = std::max(max_norm, std::abs(phi_imag[i]));
    }
    return max_norm;
}

inline Real hostRegionJouleMax(BaseParticles &particles, const std::string &variable_name, const Vecd *positions,
                               size_t total_real_particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                               AphiBenchmarkMaterialRegion region, Real body_length, Real body_height, Real body_width,
                               Real core_shell, Real dp_0)
{
    syncVariableToHost<Real>(particles, variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!team7ParticleInRegion(positions[i], layout, region))
        {
            continue;
        }
        max_value = std::max(max_value, particles.getVariableDataByName<Real>(variable_name)[i]);
    }
    return max_value;
}

struct AphiTeam7RegionalApplyMetrics
{
    Real mono_conductor_lhs_max = 0.0;
    Real contact_conductor_lhs_max = 0.0;
    Real mono_coil_lhs_max = 0.0;
    Real contact_coil_lhs_max = 0.0;
    Real mono_air_lhs_max = 0.0;
    Real contact_air_lhs_max = 0.0;
    Real mono_conductor_joule_max = 0.0;
    Real contact_conductor_joule_max = 0.0;
};

inline AphiTeam7RegionalApplyMetrics hostTeam7RegionalApplyMetricsFromMonolithic(
    SPHBody &body, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real dp_0)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AphiTeam7RegionalApplyMetrics metrics;
    metrics.mono_conductor_lhs_max =
        hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, layout,
                              AphiBenchmarkMaterialRegion::Conductor, body_length, body_height, body_width, core_shell,
                              dp_0);
    metrics.mono_coil_lhs_max = hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, layout,
                                                      AphiBenchmarkMaterialRegion::Coil, body_length, body_height,
                                                      body_width, core_shell, dp_0);
    metrics.mono_air_lhs_max = hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, layout,
                                                     AphiBenchmarkMaterialRegion::Air, body_length, body_height,
                                                     body_width, core_shell, dp_0);
    metrics.mono_conductor_joule_max =
        hostRegionJouleMax(particles, joule_names.joule_heat_source, positions, total_real_particles, layout,
                           AphiBenchmarkMaterialRegion::Conductor, body_length, body_height, body_width, core_shell,
                           dp_0);
    return metrics;
}

inline AphiTeam7RegionalApplyMetrics hostTeam7RegionalApplyMetricsFromContactCase(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
    const AphiJouleHeatingFieldNames &joule_names, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real dp_0)
{
    AphiTeam7RegionalApplyMetrics metrics;
    for (auto *body_ptr : {&case_setup.air_body, &case_setup.coil_body, &case_setup.plate_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        metrics.contact_conductor_lhs_max =
            std::max(metrics.contact_conductor_lhs_max,
                     hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, case_setup.layout,
                                           AphiBenchmarkMaterialRegion::Conductor, body_length, body_height,
                                           body_width, core_shell, dp_0));
        metrics.contact_coil_lhs_max =
            std::max(metrics.contact_coil_lhs_max,
                     hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, case_setup.layout,
                                           AphiBenchmarkMaterialRegion::Coil, body_length, body_height, body_width,
                                           core_shell, dp_0));
        metrics.contact_air_lhs_max =
            std::max(metrics.contact_air_lhs_max,
                     hostRegionLhsBlockMax(particles, names.lhs, positions, total_real_particles, case_setup.layout,
                                           AphiBenchmarkMaterialRegion::Air, body_length, body_height, body_width,
                                           core_shell, dp_0));
        metrics.contact_conductor_joule_max =
            std::max(metrics.contact_conductor_joule_max,
                     hostRegionJouleMax(particles, joule_names.joule_heat_source, positions, total_real_particles,
                                        case_setup.layout, AphiBenchmarkMaterialRegion::Conductor, body_length,
                                        body_height, body_width, core_shell, dp_0));
    }
    return metrics;
}

struct AphiTeam7PlateObservables
{
    Real plate_joule_power = 0.0;
    Real plate_joule_max = 0.0;
    Real plate_Joule_L2 = 0.0;
    Real plate_j_L2 = 0.0;
    Real plate_E_L2 = 0.0;
    Real plate_E_real_max = 0.0;
    Real plate_E_imag_max = 0.0;
    Real plate_J_real_max = 0.0;
    Real plate_J_imag_max = 0.0;
    Real plate_solution_block_max = 0.0;
};

inline Real hostRegionVolWeightedVecdL2(BaseParticles &particles, const std::string &variable_name,
                                       const Vecd *positions, size_t total_real_particles,
                                       const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                       AphiBenchmarkMaterialRegion region)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!team7ParticleInRegion(positions[i], layout, region))
        {
            continue;
        }
        sum_squared += vol[i] * values[i].squaredNorm();
    }
    return std::sqrt(sum_squared);
}

inline Real hostRegionScalarMax(BaseParticles &particles, const std::string &variable_name, const Vecd *positions,
                                size_t total_real_particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                AphiBenchmarkMaterialRegion region)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!team7ParticleInRegion(positions[i], layout, region))
        {
            continue;
        }
        max_value = std::max(max_value, values[i]);
    }
    return max_value;
}

inline AphiTeam7PlateObservables hostTeam7PlateObservablesFromBody(
    SPHBody &body, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real body_length, Real body_height, Real body_width,
    Real core_shell)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    const AphiRegionalElectromagneticObservables region_observables = hostRegionalElectromagneticObservablesFromTeam7Region(
        particles, names, joule_names, positions, total_real_particles, layout, AphiBenchmarkMaterialRegion::Conductor,
        body_length, body_height, body_width, core_shell);

    AphiTeam7PlateObservables observables;
    observables.plate_joule_power = region_observables.Joule_power;
    observables.plate_joule_max = region_observables.Joule_max;
    observables.plate_Joule_L2 = region_observables.Joule_L2;
    observables.plate_E_L2 = region_observables.E_L2;
    observables.plate_E_real_max = region_observables.E_real_max;
    observables.plate_E_imag_max = region_observables.E_imag_max;
    observables.plate_J_real_max = region_observables.J_real_max;
    observables.plate_J_imag_max = region_observables.J_imag_max;
    observables.plate_j_L2 = region_observables.J_L2;
    observables.plate_solution_block_max = region_observables.solution_block_max;
    return observables;
}

inline AphiTeam7PlateObservables hostTeam7PlateObservablesFromContactCase(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
    const AphiJouleHeatingFieldNames &joule_names, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    return hostTeam7PlateObservablesFromBody(case_setup.plate_body, names, joule_names, case_setup.layout, body_length,
                                            body_height, body_width, core_shell);
}

inline void runTeam7ContactJoulePipeline(AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
                                         const AphiJouleHeatingFieldNames &joule_names, Real omega)
{
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> air_grad_phi(
        DynamicsArgs(case_setup.air_inner(), names.solution, joule_names),
        DynamicsArgs(case_setup.air_contact(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> coil_grad_phi(
        DynamicsArgs(case_setup.coil_inner(), names.solution, joule_names),
        DynamicsArgs(case_setup.coil_to_air(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> plate_grad_phi(
        DynamicsArgs(case_setup.plate_inner(), names.solution, joule_names),
        DynamicsArgs(case_setup.plate_to_air(), names.solution, joule_names));

    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> air_electric_field(
        case_setup.air_body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> coil_electric_field(
        case_setup.coil_body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> plate_electric_field(
        case_setup.plate_body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> air_joule(case_setup.air_body, names.material,
                                                                               joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> coil_joule(case_setup.coil_body, names.material,
                                                                                joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> plate_joule(case_setup.plate_body, names.material,
                                                                                 joule_names);

    air_grad_phi.exec();
    coil_grad_phi.exec();
    plate_grad_phi.exec();
    air_electric_field.exec();
    coil_electric_field.exec();
    plate_electric_field.exec();
    air_joule.exec();
    coil_joule.exec();
    plate_joule.exec();
}

struct AphiTeam7ContactBlockGsResult
{
    AphiMatrixFreeSolverResult last_solver_result{};
    Real global_true_rel = 0.0;
    UnsignedInt block_gs_sweeps = 0;
    AphiTeam7PlateObservables plate{};
};

inline Real threeBodyGlobalTrueRelativeResidual(AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names)
{
    Real sum_squared = 0.0;
    Real rhs_sum_squared = 0.0;
    for (auto *body_ptr : {&case_setup.air_body, &case_setup.coil_body, &case_setup.plate_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        sum_squared += std::pow(hostBlockNorm(particles, names.true_residual, total_real_particles), 2);
        rhs_sum_squared += std::pow(hostBlockNorm(particles, names.rhs, total_real_particles), 2);
    }
    return std::sqrt(sum_squared) / (std::sqrt(rhs_sum_squared) + TinyReal);
}

inline void computeThreeBodyTrueResiduals(AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
                                          const AphiLhsAssemblyOptions &options)
{
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_air(
        case_setup.air_body, case_setup.air_inner(), case_setup.air_contact(), names.solution, names.lhs, names.material,
        options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_coil(
        case_setup.coil_body, case_setup.coil_inner(), case_setup.coil_to_air(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_plate(
        case_setup.plate_body, case_setup.plate_inner(), case_setup.plate_to_air(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_air_true_residual(
        case_setup.air_body, names.true_residual, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_coil_true_residual(
        case_setup.coil_body, names.true_residual, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_plate_true_residual(
        case_setup.plate_body, names.true_residual, names.rhs, names.lhs);

    apply_air.exec();
    apply_coil.exec();
    apply_plate.exec();
    compute_air_true_residual.exec();
    compute_coil_true_residual.exec();
    compute_plate_true_residual.exec();
}

struct AphiTeam7ContactCoupledGmresResult
{
    AphiGMRESResult solver_result{};
    Real global_true_rel = 0.0;
    Real air_true_rel = 0.0;
    Real coil_true_rel = 0.0;
    Real plate_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    Real air_rhs_norm = 0.0;
    Real coil_rhs_norm = 0.0;
    Real plate_rhs_norm = 0.0;
    Real air_solution_norm = 0.0;
    Real coil_solution_norm = 0.0;
    Real plate_solution_norm = 0.0;
    AphiTeam7PlateObservables plate{};
};

inline StdVec<AphiMultiBodyContactEntry> buildTeam7MultiBodyEntries(AphiTeam7ThreeBodyContactCase &case_setup)
{
    return StdVec<AphiMultiBodyContactEntry>{
        AphiMultiBodyContactEntry{case_setup.air_body, &case_setup.air_contact(), &case_setup.air_inner()},
        AphiMultiBodyContactEntry{case_setup.coil_body, &case_setup.coil_to_air(), &case_setup.coil_inner()},
        AphiMultiBodyContactEntry{case_setup.plate_body, &case_setup.plate_to_air(), &case_setup.plate_inner()}};
}

inline AphiTeam7ContactCoupledGmresResult runTeam7ThreeBodyCoupledContactGmres(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const AphiLhsAssemblyOptions &options, const AphiMatrixFreeSolverOptions &solver_options, Real body_length,
    Real body_height, Real body_width, Real core_shell, const Vecd &coil_current_real, const Vecd &coil_current_imag,
    Real impressed_current_amplitude, bool source_only_coil = false)
{
    AphiTeam7ContactCoupledGmresResult result;

    initializeTeam7ThreeBodyFields(case_setup, names, source_only_coil);
    const AphiGMRESWorkspaceNames gmres_workspace =
        buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension);
    for (auto *body_ptr : {&case_setup.air_body, &case_setup.coil_body, &case_setup.plate_body})
    {
        RegisterAphiGMRESWorkspaceCK register_workspace(*body_ptr, solver_options.gmres.restart_dimension);
        (void)register_workspace;
    }
    RegisterAphiJouleHeatingFieldsCK register_air_joule(case_setup.air_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_coil_joule(case_setup.coil_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_plate_joule(case_setup.plate_body, joule_names);
    (void)register_air_joule;
    (void)register_coil_joule;
    (void)register_plate_joule;

    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_rhs(case_setup.air_body, names.rhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_rhs(case_setup.plate_body, names.rhs);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_coil_source(
        case_setup.coil_body, names.rhs, case_setup.layout.coil, coil_current_real, coil_current_imag,
        impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_solution(case_setup.air_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_coil_solution(case_setup.coil_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_solution(case_setup.plate_body, names.solution);

    zero_air_rhs.exec();
    zero_plate_rhs.exec();
    assign_coil_source.exec();
    zero_air_solution.exec();
    zero_coil_solution.exec();
    zero_plate_solution.exec();
    case_setup.updateRelations();

    const StdVec<AphiMultiBodyContactEntry> entries = buildTeam7MultiBodyEntries(case_setup);
    AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy> solver(entries, names, gmres_workspace, options,
                                                                  solver_options.gmres);
    result.solver_result = solver.solve();
    case_setup.updateRelations();

    computeThreeBodyTrueResiduals(case_setup, names, options);
    const AphiBodywiseTrueRelativeBreakdown bodywise = buildBodywiseTrueRelativeBreakdown(
        {&case_setup.air_body, &case_setup.coil_body, &case_setup.plate_body}, names);
    result.global_true_rel = bodywise.global_true_rel;
    result.max_bodywise_true_rel = bodywise.max_bodywise_true_rel;
    result.air_true_rel = bodywise.bodies[0].true_rel;
    result.coil_true_rel = bodywise.bodies[1].true_rel;
    result.plate_true_rel = bodywise.bodies[2].true_rel;
    result.air_rhs_norm = bodywise.bodies[0].rhs_norm;
    result.coil_rhs_norm = bodywise.bodies[1].rhs_norm;
    result.plate_rhs_norm = bodywise.bodies[2].rhs_norm;
    result.air_solution_norm = bodywise.bodies[0].solution_norm;
    result.coil_solution_norm = bodywise.bodies[1].solution_norm;
    result.plate_solution_norm = bodywise.bodies[2].solution_norm;

    runTeam7ContactJoulePipeline(case_setup, names, joule_names, options.omega);
    result.plate = hostTeam7PlateObservablesFromContactCase(case_setup, names, joule_names, body_length, body_height,
                                                            body_width, core_shell);
    return result;
}

inline AphiTeam7ContactBlockGsResult runTeam7ThreeBodyContactBlockGs(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const AphiLhsAssemblyOptions &options, const AphiMatrixFreeSolverOptions &solver_options, Real body_length,
    Real body_height, Real body_width, Real core_shell, const Vecd &coil_current_real, const Vecd &coil_current_imag,
    Real impressed_current_amplitude, UnsignedInt max_block_gs_sweeps)
{
    AphiTeam7ContactBlockGsResult result;

    initializeTeam7ThreeBodyFields(case_setup, names);
    RegisterAphiGMRESWorkspaceCK register_air_workspace(case_setup.air_body, solver_options.gmres.restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_coil_workspace(case_setup.coil_body, solver_options.gmres.restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_plate_workspace(case_setup.plate_body, solver_options.gmres.restart_dimension);
    RegisterAphiJouleHeatingFieldsCK register_air_joule(case_setup.air_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_coil_joule(case_setup.coil_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_plate_joule(case_setup.plate_body, joule_names);
    (void)register_air_workspace;
    (void)register_coil_workspace;
    (void)register_plate_workspace;
    (void)register_air_joule;
    (void)register_coil_joule;
    (void)register_plate_joule;

    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_rhs(case_setup.air_body, names.rhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_rhs(case_setup.plate_body, names.rhs);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_coil_source(
        case_setup.coil_body, names.rhs, case_setup.layout.coil, coil_current_real, coil_current_imag,
        impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_solution(case_setup.air_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_coil_solution(case_setup.coil_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_solution(case_setup.plate_body, names.solution);

    zero_air_rhs.exec();
    zero_plate_rhs.exec();
    assign_coil_source.exec();
    zero_air_solution.exec();
    zero_coil_solution.exec();
    zero_plate_solution.exec();
    case_setup.updateRelations();

    const AphiGMRESWorkspaceNames gmres_workspace =
        buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension);

    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> air_solver(
        case_setup.air_body, case_setup.air_inner(), case_setup.air_contact(), names, options, solver_options);
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> coil_solver(
        case_setup.coil_body, case_setup.coil_inner(), case_setup.coil_to_air(), names, options, solver_options);
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> plate_solver(
        case_setup.plate_body, case_setup.plate_inner(), case_setup.plate_to_air(), names, options, solver_options);

    for (UnsignedInt sweep = 0; sweep < max_block_gs_sweeps; ++sweep)
    {
        zeroGmresWorkspaceOnBody(case_setup.air_body, gmres_workspace);
        zeroGmresWorkspaceOnBody(case_setup.coil_body, gmres_workspace);
        zeroGmresWorkspaceOnBody(case_setup.plate_body, gmres_workspace);

        result.last_solver_result = coil_solver.solve();
        case_setup.updateRelations();
        result.last_solver_result = air_solver.solve();
        case_setup.updateRelations();
        result.last_solver_result = plate_solver.solve();
        case_setup.updateRelations();

        computeThreeBodyTrueResiduals(case_setup, names, options);
        result.global_true_rel = threeBodyGlobalTrueRelativeResidual(case_setup, names);
        result.block_gs_sweeps = sweep + 1;
    }

    runTeam7ContactJoulePipeline(case_setup, names, joule_names, options.omega);
    result.plate = hostTeam7PlateObservablesFromContactCase(case_setup, names, joule_names, body_length, body_height,
                                                            body_width, core_shell);
    return result;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_CONTACT_TEST_HELPERS_H
