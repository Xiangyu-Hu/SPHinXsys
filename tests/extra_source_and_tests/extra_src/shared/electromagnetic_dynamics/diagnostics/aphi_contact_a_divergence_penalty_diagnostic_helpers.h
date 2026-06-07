#ifndef APHI_CONTACT_A_DIVERGENCE_PENALTY_DIAGNOSTIC_HELPERS_H
#define APHI_CONTACT_A_DIVERGENCE_PENALTY_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/aphi_block_zero_ck.hpp"
#include "electromagnetic_dynamics/aphi_contact_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactADivergencePenaltyEquivalenceMetrics
{
    size_t div_a_matched = 0;
    size_t div_a_missing = 0;
    Real div_a_max_abs_diff = 0.0;
    size_t grad_div_a_stencil_safe_matched = 0;
    size_t grad_div_a_stencil_safe_missing = 0;
    Real grad_div_a_stencil_safe_max_abs_diff = 0.0;
    size_t grad_div_a_interface_band_matched = 0;
    Real grad_div_a_interface_band_max_abs_diff = 0.0;
    size_t lhs_a_stencil_safe_matched = 0;
    size_t lhs_a_stencil_safe_missing = 0;
    Real lhs_a_stencil_safe_max_abs_diff = 0.0;
};

inline Real maxAbsVecdMapDifferenceOnStencilSafeCore(const std::map<AphiContactPositionKey, Vecd> &reference,
                                                       const std::map<AphiContactPositionKey, Vecd> &candidate,
                                                       Real x_interface, Real dp_0, size_t &matched_particles,
                                                       size_t &missing_particles, size_t &interface_band_matched,
                                                       Real &interface_band_max_abs_diff)
{
    matched_particles = 0;
    missing_particles = 0;
    interface_band_matched = 0;
    interface_band_max_abs_diff = 0.0;
    Real max_diff = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = candidate.find(entry.first);
        if (it == candidate.end())
        {
            missing_particles += 1;
            continue;
        }
        const Real diff = (entry.second - it->second).norm();
        if (isGradStencilSafeFromInterface(entry.first.position, x_interface, dp_0))
        {
            matched_particles += 1;
            max_diff = std::max(max_diff, diff);
            continue;
        }
        interface_band_matched += 1;
        interface_band_max_abs_diff = std::max(interface_band_max_abs_diff, diff);
    }

    return max_diff;
}

inline Real maxAbsBlockDifferenceOnStencilSafeCore(const AphiBlockMapByPosition &reference,
                                                   const AphiBlockMapByPosition &candidate, Real x_interface, Real dp_0,
                                                   size_t &matched_particles, size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real max_diff = 0.0;

    for (const auto &entry : reference)
    {
        if (!isGradStencilSafeFromInterface(entry.first.position, x_interface, dp_0))
        {
            continue;
        }
        const auto it = candidate.find(entry.first);
        if (it == candidate.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        const AphiContactBlockByPosition &lhs = entry.second;
        const AphiContactBlockByPosition &rhs = it->second;
        max_diff = std::max(max_diff, (lhs.a_real - rhs.a_real).norm());
        max_diff = std::max(max_diff, (lhs.a_imag - rhs.a_imag).norm());
        max_diff = std::max(max_diff, std::abs(lhs.phi_real - rhs.phi_real));
        max_diff = std::max(max_diff, std::abs(lhs.phi_imag - rhs.phi_imag));
    }

    return max_diff;
}

inline void collectPenaltyMapsFromBody(std::map<AphiContactPositionKey, Real> &div_a_map,
                                       std::map<AphiContactPositionKey, Vecd> &grad_div_a_map,
                                       AphiBlockMapByPosition &lhs_map, BaseParticles &particles,
                                       const AphiVariableNames &names, const AphiADivergencePenaltyScratchNames &scratch,
                                       Real body_length, Real body_height, Real body_width, Real core_shell,
                                       Real x_interface, Real dp_0)
{
    collectCoreRealByPosition(div_a_map, particles, scratch.div_a_real, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
    collectCoreVecdByPosition(grad_div_a_map, particles, scratch.grad_div_a_real, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
    collectCoreBlockByPosition(lhs_map, particles, names.lhs, body_length, body_height, body_width, core_shell,
                               x_interface, dp_0);
}

inline void runMonolithicADivergencePenaltyApply(
    std::map<AphiContactPositionKey, Real> &div_a_map, std::map<AphiContactPositionKey, Vecd> &grad_div_a_map,
    AphiBlockMapByPosition &lhs_map, int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
    Real boundary_width, Real core_shell, Real x_interface, Real sigma_left, Real sigma_right, Real nu,
    Real lambda_a, std::map<AphiContactPositionKey, Real> *interior_div_a_map = nullptr)
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));
    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    SolidBody body(sph_system, makeShared<AphiHalfSpaceBoxShape>("MonolithicBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    Inner<> inner_ck(body);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    const AphiADivergencePenaltyScratchNames scratch = aphiDefaultADivergencePenaltyScratchNames();
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignPiecewiseSigmaHalfSpaceCK> assign_sigma(
        body, x_interface, sigma_left, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(body, names);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_lhs(body, names.lhs);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    AphiPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> penalty_pipeline(
        body, inner_ck, names.solution, names.lhs, lambda_a, scratch);

    initialize_aphi.exec();
    set_material.exec();
    assign_sigma.exec();
    assign_fields.exec();
    zero_lhs.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    penalty_pipeline.exec();

    collectPenaltyMapsFromBody(div_a_map, grad_div_a_map, lhs_map, body.getBaseParticles(), names, scratch, body_length,
                               body_height, body_width, core_shell, x_interface, dp_0);
    if (interior_div_a_map != nullptr)
    {
        collectMatchedRealByPosition(*interior_div_a_map, body.getBaseParticles(), scratch.div_a_real, body_length,
                                     body_height, body_width, core_shell, x_interface, dp_0, false, true);
    }
}

inline void runSplitContactADivergencePenaltyApply(
    std::map<AphiContactPositionKey, Real> &div_a_map, std::map<AphiContactPositionKey, Vecd> &grad_div_a_map,
    AphiBlockMapByPosition &lhs_map, int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
    Real boundary_width, Real core_shell, Real x_interface, Real sigma_left, Real sigma_right, Real nu,
    Real lambda_a, std::map<AphiContactPositionKey, Real> *interior_div_a_map = nullptr)
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));
    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd left_center(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd right_center(0.75 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);

    SolidBody left_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("LeftBody", left_center, halfsize));
    SolidBody right_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("RightBody", right_center, halfsize));
    for (auto *body_ptr : {&left_body, &right_body})
    {
        body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        body_ptr->defineMaterial<Solid>();
        body_ptr->defineBodyLevelSetShape();
        body_ptr->generateParticles<BaseParticles, Lattice>();
    }

    Inner<> left_inner(left_body);
    Inner<> right_inner(right_body);
    Contact<> left_to_right(left_body, {&right_body});
    Contact<> right_to_left(right_body, {&left_body});

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    const AphiADivergencePenaltyScratchNames scratch = aphiDefaultADivergencePenaltyScratchNames();
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, sigma_right, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(left_body, sigma_left, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(right_body, sigma_right, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_left_fields(left_body, names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_right_fields(right_body, names);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_lhs(left_body, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_lhs(right_body, names.lhs);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_to_right);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_to_left);

    AphiContactPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> left_penalty(
        left_body, left_inner, left_to_right, names.solution, names.lhs, lambda_a, scratch);
    AphiContactPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> right_penalty(
        right_body, right_inner, right_to_left, names.solution, names.lhs, lambda_a, scratch);

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    assign_left_fields.exec();
    assign_right_fields.exec();
    zero_left_lhs.exec();
    zero_right_lhs.exec();
    update_left_cell_linked_list.exec();
    update_right_cell_linked_list.exec();
    update_left_inner.exec();
    update_right_inner.exec();
    update_left_contact.exec();
    update_right_contact.exec();
    left_penalty.exec();
    right_penalty.exec();

    for (auto *body_ptr : {&left_body, &right_body})
    {
        std::map<AphiContactPositionKey, Real> body_div_a;
        std::map<AphiContactPositionKey, Vecd> body_grad_div_a;
        AphiBlockMapByPosition body_lhs;
        collectPenaltyMapsFromBody(body_div_a, body_grad_div_a, body_lhs, body_ptr->getBaseParticles(), names, scratch,
                                   body_length, body_height, body_width, core_shell, x_interface, dp_0);
        for (const auto &entry : body_div_a)
        {
            div_a_map[entry.first] = entry.second;
        }
        for (const auto &entry : body_grad_div_a)
        {
            grad_div_a_map[entry.first] = entry.second;
        }
        for (const auto &entry : body_lhs)
        {
            lhs_map[entry.first] = entry.second;
        }
        if (interior_div_a_map != nullptr)
        {
            std::map<AphiContactPositionKey, Real> body_interior_div_a;
            collectMatchedRealByPosition(body_interior_div_a, body_ptr->getBaseParticles(), scratch.div_a_real,
                                         body_length, body_height, body_width, core_shell, x_interface, dp_0, false,
                                         true);
            for (const auto &entry : body_interior_div_a)
            {
                (*interior_div_a_map)[entry.first] = entry.second;
            }
        }
    }
}

inline AphiContactADivergencePenaltyEquivalenceMetrics runContactADivergencePenaltyEquivalenceMetrics(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = dp_0;
    const Real x_interface = 0.5;
    const Real sigma_left = 10.0;
    const Real sigma_right = 1.0;
    const Real nu = 1.5;
    const Real lambda_a = 12.5;

    std::map<AphiContactPositionKey, Real> mono_div_a;
    std::map<AphiContactPositionKey, Vecd> mono_grad_div_a;
    AphiBlockMapByPosition mono_lhs;
    std::map<AphiContactPositionKey, Real> contact_div_a;
    std::map<AphiContactPositionKey, Vecd> contact_grad_div_a;
    AphiBlockMapByPosition contact_lhs;
    runMonolithicADivergencePenaltyApply(mono_div_a, mono_grad_div_a, mono_lhs, ac, av, dp_0, body_length, body_height,
                                         body_width, boundary_width, core_shell, x_interface, sigma_left, sigma_right,
                                         nu, lambda_a);
    runSplitContactADivergencePenaltyApply(contact_div_a, contact_grad_div_a, contact_lhs, ac, av, dp_0, body_length,
                                           body_height, body_width, boundary_width, core_shell, x_interface,
                                           sigma_left, sigma_right, nu, lambda_a);

    AphiContactADivergencePenaltyEquivalenceMetrics metrics;
    metrics.div_a_max_abs_diff =
        maxAbsScalarDifference(mono_div_a, contact_div_a, metrics.div_a_matched, metrics.div_a_missing);
    metrics.grad_div_a_stencil_safe_max_abs_diff = maxAbsVecdMapDifferenceOnStencilSafeCore(
        mono_grad_div_a, contact_grad_div_a, x_interface, dp_0, metrics.grad_div_a_stencil_safe_matched,
        metrics.grad_div_a_stencil_safe_missing, metrics.grad_div_a_interface_band_matched,
        metrics.grad_div_a_interface_band_max_abs_diff);
    metrics.lhs_a_stencil_safe_max_abs_diff = maxAbsBlockDifferenceOnStencilSafeCore(
        mono_lhs, contact_lhs, x_interface, dp_0, metrics.lhs_a_stencil_safe_matched, metrics.lhs_a_stencil_safe_missing);
    return metrics;
}

inline void printContactADivergencePenaltyEquivalenceMetrics(const std::string &prefix,
                                                             const AphiContactADivergencePenaltyEquivalenceMetrics &m)
{
    std::cout << prefix << " div_a_matched=" << m.div_a_matched << " div_a_missing=" << m.div_a_missing
              << " div_a_max_abs_diff=" << m.div_a_max_abs_diff
              << " grad_div_a_stencil_safe_matched=" << m.grad_div_a_stencil_safe_matched
              << " grad_div_a_stencil_safe_missing=" << m.grad_div_a_stencil_safe_missing
              << " grad_div_a_stencil_safe_max_abs_diff=" << m.grad_div_a_stencil_safe_max_abs_diff
              << " grad_div_a_interface_band_matched=" << m.grad_div_a_interface_band_matched
              << " grad_div_a_interface_band_max_abs_diff=" << m.grad_div_a_interface_band_max_abs_diff
              << " lhs_a_stencil_safe_matched=" << m.lhs_a_stencil_safe_matched
              << " lhs_a_stencil_safe_missing=" << m.lhs_a_stencil_safe_missing
              << " lhs_a_stencil_safe_max_abs_diff=" << m.lhs_a_stencil_safe_max_abs_diff << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_A_DIVERGENCE_PENALTY_DIAGNOSTIC_HELPERS_H
