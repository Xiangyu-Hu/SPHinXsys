/**
 * @file test_3d_ophelie_phi_laplace_constant.cpp
 * @brief Constant PhiImag -> pairwise Laplace operator ~ 0; diagonal positive.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

int main(int, char *[])
{
    OphelieParameters params;
    params.sigma_glass_ = 16.0;
    const Real phi_constant = 2.5;

    const Real dp = 0.08;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        phi_imag[i] = phi_constant;
    }
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);

    StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_lhs(glass_body, glass_names.phi_lhs_imag);
    StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_diag(glass_body, glass_names.phi_laplace_diag);
    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> apply_laplace(
        *glass_inner, glass_names.phi_imag, glass_names.sigma, glass_names.phi_lhs_imag,
        params.pair_weight_regularization_);
    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_diag(
        *glass_inner, glass_names.sigma, glass_names.phi_laplace_diag, params.pair_weight_regularization_);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);

    zero_lhs.exec();
    zero_diag.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    apply_laplace.exec();
    compute_diag.exec();

    syncVariableToHost<Real>(particles, glass_names.phi_lhs_imag);
    syncVariableToHost<Real>(particles, glass_names.phi_laplace_diag);
    const Real *laplace = particles.getVariableDataByName<Real>(glass_names.phi_lhs_imag);
    const Real *diag = particles.getVariableDataByName<Real>(glass_names.phi_laplace_diag);

    Real max_laplace = 0.0;
    Real min_diag = diag[0];
    for (size_t i = 0; i < n; ++i)
    {
        max_laplace = std::max(max_laplace, std::abs(laplace[i]));
        min_diag = std::min(min_diag, diag[i]);
    }
    const Real laplace_rel = max_laplace / (std::abs(phi_constant) + TinyReal);

    const bool passed = n > 0 && laplace_rel < 0.05 && min_diag > 0.0;
    std::cout << "test_3d_ophelie_phi_laplace_constant n=" << n << " phi=" << phi_constant
              << " max_laplace=" << max_laplace << " laplace_rel=" << laplace_rel << " min_diag=" << min_diag
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
