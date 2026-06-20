/**
 * @file test_3d_ophelie_phi_rhs_flux_sign_audit.cpp
 * @brief Linear A(x): legacy-flux RHS ≈ −div(σA) RHS (same uncorrected div kernel, constant sigma).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
class OphelieRhsSignAuditGlassBoxShape : public ComplexShape
{
  public:
    OphelieRhsSignAuditGlassBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};
} // namespace

int main(int, char *[])
{
    OphelieParameters params;
    params.frequency_ = 300.0e3;
    params.sigma_glass_ = 16.0;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;

    const Real a_grad = 0.8;
    const Real dp = 0.06;
    const Vecd center(0.0, 0.0, 0.25);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieRhsSignAuditGlassBoxShape>("GlassBody", center, halfsize));
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

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_src_real = particles.getVariableDataByName<Vecd>(glass_names.a_src_real);
    for (size_t i = 0; i < n; ++i)
    {
        a_src_real[i] = Vecd(a_grad * pos[i][0], 0.0, 0.0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_src_real);

    StdVec<Real> rhs_div(n, Real(0));
    StdVec<Real> rhs_legacy(n, Real(0));
    OphelieParameters div_params = params;
    div_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_params);
    hostReadScalarField(particles, glass_names.phi_rhs_imag, rhs_div.data(), n);
    OphelieParameters legacy_params = params;
    legacy_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, legacy_params);
    hostReadScalarField(particles, glass_names.phi_rhs_imag, rhs_legacy.data(), n);

    const OpheliePhiRhsAlignmentMetrics alignment =
        measurePhiRhsAlignmentMetrics(particles, rhs_div.data(), rhs_legacy.data(), n);

    const bool passed = n > 0 && alignment.rhs_cosine_div_neg_legacy > 0.95 &&
                        alignment.rhs_div_vs_neg_legacy_vol < 0.35;

    std::cout << "test_3d_ophelie_phi_rhs_flux_sign_audit n=" << n << " rhs_cosine_div_legacy="
              << alignment.rhs_cosine_div_legacy << " rhs_cosine_div_neg_legacy="
              << alignment.rhs_cosine_div_neg_legacy << " rhs_div_vs_neg_legacy_vol="
              << alignment.rhs_div_vs_neg_legacy_vol << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
