/**
 * @file test_3d_ophelie_vector_divergence_mms.cpp
 * @brief Continuous vector-field divergence MMS: D_unc(A) vs D_c(A) vs div A_exact.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_vector_divergence_diagnostics.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

constexpr OphelieVectorDivergenceMmsCase kAllCases[] = {
    OphelieVectorDivergenceMmsCase::ConstantX,
    OphelieVectorDivergenceMmsCase::LinearXYZ,
    OphelieVectorDivergenceMmsCase::RotationalXY,
    OphelieVectorDivergenceMmsCase::QuadraticXYZ,
    OphelieVectorDivergenceMmsCase::SingularToroidalLegacy,
    OphelieVectorDivergenceMmsCase::SmoothToroidal,
};

inline bool isFiniteMetrics(const OphelieVectorDivergenceErrorMetrics &metrics)
{
    return std::isfinite(metrics.l2_all) && std::isfinite(metrics.l2_interior) &&
           std::isfinite(metrics.l2_boundary) && std::isfinite(metrics.linf_all) &&
           std::isfinite(metrics.flipped_l2_interior) && std::isfinite(metrics.best_sign_l2_interior) &&
           std::isfinite(metrics.sign_alpha);
}

inline bool isFiniteReport(const OphelieVectorDivergenceMmsCaseReport &report)
{
    return isFiniteMetrics(report.uncorrected) && isFiniteMetrics(report.corrected) &&
           std::isfinite(report.corrected_vs_uncorrected_l2) && std::isfinite(report.corrected_vs_uncorrected_cosine);
}

inline bool passSmokeChecks(const OphelieVectorDivergenceMmsCaseReport &report)
{
    if (!isFiniteReport(report))
    {
        return false;
    }
    if (report.kind == OphelieVectorDivergenceMmsCase::ConstantX)
    {
        return report.uncorrected.l2_all < Real(1.0e-5) && report.corrected.l2_all < Real(1.0e-5);
    }
    return true;
}

inline bool passValidationChecks(const OphelieVectorDivergenceMmsCaseReport &report)
{
    const Real linear_tol = Real(0.2);
    const Real rotational_tol = Real(0.05);
    const Real smooth_toroidal_tol = Real(0.1);
    switch (report.kind)
    {
    case OphelieVectorDivergenceMmsCase::LinearXYZ:
        return report.uncorrected.best_sign_l2_interior < linear_tol &&
               report.corrected.best_sign_l2_interior < linear_tol;
    case OphelieVectorDivergenceMmsCase::RotationalXY:
        return report.uncorrected.l2_interior < rotational_tol && report.corrected.l2_interior < rotational_tol;
    case OphelieVectorDivergenceMmsCase::SmoothToroidal:
        return report.uncorrected.l2_interior < smooth_toroidal_tol && report.corrected.l2_interior < smooth_toroidal_tol;
    default:
        return true;
    }
}

inline bool runAtDp(Real dp, const OphelieTestCliOptions &cli_options, size_t &n_particles, bool &validation_passed)
{
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.40, 0.40, 0.40);
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
    n_particles = glass_body.getBaseParticles().TotalRealParticles();

    bool smoke_passed = n_particles > 0;
    for (const OphelieVectorDivergenceMmsCase kind : kAllCases)
    {
        const OphelieVectorDivergenceMmsCaseReport report =
            evaluateOphelieVectorDivergenceMmsCaseBothOperators<MainExecutionPolicy>(
                glass_body, *glass_inner, glass_names, kind, center, halfsize, dp);
        logOphelieVectorDivergenceMmsCaseReport(dp, report);
        appendOphelieVectorDivergenceMmsCsv(cli_options.vector_divergence_mms_csv_path, dp, report);
        smoke_passed = smoke_passed && passSmokeChecks(report);
        validation_passed = validation_passed && passValidationChecks(report);
    }
    return smoke_passed;
}

} // namespace

int main(int ac, char *av[])
{
    OphelieParameters params;
    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);

    std::vector<Real> dp_list;
    if (cli_options.vector_divergence_mms_dp > 0.0)
    {
        dp_list.push_back(cli_options.vector_divergence_mms_dp);
    }
    else if (cli_options.vector_divergence_mms_scan)
    {
        dp_list = {0.04, 0.03, 0.02, 0.015};
    }
    else
    {
        dp_list = {0.06};
    }

    bool smoke_passed = true;
    bool validation_passed = true;
    size_t n_particles = 0;
    for (const Real dp : dp_list)
    {
        smoke_passed = runAtDp(dp, cli_options, n_particles, validation_passed) && smoke_passed;
    }

    std::cout << "test_3d_ophelie_vector_divergence_mms n=" << n_particles << " dp_count=" << dp_list.size()
              << " smoke_passed=" << (smoke_passed ? 1 : 0) << " validation_passed=" << (validation_passed ? 1 : 0)
              << std::endl;
    return smoke_passed ? 0 : 1;
}
