/**
 * Stage 10.11 follow-up: B=curl A dual-track diagnostic (B-corrected grad vs pairwise uncorrected grad).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

struct CurlBDualTrackRow
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real b_corrected_core_rel = 0.0;
    Real b_pairwise_core_rel = 0.0;
};

CurlBDualTrackRow runCurlBDualTrackRow(int ac, char *av[], Real dp_0, Real body_length, Real body_height,
                                       Real body_width, Real boundary_width, Real core_shell, Real sigma, Real nu,
                                       AphiDivFreeValidationFieldKind field_kind)
{
    CurlBDualTrackRow row;
    row.field_kind = field_kind;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_field(
        test_body.body, names.solution, field_kind);

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();

    const std::string b_corrected_real = "BCorrectedReal";
    const std::string b_corrected_imag = "BCorrectedImag";
    const std::string b_pairwise_real = "BPairwiseReal";
    const std::string b_pairwise_imag = "BPairwiseImag";

    execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, b_corrected_real, b_corrected_imag,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, b_pairwise_real, b_pairwise_imag,
                                 AphiBCurlDiagnosticMode::PairwiseUncorrectedGrad);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    row.b_corrected_core_rel =
        hostCoreDivFreeBRelativeError(particles, b_corrected_real, b_corrected_imag, positions, total_real_particles,
                                      body_length, body_height, body_width, core_shell, field_kind);
    row.b_pairwise_core_rel =
        hostCoreDivFreeBRelativeError(particles, b_pairwise_real, b_pairwise_imag, positions, total_real_particles,
                                      body_length, body_height, body_width, core_shell, field_kind);
    return row;
}

void printCurlBDualTrackRow(const std::string &prefix, Real dp, Real core_shell, const CurlBDualTrackRow &row)
{
    const Real ratio = row.b_corrected_core_rel / (row.b_pairwise_core_rel + TinyReal);
    std::cout << prefix << " field=" << divFreeValidationFieldName(row.field_kind) << " dp=" << dp
              << " core_shell=" << core_shell << " b_corrected_core_rel=" << row.b_corrected_core_rel
              << " b_pairwise_core_rel=" << row.b_pairwise_core_rel << " b_corrected_vs_pairwise_ratio=" << ratio
              << std::endl;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real max_linear2d_b_corrected_rel = 0.15;
    const Real min_oscillatory_b_corrected_rel = 0.02;
    const Real max_oscillatory_b_corrected_rel = 0.08;
    const Real max_b_corrected_vs_pairwise_ratio = 0.05;

    const CurlBDualTrackRow linear2d = runCurlBDualTrackRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu,
        AphiDivFreeValidationFieldKind::Linear2D);
    const CurlBDualTrackRow az2d = runCurlBDualTrackRow(ac, av, dp_0, body_length, body_height, body_width,
                                                        boundary_width, core_shell, sigma, nu,
                                                        AphiDivFreeValidationFieldKind::Az2D);
    const CurlBDualTrackRow crosssine = runCurlBDualTrackRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu,
        AphiDivFreeValidationFieldKind::CrossSine3D);

    printCurlBDualTrackRow("test_3d_aphi_ck_curl_b_dual_track_diagnostic", dp_0, core_shell, linear2d);
    printCurlBDualTrackRow("test_3d_aphi_ck_curl_b_dual_track_diagnostic", dp_0, core_shell, az2d);
    printCurlBDualTrackRow("test_3d_aphi_ck_curl_b_dual_track_diagnostic", dp_0, core_shell, crosssine);

    const auto oscillatory_dual_track_ok = [&](const CurlBDualTrackRow &row) {
        const Real ratio = row.b_corrected_core_rel / (row.b_pairwise_core_rel + TinyReal);
        return row.b_corrected_core_rel >= min_oscillatory_b_corrected_rel &&
               row.b_corrected_core_rel <= max_oscillatory_b_corrected_rel && row.b_pairwise_core_rel > row.b_corrected_core_rel &&
               ratio <= max_b_corrected_vs_pairwise_ratio;
    };

    const bool linear2d_ok = linear2d.b_corrected_core_rel <= max_linear2d_b_corrected_rel;
    const bool az2d_ok = oscillatory_dual_track_ok(az2d);
    const bool crosssine_ok = oscillatory_dual_track_ok(crosssine);
    const bool passed = linear2d_ok && az2d_ok && crosssine_ok;

    std::cout << "test_3d_aphi_ck_curl_b_dual_track_diagnostic linear2d_ok=" << (linear2d_ok ? 1 : 0)
              << " az2d_ok=" << (az2d_ok ? 1 : 0) << " crosssine_ok=" << (crosssine_ok ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
