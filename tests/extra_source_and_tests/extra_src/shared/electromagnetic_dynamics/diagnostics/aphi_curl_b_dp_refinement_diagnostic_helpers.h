#ifndef APHI_CURL_B_DP_REFINEMENT_DIAGNOSTIC_HELPERS_H
#define APHI_CURL_B_DP_REFINEMENT_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiCurlBDpRefinementRow
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real core_shell = 0.0;
    size_t total_real_particles = 0;
    Real b_error_vs_exact = 0.0;
    Real E_combined_error_vs_exact = 0.0;
    Real joule_error_vs_exact = 0.0;
    Real core_div_a_grad_den_rel = 0.0;
    Real boundary_div_a_grad_den_rel = 0.0;
};

inline Vecd divFreeValidationEExactRealField(AphiDivFreeValidationFieldKind kind, Real omega, Real x, Real y, Real z)
{
    (void)kind;
    (void)omega;
    (void)x;
    (void)y;
    (void)z;
    return Vecd(0.0, 0.0, 0.0);
}

inline Vecd divFreeValidationEExactImagField(AphiDivFreeValidationFieldKind kind, Real omega, Real x, Real y, Real z)
{
    const Vecd a_real = divFreeValidationARealField(kind, x, y, z);
    (void)z;
    return -omega * a_real;
}

inline Real hostCoreDivFreeECombinedRelativeError(
    BaseParticles &particles, const AphiJouleHeatingFieldNames &joule_fields, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell, Real omega,
    AphiDivFreeValidationFieldKind field_kind)
{
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_real);
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *e_real = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_real);
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real error_squared = 0.0;
    Real reference_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Vecd exact_e_re =
            divFreeValidationEExactRealField(field_kind, omega, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd exact_e_im =
            divFreeValidationEExactImagField(field_kind, omega, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd err_re = e_real[i] - exact_e_re;
        const Vecd err_im = e_imag[i] - exact_e_im;
        error_squared += vol[i] * (err_re.squaredNorm() + err_im.squaredNorm());
        reference_squared += vol[i] * (exact_e_re.squaredNorm() + exact_e_im.squaredNorm());
    }
    return std::sqrt(error_squared) / (std::sqrt(reference_squared) + TinyReal);
}

inline AphiCurlBDpRefinementRow runCurlBDpRefinementRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, AphiDivFreeValidationFieldKind field_kind)
{
    AphiCurlBDpRefinementRow row;
    row.field_kind = field_kind;
    row.dp = dp_0;
    row.core_shell = core_shell;

    const Real boundary_width = 3.0 * dp_0;
    const std::string b_real_name = "DiagnosticBReal";
    const std::string b_imag_name = "DiagnosticBImag";

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_field(
        test_body.body, names.solution, field_kind);

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();

    execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, b_real_name, b_imag_name,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::BCorrectedTrace);
    execInnerJoulePostProcess(test_body, names, omega, joule_fields);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    row.total_real_particles = total_real_particles;

    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };

    row.b_error_vs_exact =
        hostCoreDivFreeBRelativeError(particles, b_real_name, b_imag_name, positions, total_real_particles,
                                      body_length, body_height, body_width, core_shell, field_kind);
    row.E_combined_error_vs_exact =
        hostCoreDivFreeECombinedRelativeError(particles, joule_fields, positions, total_real_particles, body_length,
                                              body_height, body_width, core_shell, omega, field_kind);
    row.core_div_a_grad_den_rel =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_core,
                                     AphiDivADiagnosticMode::BCorrectedTrace)
            .div_a_relative;
    row.boundary_div_a_grad_den_rel =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_boundary,
                                     AphiDivADiagnosticMode::BCorrectedTrace)
            .div_a_relative;

    const AphiRegionalElectromagneticObservables observables =
        hostCoreObservablesFromJouleFields(particles, names, joule_fields, positions, total_real_particles,
                                           body_length, body_height, body_width, core_shell, sigma);
    Real analytic_joule_power = 0.0;
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Vecd exact_e_re =
            divFreeValidationEExactRealField(field_kind, omega, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd exact_e_im =
            divFreeValidationEExactImagField(field_kind, omega, positions[i][0], positions[i][1], positions[i][2]);
        analytic_joule_power += 0.5 * sigma * vol[i] * (exact_e_re.squaredNorm() + exact_e_im.squaredNorm());
    }
    row.joule_error_vs_exact =
        std::abs(observables.Joule_power - analytic_joule_power) / (std::abs(analytic_joule_power) + TinyReal);

    return row;
}

inline void printCurlBDpRefinementRow(const std::string &prefix, const AphiCurlBDpRefinementRow &row)
{
    std::cout << prefix << " field=" << divFreeValidationFieldName(row.field_kind) << " dp=" << row.dp
              << " N=" << row.total_real_particles << " core_shell=" << row.core_shell
              << " B_err_vs_exact=" << row.b_error_vs_exact << " E_combined_err_vs_exact=" << row.E_combined_error_vs_exact
              << " Joule_err_vs_exact=" << row.joule_error_vs_exact << " core_divA_gradDen_rel=" << row.core_div_a_grad_den_rel
              << " boundary_divA_gradDen_rel=" << row.boundary_div_a_grad_den_rel << std::endl;
}

inline const AphiCurlBDpRefinementRow *findCurlBDpRefinementRowByDp(const std::vector<AphiCurlBDpRefinementRow> &rows,
                                                                  Real dp)
{
    for (const AphiCurlBDpRefinementRow &row : rows)
    {
        if (std::abs(row.dp - dp) <= 1.0e-12)
        {
            return &row;
        }
    }
    return nullptr;
}

inline bool curlBDpRefinementOscillatoryConverges(const std::vector<AphiCurlBDpRefinementRow> &rows, Real max_fine_b_error)
{
    if (rows.size() < 2)
    {
        return false;
    }
    const size_t medium = rows.size() - 2;
    const size_t fine = rows.size() - 1;
    return rows[fine].b_error_vs_exact + TinyReal < rows[medium].b_error_vs_exact &&
           rows[fine].b_error_vs_exact <= max_fine_b_error;
}

/** Stage 10.13-P3: Az2D B=curlA dp gate on {0.10, 0.075, 0.05}. */
inline bool curlBDpRefinementAz2DStage1013Gate(const std::vector<AphiCurlBDpRefinementRow> &rows,
                                               Real max_dp005_b_error = 0.02, Real max_dp075_over_dp005_ratio = 1.2)
{
    const AphiCurlBDpRefinementRow *const row_dp10 = findCurlBDpRefinementRowByDp(rows, 0.10);
    const AphiCurlBDpRefinementRow *const row_dp075 = findCurlBDpRefinementRowByDp(rows, 0.075);
    const AphiCurlBDpRefinementRow *const row_dp05 = findCurlBDpRefinementRowByDp(rows, 0.05);
    if (row_dp10 == nullptr || row_dp075 == nullptr || row_dp05 == nullptr)
    {
        return false;
    }
    return row_dp075->b_error_vs_exact + TinyReal < row_dp10->b_error_vs_exact &&
           row_dp05->b_error_vs_exact <= row_dp075->b_error_vs_exact * max_dp075_over_dp005_ratio &&
           row_dp05->b_error_vs_exact < max_dp005_b_error;
}

inline bool curlBDpRefinementLinear2DStable(const std::vector<AphiCurlBDpRefinementRow> &rows, Real max_b_error)
{
    for (const AphiCurlBDpRefinementRow &row : rows)
    {
        if (row.b_error_vs_exact > max_b_error)
        {
            return false;
        }
    }
    return true;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CURL_B_DP_REFINEMENT_DIAGNOSTIC_HELPERS_H
