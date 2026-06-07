#ifndef APHI_H_CURL_H_J_OBSERVABLE_DIAGNOSTIC_HELPERS_H
#define APHI_H_CURL_H_J_OBSERVABLE_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dp_refinement_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

class AphiScaleVecdFieldByRealCK : public LocalDynamics
{
  public:
    AphiScaleVecdFieldByRealCK(SPHBody &sph_body, const std::string &input_name, const std::string &output_name,
                               Real scale)
        : LocalDynamics(sph_body), scale_(scale),
          dv_input_(particles_->template getVariableByName<Vecd>(input_name)),
          dv_output_(particles_->template getVariableByName<Vecd>(output_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : scale_(encloser.scale_), input_(encloser.dv_input_->DelegatedData(ex_policy)),
              output_(encloser.dv_output_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            output_[index_i] = scale_ * input_[index_i];
        }

      protected:
        Real scale_;
        Vecd *input_;
        Vecd *output_;
    };

  protected:
    Real scale_;
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Vecd> *dv_output_;
};

inline Vecd divFreeValidationHRealField(AphiDivFreeValidationFieldKind kind, Real nu, Real x, Real y, Real z)
{
    return nu * divFreeValidationBRealField(kind, x, y, z);
}

inline Vecd divFreeValidationJExactImagField(AphiDivFreeValidationFieldKind kind, Real sigma, Real omega, Real x, Real y,
                                             Real z)
{
    return sigma * divFreeValidationEExactImagField(kind, omega, x, y, z);
}

inline Vecd divFreeValidationCurlHRealField(AphiDivFreeValidationFieldKind kind, Real nu, Real x, Real y, Real z)
{
    const Real pi = Pi;
    const Real sx = std::sin(pi * x);
    const Real sy = std::sin(pi * y);
    if (kind == AphiDivFreeValidationFieldKind::Linear2D)
    {
        return Vecd::Zero();
    }
    if (kind == AphiDivFreeValidationFieldKind::Az2D)
    {
        return Vecd(0.0, 0.0, 2.0 * nu * pi * pi * sx * sy);
    }
    (void)z;
    return Vecd::Zero();
}

inline void execBodyCurlFromVectorFieldDiagnostic(SPHBody &body, Inner<> &inner, const std::string &field_real,
                                                  const std::string &field_imag, const std::string &curl_real,
                                                  const std::string &curl_imag, AphiBCurlDiagnosticMode mode)
{
    const std::string grad_real_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? field_real + "Gradient" : field_real + "GradientPairwise";
    const std::string grad_imag_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? field_imag + "Gradient" : field_imag + "GradientPairwise";

    BaseParticles &particles = body.getBaseParticles();
    if (mode == AphiBCurlDiagnosticMode::PairwiseUncorrectedGrad)
    {
        particles.registerStateVariable<Matd>(grad_real_name, ZeroData<Matd>::value);
        particles.registerStateVariable<Matd>(grad_imag_name, ZeroData<Matd>::value);
        InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorGradientCK<Inner<>>> field_real_gradient(
            inner, field_real, grad_real_name);
        InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorGradientCK<Inner<>>> field_imag_gradient(
            inner, field_imag, grad_imag_name);
        field_real_gradient.exec();
        field_imag_gradient.exec();
    }
    else
    {
        InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(
            inner);
        InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> field_real_gradient(inner, field_real);
        InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> field_imag_gradient(inner, field_imag);
        linear_correction_matrix.exec();
        field_real_gradient.exec();
        field_imag_gradient.exec();
    }

    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_real_kernel(body, grad_real_name, curl_real);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_imag_kernel(body, grad_imag_name, curl_imag);
    curl_real_kernel.exec();
    curl_imag_kernel.exec();
}

inline Real hostCoreVectorFieldRelativeError(
    BaseParticles &particles, const std::string &field_real_name, const std::string &field_imag_name,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell,
    const std::function<Vecd(AphiDivFreeValidationFieldKind, Real, Real, Real)> &exact_real,
    const std::function<Vecd(AphiDivFreeValidationFieldKind, Real, Real, Real)> &exact_imag,
    AphiDivFreeValidationFieldKind field_kind)
{
    syncVariableToHost<Vecd>(particles, field_real_name);
    syncVariableToHost<Vecd>(particles, field_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *field_real = particles.getVariableDataByName<Vecd>(field_real_name);
    const Vecd *field_imag = particles.getVariableDataByName<Vecd>(field_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real error_squared = 0.0;
    Real reference_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Vecd exact_re = exact_real(field_kind, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd exact_im = exact_imag(field_kind, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd err_re = field_real[i] - exact_re;
        const Vecd err_im = field_imag[i] - exact_im;
        error_squared += vol[i] * (err_re.squaredNorm() + err_im.squaredNorm());
        reference_squared += vol[i] * (exact_re.squaredNorm() + exact_im.squaredNorm());
    }
    return std::sqrt(error_squared) / (std::sqrt(reference_squared) + TinyReal);
}

struct AphiHCurlHJObservableRow
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real core_shell = 0.0;
    Real b_error_vs_exact = 0.0;
    Real h_error_vs_exact = 0.0;
    Real h_to_b_error_ratio = 0.0;
    Real curl_h_error_vs_exact = 0.0;
    Real j_error_vs_exact = 0.0;
    Real curl_h_to_j_error_ratio = 0.0;
    Real discrete_curl_h_to_j_l2_ratio = 0.0;
};

inline Real hostCoreCombinedVectorFieldL2(BaseParticles &particles, const std::string &field_real_name,
                                          const std::string &field_imag_name, const Vecd *positions,
                                          size_t total_real_particles, Real body_length, Real body_height,
                                          Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, field_real_name);
    syncVariableToHost<Vecd>(particles, field_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *field_real = particles.getVariableDataByName<Vecd>(field_real_name);
    const Vecd *field_imag = particles.getVariableDataByName<Vecd>(field_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        sum_squared += vol[i] * (field_real[i].squaredNorm() + field_imag[i].squaredNorm());
    }
    return std::sqrt(sum_squared);
}

inline AphiHCurlHJObservableRow runHCurlHJObservableRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, AphiDivFreeValidationFieldKind field_kind)
{
    AphiHCurlHJObservableRow row;
    row.field_kind = field_kind;
    row.dp = dp_0;
    row.core_shell = core_shell;

    const Real boundary_width = 3.0 * dp_0;
    const std::string b_real_name = "DiagnosticBReal";
    const std::string b_imag_name = "DiagnosticBImag";
    const std::string h_real_name = "DiagnosticHReal";
    const std::string h_imag_name = "DiagnosticHImag";
    const std::string curl_h_real_name = "DiagnosticCurlHReal";
    const std::string curl_h_imag_name = "DiagnosticCurlHImag";

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    BaseParticles &particles = test_body.body.getBaseParticles();
    particles.registerStateVariable<Vecd>(b_real_name, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(b_imag_name, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(h_real_name, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(h_imag_name, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(curl_h_real_name, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(curl_h_imag_name, ZeroData<Vecd>::value);

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
    StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByRealCK> scale_h_from_b_real(
        test_body.body, b_real_name, h_real_name, nu);
    StateDynamics<MainExecutionPolicy, AphiScaleVecdFieldByRealCK> scale_h_from_b_imag(
        test_body.body, b_imag_name, h_imag_name, nu);
    scale_h_from_b_real.exec();
    scale_h_from_b_imag.exec();
    execBodyCurlFromVectorFieldDiagnostic(test_body.body, test_body.inner(), h_real_name, h_imag_name, curl_h_real_name,
                                          curl_h_imag_name, AphiBCurlDiagnosticMode::BCorrectedGrad);
    execInnerJoulePostProcess(test_body, names, omega, joule_fields);

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    row.b_error_vs_exact =
        hostCoreDivFreeBRelativeError(particles, b_real_name, b_imag_name, positions, total_real_particles, body_length,
                                      body_height, body_width, core_shell, field_kind);
    row.h_error_vs_exact = hostCoreVectorFieldRelativeError(
        particles, h_real_name, h_imag_name, positions, total_real_particles, body_length, body_height, body_width,
        core_shell,
        [&](AphiDivFreeValidationFieldKind kind, Real x, Real y, Real z) {
            return divFreeValidationHRealField(kind, nu, x, y, z);
        },
        [&](AphiDivFreeValidationFieldKind, Real, Real, Real) { return Vecd::Zero(); }, field_kind);
    row.h_to_b_error_ratio = row.h_error_vs_exact / (row.b_error_vs_exact + TinyReal);

    row.curl_h_error_vs_exact = hostCoreVectorFieldRelativeError(
        particles, curl_h_real_name, curl_h_imag_name, positions, total_real_particles, body_length, body_height,
        body_width, core_shell,
        [&](AphiDivFreeValidationFieldKind kind, Real x, Real y, Real z) {
            return divFreeValidationCurlHRealField(kind, nu, x, y, z);
        },
        [&](AphiDivFreeValidationFieldKind, Real, Real, Real) { return Vecd::Zero(); }, field_kind);
    row.j_error_vs_exact = hostCoreVectorFieldRelativeError(
        particles, joule_fields.current_density_real, joule_fields.current_density_imag, positions, total_real_particles,
        body_length, body_height, body_width, core_shell,
        [&](AphiDivFreeValidationFieldKind, Real, Real, Real) { return Vecd::Zero(); },
        [&](AphiDivFreeValidationFieldKind kind, Real x, Real y, Real z) {
            return divFreeValidationJExactImagField(kind, sigma, omega, x, y, z);
        },
        field_kind);
    row.curl_h_to_j_error_ratio = row.curl_h_error_vs_exact / (row.j_error_vs_exact + TinyReal);
    row.discrete_curl_h_to_j_l2_ratio =
        hostCoreCombinedVectorFieldL2(particles, curl_h_real_name, curl_h_imag_name, positions, total_real_particles,
                                    body_length, body_height, body_width, core_shell) /
        (hostCoreCombinedVectorFieldL2(particles, joule_fields.current_density_real, joule_fields.current_density_imag,
                                      positions, total_real_particles, body_length, body_height, body_width,
                                      core_shell) +
         TinyReal);
    return row;
}

inline void printHCurlHJObservableRow(const char *test_name, const AphiHCurlHJObservableRow &row)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(row.field_kind) << " dp=" << row.dp
              << " core_shell=" << row.core_shell << " B_err=" << row.b_error_vs_exact << " H_err=" << row.h_error_vs_exact
              << " H_to_B_err_ratio=" << row.h_to_b_error_ratio << " curlH_err=" << row.curl_h_error_vs_exact
              << " J_err=" << row.j_error_vs_exact << " curlH_to_J_err_ratio=" << row.curl_h_to_j_error_ratio
              << " discrete_curlH_to_J_l2_ratio=" << row.discrete_curl_h_to_j_l2_ratio << std::endl;
}

inline bool hCurlHJObservableLinear2DGateOk(const AphiHCurlHJObservableRow &row, Real max_b_error, Real max_h_error)
{
    if (row.b_error_vs_exact > max_b_error || row.h_error_vs_exact > max_h_error)
    {
        return false;
    }
    // H = nu*B is exact scaling; tiny B/H rel errors can make ratio deviate slightly from 1.
    const Real ratio_deviation = std::abs(row.h_to_b_error_ratio - 1.0);
    return ratio_deviation <= 0.15 || (row.b_error_vs_exact < 1.0e-5 && row.h_error_vs_exact < 1.0e-5);
}

inline bool hCurlHJObservableOscillatoryGateOk(const AphiHCurlHJObservableRow &row, Real max_b_error, Real max_h_error)
{
    return row.b_error_vs_exact <= max_b_error && row.h_error_vs_exact <= max_h_error &&
           std::abs(row.h_to_b_error_ratio - 1.0) <= 0.05;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_H_CURL_H_J_OBSERVABLE_DIAGNOSTIC_HELPERS_H
