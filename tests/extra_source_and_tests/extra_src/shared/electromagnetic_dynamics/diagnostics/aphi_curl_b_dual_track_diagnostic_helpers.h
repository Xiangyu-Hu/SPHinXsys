#ifndef APHI_CURL_B_DUAL_TRACK_DIAGNOSTIC_HELPERS_H
#define APHI_CURL_B_DUAL_TRACK_DIAGNOSTIC_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_curl_a_diagnostic_ck.hpp"
#include "general_gradient.h"
#include "kernel_correction_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

enum class AphiBCurlDiagnosticMode
{
    BCorrectedGrad,
    PairwiseUncorrectedGrad
};

template <typename InnerRelationType>
inline void execBodyCurlBFromADiagnostic(SPHBody &body, InnerRelationType &inner, const AphiVariableNames &names,
                                         const std::string &b_real_name, const std::string &b_imag_name,
                                         AphiBCurlDiagnosticMode mode)
{
    const std::string grad_a_real_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? names.solution.a_real + "Gradient"
                                                        : names.solution.a_real + "GradientPairwise";
    const std::string grad_a_imag_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? names.solution.a_imag + "Gradient"
                                                        : names.solution.a_imag + "GradientPairwise";

    BaseParticles &particles = body.getBaseParticles();
    if (mode == AphiBCurlDiagnosticMode::PairwiseUncorrectedGrad)
    {
        particles.registerStateVariable<Matd>(grad_a_real_name, ZeroData<Matd>::value);
        particles.registerStateVariable<Matd>(grad_a_imag_name, ZeroData<Matd>::value);
        InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorGradientCK<Inner<>>> a_real_gradient(
            inner, names.solution.a_real, grad_a_real_name);
        InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorGradientCK<Inner<>>> a_imag_gradient(
            inner, names.solution.a_imag, grad_a_imag_name);
        a_real_gradient.exec();
        a_imag_gradient.exec();
    }
    else
    {
        InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(inner);
        InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(inner, names.solution.a_real);
        InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(inner, names.solution.a_imag);
        linear_correction_matrix.exec();
        a_real_gradient.exec();
        a_imag_gradient.exec();
    }

    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_real(body, grad_a_real_name, b_real_name);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_imag(body, grad_a_imag_name, b_imag_name);
    curl_a_real.exec();
    curl_a_imag.exec();
}

inline void execBodyCurlBFromADiagnostic(SPHBody &body, Inner<> &inner, const AphiVariableNames &names,
                                         const std::string &b_real_name, const std::string &b_imag_name,
                                         AphiBCurlDiagnosticMode mode)
{
    execBodyCurlBFromADiagnostic<Inner<>>(body, inner, names, b_real_name, b_imag_name, mode);
}

/** Adaptive air inner: Inner<Relation<AdaptiveBody>> with LinearGradient<Inner<Vecd, Relation<...>>>. */
template <typename AdaptiveAirInnerType, typename AdaptiveAirBodyType>
inline void execBodyCurlBFromAdaptiveInnerDiagnostic(SPHBody &body, AdaptiveAirInnerType &inner,
                                                     const AphiVariableNames &names, const std::string &b_real_name,
                                                     const std::string &b_imag_name, AphiBCurlDiagnosticMode mode)
{
    const std::string grad_a_real_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? names.solution.a_real + "Gradient"
                                                        : names.solution.a_real + "GradientPairwise";
    const std::string grad_a_imag_name =
        mode == AphiBCurlDiagnosticMode::BCorrectedGrad ? names.solution.a_imag + "Gradient"
                                                        : names.solution.a_imag + "GradientPairwise";
    (void)mode;
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearCorrectionMatrix<Inner<WithUpdate, Relation<AdaptiveAirBodyType>>>>
        linear_correction_matrix(inner);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd, Relation<AdaptiveAirBodyType>>>>
        a_real_gradient(inner, names.solution.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd, Relation<AdaptiveAirBodyType>>>>
        a_imag_gradient(inner, names.solution.a_imag);
    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_real(body, grad_a_real_name, b_real_name);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_imag(body, grad_a_imag_name, b_imag_name);
    curl_a_real.exec();
    curl_a_imag.exec();
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CURL_B_DUAL_TRACK_DIAGNOSTIC_HELPERS_H
