#ifndef APHI_CURL_A_DIAGNOSTIC_CK_HPP
#define APHI_CURL_A_DIAGNOSTIC_CK_HPP

#include "electromagnetic_dynamics/diagnostics/aphi_curl_a_diagnostic_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiVectorGradientCurlCK::AphiVectorGradientCurlCK(SPHBody &sph_body, const std::string &gradient_name,
                                                        const std::string &curl_name)
    : LocalDynamics(sph_body), dv_gradient_(particles_->template getVariableByName<Matd>(gradient_name)),
      dv_curl_(particles_->registerStateVariable<Vecd>(curl_name, ZeroData<Vecd>::value))
{
    particles_->addVariableToWrite<Vecd>(curl_name);
}

template <class ExecutionPolicy, class EncloserType>
inline AphiVectorGradientCurlCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : gradient_(encloser.dv_gradient_->DelegatedData(ex_policy)),
      curl_(encloser.dv_curl_->DelegatedData(ex_policy))
{
}

inline void AphiVectorGradientCurlCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Matd &grad_a = gradient_[index_i];
    curl_[index_i] = Vecd(grad_a(1, 2) - grad_a(2, 1), grad_a(2, 0) - grad_a(0, 2), grad_a(0, 1) - grad_a(1, 0));
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CURL_A_DIAGNOSTIC_CK_HPP
