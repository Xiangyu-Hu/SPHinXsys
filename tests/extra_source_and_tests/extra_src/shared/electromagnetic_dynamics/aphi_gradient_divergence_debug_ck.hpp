#ifndef APHI_GRADIENT_DIVERGENCE_DEBUG_CK_HPP
#define APHI_GRADIENT_DIVERGENCE_DEBUG_CK_HPP

#include "electromagnetic_dynamics/aphi_gradient_divergence_debug_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiVectorGradientDivergenceCK::AphiVectorGradientDivergenceCK(SPHBody &sph_body, const std::string &gradient_name,
                                                                      const std::string &divergence_name)
    : LocalDynamics(sph_body), dv_gradient_(particles_->template getVariableByName<Matd>(gradient_name)),
      dv_divergence_(particles_->registerStateVariable<Real>(divergence_name, Real(0)))
{
    particles_->addVariableToWrite<Real>(divergence_name);
}

template <class ExecutionPolicy, class EncloserType>
inline AphiVectorGradientDivergenceCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : gradient_(encloser.dv_gradient_->DelegatedData(ex_policy)),
      divergence_(encloser.dv_divergence_->DelegatedData(ex_policy))
{
}

inline void AphiVectorGradientDivergenceCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    divergence_[index_i] = gradient_[index_i].trace();
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GRADIENT_DIVERGENCE_DEBUG_CK_HPP
