#ifndef APHI_GRADIENT_DIVERGENCE_DEBUG_CK_H
#define APHI_GRADIENT_DIVERGENCE_DEBUG_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace electromagnetics
{

class AphiVectorGradientDivergenceCK : public LocalDynamics
{
  public:
    explicit AphiVectorGradientDivergenceCK(SPHBody &sph_body, const std::string &gradient_name,
                                            const std::string &divergence_name);
    virtual ~AphiVectorGradientDivergenceCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Matd *gradient_;
        Real *divergence_;
    };

  protected:
    DiscreteVariable<Matd> *dv_gradient_;
    DiscreteVariable<Real> *dv_divergence_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GRADIENT_DIVERGENCE_DEBUG_CK_H
