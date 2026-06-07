#ifndef APHI_CURL_A_DIAGNOSTIC_CK_H
#define APHI_CURL_A_DIAGNOSTIC_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace electromagnetics
{

/** B = curl A from SPH gradient tensor: curl_m = gradA[n,m] antisymmetric part. */
class AphiVectorGradientCurlCK : public LocalDynamics
{
  public:
    explicit AphiVectorGradientCurlCK(SPHBody &sph_body, const std::string &gradient_name, const std::string &curl_name);
    virtual ~AphiVectorGradientCurlCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Matd *gradient_;
        Vecd *curl_;
    };

  protected:
    DiscreteVariable<Matd> *dv_gradient_;
    DiscreteVariable<Vecd> *dv_curl_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CURL_A_DIAGNOSTIC_CK_H
