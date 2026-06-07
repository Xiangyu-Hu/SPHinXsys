#ifndef APHI_PHI_GAUGE_PENALTY_CK_H
#define APHI_PHI_GAUGE_PENALTY_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

/** Stage 7: local phi gauge penalty on LhsPhi: lhs_phi += penalty * input_phi. */
class AphiPhiGaugePenaltyCK : public LocalDynamics
{
  public:
    explicit AphiPhiGaugePenaltyCK(SPHBody &sph_body, const AphiBlockNames &input_block, const AphiBlockNames &output_block,
                                   Real phi_gauge_penalty);
    virtual ~AphiPhiGaugePenaltyCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real phi_gauge_penalty_;
        Real *in_phi_real_;
        Real *in_phi_imag_;
        Real *out_phi_real_;
        Real *out_phi_imag_;
    };

  protected:
    Real phi_gauge_penalty_;
    DiscreteVariable<Real> *dv_in_phi_real_;
    DiscreteVariable<Real> *dv_in_phi_imag_;
    DiscreteVariable<Real> *dv_out_phi_real_;
    DiscreteVariable<Real> *dv_out_phi_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PHI_GAUGE_PENALTY_CK_H
