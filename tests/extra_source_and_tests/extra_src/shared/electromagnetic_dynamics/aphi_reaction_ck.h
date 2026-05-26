#ifndef APHI_REACTION_CK_H
#define APHI_REACTION_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

/** Local A-equation reaction: LhsAReal += -omega*sigma*AIm, LhsAImag += +omega*sigma*ARe. */
class AphiReactionCK : public LocalDynamics
{
  public:
    explicit AphiReactionCK(SPHBody &sph_body, Real omega, const AphiVariableNames &variable_names);
    virtual ~AphiReactionCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *input_a_real_;
        Vecd *input_a_imag_;
        Vecd *output_a_real_;
        Vecd *output_a_imag_;
        Real *sigma_;
        Real omega_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_input_a_real_;
    DiscreteVariable<Vecd> *dv_input_a_imag_;
    DiscreteVariable<Vecd> *dv_output_a_real_;
    DiscreteVariable<Vecd> *dv_output_a_imag_;
    DiscreteVariable<Real> *dv_sigma_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_REACTION_CK_H
