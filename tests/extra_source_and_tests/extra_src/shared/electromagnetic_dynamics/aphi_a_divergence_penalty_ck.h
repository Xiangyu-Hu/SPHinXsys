#ifndef APHI_A_DIVERGENCE_PENALTY_CK_H
#define APHI_A_DIVERGENCE_PENALTY_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline std::string aphiDivAFieldName(const std::string &a_base_name) { return a_base_name + "DivA"; }

inline std::string aphiGradDivAFieldName(const std::string &a_base_name) { return aphiDivAFieldName(a_base_name) + "Gradient"; }

/** Shared scratch fields for penalty apply (one matvec at a time; not tied to input block prefix). */
struct AphiADivergencePenaltyScratchNames
{
    std::string div_a_real = "PenaltyScratchDivAReal";
    std::string div_a_imag = "PenaltyScratchDivAImag";
    std::string grad_div_a_real = "PenaltyScratchGradDivAReal";
    std::string grad_div_a_imag = "PenaltyScratchGradDivAImag";
};

inline const AphiADivergencePenaltyScratchNames &aphiDefaultADivergencePenaltyScratchNames()
{
    static const AphiADivergencePenaltyScratchNames names{};
    return names;
}

/** Stage 10.7: Coulomb gauge regularization on LhsA: lhs_a -= penalty * grad(div A). */
class AphiGradDivAPenaltyCK : public LocalDynamics
{
  public:
    explicit AphiGradDivAPenaltyCK(SPHBody &sph_body, const std::string &grad_div_a_real_name,
                                   const std::string &grad_div_a_imag_name, const AphiBlockNames &output_block,
                                   Real a_divergence_penalty);
    virtual ~AphiGradDivAPenaltyCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real a_divergence_penalty_;
        Vecd *grad_div_a_real_;
        Vecd *grad_div_a_imag_;
        Vecd *out_a_real_;
        Vecd *out_a_imag_;
    };

  protected:
    Real a_divergence_penalty_;
    DiscreteVariable<Vecd> *dv_grad_div_a_real_;
    DiscreteVariable<Vecd> *dv_grad_div_a_imag_;
    DiscreteVariable<Vecd> *dv_out_a_real_;
    DiscreteVariable<Vecd> *dv_out_a_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_A_DIVERGENCE_PENALTY_CK_H
