#ifndef APHI_BLOCK_ZERO_CK_H
#define APHI_BLOCK_ZERO_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

/** Zero an A-phi block (solution, lhs, rhs, etc.). Stage 3B debug op only. */
class AphiZeroBlockCK : public LocalDynamics
{
  public:
    explicit AphiZeroBlockCK(SPHBody &sph_body, const AphiBlockNames &block_names);
    virtual ~AphiZeroBlockCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

/** Copy src block to dst block. Stage 3B discrete-RHS test helper only. */
class AphiCopyBlockCK : public LocalDynamics
{
  public:
    explicit AphiCopyBlockCK(SPHBody &sph_body, const AphiBlockNames &dst_names, const AphiBlockNames &src_names);
    virtual ~AphiCopyBlockCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *dst_a_real_;
        Vecd *dst_a_imag_;
        Real *dst_phi_real_;
        Real *dst_phi_imag_;
        Vecd *src_a_real_;
        Vecd *src_a_imag_;
        Real *src_phi_real_;
        Real *src_phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_dst_a_real_;
    DiscreteVariable<Vecd> *dv_dst_a_imag_;
    DiscreteVariable<Real> *dv_dst_phi_real_;
    DiscreteVariable<Real> *dv_dst_phi_imag_;
    DiscreteVariable<Vecd> *dv_src_a_real_;
    DiscreteVariable<Vecd> *dv_src_a_imag_;
    DiscreteVariable<Real> *dv_src_phi_real_;
    DiscreteVariable<Real> *dv_src_phi_imag_;
};

/** Residual = RHS - LHS into variable_names.residual. */
class AphiComputeResidualCK : public LocalDynamics
{
  public:
    explicit AphiComputeResidualCK(SPHBody &sph_body, const AphiVariableNames &variable_names);
    virtual ~AphiComputeResidualCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *residual_a_real_;
        Vecd *residual_a_imag_;
        Real *residual_phi_real_;
        Real *residual_phi_imag_;
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        Real *rhs_phi_real_;
        Real *rhs_phi_imag_;
        Vecd *lhs_a_real_;
        Vecd *lhs_a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_residual_a_real_;
    DiscreteVariable<Vecd> *dv_residual_a_imag_;
    DiscreteVariable<Real> *dv_residual_phi_real_;
    DiscreteVariable<Real> *dv_residual_phi_imag_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    DiscreteVariable<Real> *dv_rhs_phi_real_;
    DiscreteVariable<Real> *dv_rhs_phi_imag_;
    DiscreteVariable<Vecd> *dv_lhs_a_real_;
    DiscreteVariable<Vecd> *dv_lhs_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;
};

/** output = rhs_block - lhs_block. Used for true residual diagnostics. */
class AphiComputeBlockResidualCK : public LocalDynamics
{
  public:
    explicit AphiComputeBlockResidualCK(SPHBody &sph_body, const AphiBlockNames &output_names,
                                        const AphiBlockNames &rhs_names, const AphiBlockNames &lhs_names);
    virtual ~AphiComputeBlockResidualCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *output_a_real_;
        Vecd *output_a_imag_;
        Real *output_phi_real_;
        Real *output_phi_imag_;
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        Real *rhs_phi_real_;
        Real *rhs_phi_imag_;
        Vecd *lhs_a_real_;
        Vecd *lhs_a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_output_a_real_;
    DiscreteVariable<Vecd> *dv_output_a_imag_;
    DiscreteVariable<Real> *dv_output_phi_real_;
    DiscreteVariable<Real> *dv_output_phi_imag_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    DiscreteVariable<Real> *dv_rhs_phi_real_;
    DiscreteVariable<Real> *dv_rhs_phi_imag_;
    DiscreteVariable<Vecd> *dv_lhs_a_real_;
    DiscreteVariable<Vecd> *dv_lhs_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_ZERO_CK_H
