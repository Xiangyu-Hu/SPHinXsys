#ifndef APHI_BLOCK_VECTOR_OPS_CK_H
#define APHI_BLOCK_VECTOR_OPS_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

/** dst += alpha * src on all four A-phi block components. */
class AphiBlockAXPYCK : public LocalDynamics
{
  public:
    explicit AphiBlockAXPYCK(SPHBody &sph_body, const AphiBlockNames &dst_names, Real alpha,
                             const AphiBlockNames &src_names);
    virtual ~AphiBlockAXPYCK() = default;

    void setAlpha(Real alpha) { alpha_ = alpha; }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real alpha_;
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
    Real alpha_;
    DiscreteVariable<Vecd> *dv_dst_a_real_;
    DiscreteVariable<Vecd> *dv_dst_a_imag_;
    DiscreteVariable<Real> *dv_dst_phi_real_;
    DiscreteVariable<Real> *dv_dst_phi_imag_;
    DiscreteVariable<Vecd> *dv_src_a_real_;
    DiscreteVariable<Vecd> *dv_src_a_imag_;
    DiscreteVariable<Real> *dv_src_phi_real_;
    DiscreteVariable<Real> *dv_src_phi_imag_;
};

/** dst = coeff_x * block_x + coeff_y * block_y. */
class AphiBlockLinearCombinationCK : public LocalDynamics
{
  public:
    explicit AphiBlockLinearCombinationCK(SPHBody &sph_body, const AphiBlockNames &dst_names, Real coeff_x, Real coeff_y,
                                          const AphiBlockNames &block_x, const AphiBlockNames &block_y);
    virtual ~AphiBlockLinearCombinationCK() = default;

    void setCoefficients(Real coeff_x, Real coeff_y)
    {
        coeff_x_ = coeff_x;
        coeff_y_ = coeff_y;
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real coeff_x_;
        Real coeff_y_;
        Vecd *dst_a_real_;
        Vecd *dst_a_imag_;
        Real *dst_phi_real_;
        Real *dst_phi_imag_;
        Vecd *x_a_real_;
        Vecd *x_a_imag_;
        Real *x_phi_real_;
        Real *x_phi_imag_;
        Vecd *y_a_real_;
        Vecd *y_a_imag_;
        Real *y_phi_real_;
        Real *y_phi_imag_;
    };

  protected:
    Real coeff_x_;
    Real coeff_y_;
    DiscreteVariable<Vecd> *dv_dst_a_real_;
    DiscreteVariable<Vecd> *dv_dst_a_imag_;
    DiscreteVariable<Real> *dv_dst_phi_real_;
    DiscreteVariable<Real> *dv_dst_phi_imag_;
    DiscreteVariable<Vecd> *dv_x_a_real_;
    DiscreteVariable<Vecd> *dv_x_a_imag_;
    DiscreteVariable<Real> *dv_x_phi_real_;
    DiscreteVariable<Real> *dv_x_phi_imag_;
    DiscreteVariable<Vecd> *dv_y_a_real_;
    DiscreteVariable<Vecd> *dv_y_a_imag_;
    DiscreteVariable<Real> *dv_y_phi_real_;
    DiscreteVariable<Real> *dv_y_phi_imag_;
};

/** X += alpha * P + omega * S. */
class AphiBlockBiCGStabUpdateSolutionCK : public LocalDynamics
{
  public:
    explicit AphiBlockBiCGStabUpdateSolutionCK(SPHBody &sph_body, const AphiBlockNames &solution_names, Real alpha,
                                                 Real omega, const AphiBlockNames &search_names,
                                                 const AphiBlockNames &s_names);
    virtual ~AphiBlockBiCGStabUpdateSolutionCK() = default;

    void setScalars(Real alpha, Real omega)
    {
        alpha_ = alpha;
        omega_ = omega;
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real alpha_;
        Real omega_;
        Vecd *x_a_real_;
        Vecd *x_a_imag_;
        Real *x_phi_real_;
        Real *x_phi_imag_;
        Vecd *p_a_real_;
        Vecd *p_a_imag_;
        Real *p_phi_real_;
        Real *p_phi_imag_;
        Vecd *s_a_real_;
        Vecd *s_a_imag_;
        Real *s_phi_real_;
        Real *s_phi_imag_;
    };

  protected:
    Real alpha_;
    Real omega_;
    DiscreteVariable<Vecd> *dv_x_a_real_;
    DiscreteVariable<Vecd> *dv_x_a_imag_;
    DiscreteVariable<Real> *dv_x_phi_real_;
    DiscreteVariable<Real> *dv_x_phi_imag_;
    DiscreteVariable<Vecd> *dv_p_a_real_;
    DiscreteVariable<Vecd> *dv_p_a_imag_;
    DiscreteVariable<Real> *dv_p_phi_real_;
    DiscreteVariable<Real> *dv_p_phi_imag_;
    DiscreteVariable<Vecd> *dv_s_a_real_;
    DiscreteVariable<Vecd> *dv_s_a_imag_;
    DiscreteVariable<Real> *dv_s_phi_real_;
    DiscreteVariable<Real> *dv_s_phi_imag_;
};

/** R = S - omega * T. */
class AphiBlockBiCGStabUpdateResidualCK : public LocalDynamics
{
  public:
    explicit AphiBlockBiCGStabUpdateResidualCK(SPHBody &sph_body, const AphiBlockNames &residual_names, Real omega,
                                               const AphiBlockNames &s_names, const AphiBlockNames &t_names);
    virtual ~AphiBlockBiCGStabUpdateResidualCK() = default;

    void setOmega(Real omega) { omega_ = omega; }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real omega_;
        Vecd *r_a_real_;
        Vecd *r_a_imag_;
        Real *r_phi_real_;
        Real *r_phi_imag_;
        Vecd *s_a_real_;
        Vecd *s_a_imag_;
        Real *s_phi_real_;
        Real *s_phi_imag_;
        Vecd *t_a_real_;
        Vecd *t_a_imag_;
        Real *t_phi_real_;
        Real *t_phi_imag_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_r_a_real_;
    DiscreteVariable<Vecd> *dv_r_a_imag_;
    DiscreteVariable<Real> *dv_r_phi_real_;
    DiscreteVariable<Real> *dv_r_phi_imag_;
    DiscreteVariable<Vecd> *dv_s_a_real_;
    DiscreteVariable<Vecd> *dv_s_a_imag_;
    DiscreteVariable<Real> *dv_s_phi_real_;
    DiscreteVariable<Real> *dv_s_phi_imag_;
    DiscreteVariable<Vecd> *dv_t_a_real_;
    DiscreteVariable<Vecd> *dv_t_a_imag_;
    DiscreteVariable<Real> *dv_t_phi_real_;
    DiscreteVariable<Real> *dv_t_phi_imag_;
};

/** P = R + beta * (P - omega * V). */
class AphiBlockBiCGStabUpdateSearchCK : public LocalDynamics
{
  public:
    explicit AphiBlockBiCGStabUpdateSearchCK(SPHBody &sph_body, const AphiBlockNames &search_names, Real beta,
                                             Real omega, const AphiBlockNames &residual_names,
                                             const AphiBlockNames &v_names);
    virtual ~AphiBlockBiCGStabUpdateSearchCK() = default;

    void setScalars(Real beta, Real omega)
    {
        beta_ = beta;
        omega_ = omega;
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real beta_;
        Real omega_;
        Vecd *p_a_real_;
        Vecd *p_a_imag_;
        Real *p_phi_real_;
        Real *p_phi_imag_;
        Vecd *r_a_real_;
        Vecd *r_a_imag_;
        Real *r_phi_real_;
        Real *r_phi_imag_;
        Vecd *v_a_real_;
        Vecd *v_a_imag_;
        Real *v_phi_real_;
        Real *v_phi_imag_;
    };

  protected:
    Real beta_;
    Real omega_;
    DiscreteVariable<Vecd> *dv_p_a_real_;
    DiscreteVariable<Vecd> *dv_p_a_imag_;
    DiscreteVariable<Real> *dv_p_phi_real_;
    DiscreteVariable<Real> *dv_p_phi_imag_;
    DiscreteVariable<Vecd> *dv_r_a_real_;
    DiscreteVariable<Vecd> *dv_r_a_imag_;
    DiscreteVariable<Real> *dv_r_phi_real_;
    DiscreteVariable<Real> *dv_r_phi_imag_;
    DiscreteVariable<Vecd> *dv_v_a_real_;
    DiscreteVariable<Vecd> *dv_v_a_imag_;
    DiscreteVariable<Real> *dv_v_phi_real_;
    DiscreteVariable<Real> *dv_v_phi_imag_;
};

/** P = Z + beta * P for preconditioned CG. */
class AphiBlockCGUpdateSearchCK : public LocalDynamics
{
  public:
    explicit AphiBlockCGUpdateSearchCK(SPHBody &sph_body, const AphiBlockNames &search_names, Real beta,
                                       const AphiBlockNames &z_names);
    virtual ~AphiBlockCGUpdateSearchCK() = default;

    void setBeta(Real beta) { beta_ = beta; }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real beta_;
        Vecd *p_a_real_;
        Vecd *p_a_imag_;
        Real *p_phi_real_;
        Real *p_phi_imag_;
        Vecd *z_a_real_;
        Vecd *z_a_imag_;
        Real *z_phi_real_;
        Real *z_phi_imag_;
    };

  protected:
    Real beta_;
    DiscreteVariable<Vecd> *dv_p_a_real_;
    DiscreteVariable<Vecd> *dv_p_a_imag_;
    DiscreteVariable<Real> *dv_p_phi_real_;
    DiscreteVariable<Real> *dv_p_phi_imag_;
    DiscreteVariable<Vecd> *dv_z_a_real_;
    DiscreteVariable<Vecd> *dv_z_a_imag_;
    DiscreteVariable<Real> *dv_z_phi_real_;
    DiscreteVariable<Real> *dv_z_phi_imag_;
};

/** dst = alpha * src (overwrite, not accumulate). */
class AphiBlockScaleCopyCK : public LocalDynamics
{
  public:
    explicit AphiBlockScaleCopyCK(SPHBody &sph_body, const AphiBlockNames &dst_names, Real alpha,
                                  const AphiBlockNames &src_names);
    virtual ~AphiBlockScaleCopyCK() = default;

    void setAlpha(Real alpha) { alpha_ = alpha; }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real alpha_;
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
    Real alpha_;
    DiscreteVariable<Vecd> *dv_dst_a_real_;
    DiscreteVariable<Vecd> *dv_dst_a_imag_;
    DiscreteVariable<Real> *dv_dst_phi_real_;
    DiscreteVariable<Real> *dv_dst_phi_imag_;
    DiscreteVariable<Vecd> *dv_src_a_real_;
    DiscreteVariable<Vecd> *dv_src_a_imag_;
    DiscreteVariable<Real> *dv_src_phi_real_;
    DiscreteVariable<Real> *dv_src_phi_imag_;
};

/** Vol-weighted inner product <X,Y>. */
class AphiBlockDotProductCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    explicit AphiBlockDotProductCK(SPHBody &sph_body, const AphiBlockNames &block_x, const AphiBlockNames &block_y);
    virtual ~AphiBlockDotProductCK() = default;

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        Real reduce(size_t index_i, Real dt = 0.0);

      protected:
        Real *vol_;
        Vecd *x_a_real_;
        Vecd *x_a_imag_;
        Real *x_phi_real_;
        Real *x_phi_imag_;
        Vecd *y_a_real_;
        Vecd *y_a_imag_;
        Real *y_phi_real_;
        Real *y_phi_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_vol_;
    DiscreteVariable<Vecd> *dv_x_a_real_;
    DiscreteVariable<Vecd> *dv_x_a_imag_;
    DiscreteVariable<Real> *dv_x_phi_real_;
    DiscreteVariable<Real> *dv_x_phi_imag_;
    DiscreteVariable<Vecd> *dv_y_a_real_;
    DiscreteVariable<Vecd> *dv_y_a_imag_;
    DiscreteVariable<Real> *dv_y_phi_real_;
    DiscreteVariable<Real> *dv_y_phi_imag_;
};

/** Vol-weighted squared norm <X,X>. Norm is sqrt(exec result). */
class AphiBlockNormSquaredCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    explicit AphiBlockNormSquaredCK(SPHBody &sph_body, const AphiBlockNames &block_names);
    virtual ~AphiBlockNormSquaredCK() = default;

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        Real reduce(size_t index_i, Real dt = 0.0);

      protected:
        Real *vol_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_vol_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_VECTOR_OPS_CK_H
