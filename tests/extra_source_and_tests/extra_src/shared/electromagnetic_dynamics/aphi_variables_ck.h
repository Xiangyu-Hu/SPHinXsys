#ifndef APHI_VARIABLES_CK_H
#define APHI_VARIABLES_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{

class InitializeAphiVariablesCK : public LocalDynamics
{
  public:
    explicit InitializeAphiVariablesCK(SPHBody &sph_body, Real default_sigma = 0.0, Real default_nu = 1.0,
                                       const AphiVariableNames &variable_names = AphiVariableNames{});
    virtual ~InitializeAphiVariablesCK() = default;

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

        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        Real *rhs_phi_real_;
        Real *rhs_phi_imag_;

        Vecd *lhs_a_real_;
        Vecd *lhs_a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;

        Vecd *residual_a_real_;
        Vecd *residual_a_imag_;
        Real *residual_phi_real_;
        Real *residual_phi_imag_;

        Vecd *r_hat_a_real_;
        Vecd *r_hat_a_imag_;
        Real *r_hat_phi_real_;
        Real *r_hat_phi_imag_;

        Vecd *search_a_real_;
        Vecd *search_a_imag_;
        Real *search_phi_real_;
        Real *search_phi_imag_;

        Vecd *v_a_real_;
        Vecd *v_a_imag_;
        Real *v_phi_real_;
        Real *v_phi_imag_;

        Vecd *s_a_real_;
        Vecd *s_a_imag_;
        Real *s_phi_real_;
        Real *s_phi_imag_;

        Vecd *t_a_real_;
        Vecd *t_a_imag_;
        Real *t_phi_real_;
        Real *t_phi_imag_;

        Vecd *v_old_a_real_;
        Vecd *v_old_a_imag_;
        Real *v_old_phi_real_;
        Real *v_old_phi_imag_;

        Vecd *z_a_real_;
        Vecd *z_a_imag_;
        Real *z_phi_real_;
        Real *z_phi_imag_;

        Vecd *y_a_real_;
        Vecd *y_a_imag_;
        Real *y_phi_real_;
        Real *y_phi_imag_;

        Vecd *true_residual_a_real_;
        Vecd *true_residual_a_imag_;
        Real *true_residual_phi_real_;
        Real *true_residual_phi_imag_;

        Real *sigma_;
        Real *nu_;
        Real default_sigma_;
        Real default_nu_;
    };

  protected:
    void addBlockVariablesToWrite(const AphiBlockNames &block_names);

    Real default_sigma_;
    Real default_nu_;
    AphiVariableNames variable_names_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;

    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    DiscreteVariable<Real> *dv_rhs_phi_real_;
    DiscreteVariable<Real> *dv_rhs_phi_imag_;

    DiscreteVariable<Vecd> *dv_lhs_a_real_;
    DiscreteVariable<Vecd> *dv_lhs_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;

    DiscreteVariable<Vecd> *dv_residual_a_real_;
    DiscreteVariable<Vecd> *dv_residual_a_imag_;
    DiscreteVariable<Real> *dv_residual_phi_real_;
    DiscreteVariable<Real> *dv_residual_phi_imag_;

    DiscreteVariable<Vecd> *dv_r_hat_a_real_;
    DiscreteVariable<Vecd> *dv_r_hat_a_imag_;
    DiscreteVariable<Real> *dv_r_hat_phi_real_;
    DiscreteVariable<Real> *dv_r_hat_phi_imag_;

    DiscreteVariable<Vecd> *dv_search_a_real_;
    DiscreteVariable<Vecd> *dv_search_a_imag_;
    DiscreteVariable<Real> *dv_search_phi_real_;
    DiscreteVariable<Real> *dv_search_phi_imag_;

    DiscreteVariable<Vecd> *dv_v_a_real_;
    DiscreteVariable<Vecd> *dv_v_a_imag_;
    DiscreteVariable<Real> *dv_v_phi_real_;
    DiscreteVariable<Real> *dv_v_phi_imag_;

    DiscreteVariable<Vecd> *dv_s_a_real_;
    DiscreteVariable<Vecd> *dv_s_a_imag_;
    DiscreteVariable<Real> *dv_s_phi_real_;
    DiscreteVariable<Real> *dv_s_phi_imag_;

    DiscreteVariable<Vecd> *dv_t_a_real_;
    DiscreteVariable<Vecd> *dv_t_a_imag_;
    DiscreteVariable<Real> *dv_t_phi_real_;
    DiscreteVariable<Real> *dv_t_phi_imag_;

    DiscreteVariable<Vecd> *dv_v_old_a_real_;
    DiscreteVariable<Vecd> *dv_v_old_a_imag_;
    DiscreteVariable<Real> *dv_v_old_phi_real_;
    DiscreteVariable<Real> *dv_v_old_phi_imag_;

    DiscreteVariable<Vecd> *dv_z_a_real_;
    DiscreteVariable<Vecd> *dv_z_a_imag_;
    DiscreteVariable<Real> *dv_z_phi_real_;
    DiscreteVariable<Real> *dv_z_phi_imag_;

    DiscreteVariable<Vecd> *dv_y_a_real_;
    DiscreteVariable<Vecd> *dv_y_a_imag_;
    DiscreteVariable<Real> *dv_y_phi_real_;
    DiscreteVariable<Real> *dv_y_phi_imag_;

    DiscreteVariable<Vecd> *dv_true_residual_a_real_;
    DiscreteVariable<Vecd> *dv_true_residual_a_imag_;
    DiscreteVariable<Real> *dv_true_residual_phi_real_;
    DiscreteVariable<Real> *dv_true_residual_phi_imag_;

    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
};

class SetAphiMaterialPropertiesCK : public LocalDynamics
{
  public:
    explicit SetAphiMaterialPropertiesCK(SPHBody &sph_body, Real sigma, Real nu,
                                         const AphiMaterialNames &material_names = AphiMaterialNames{});
    virtual ~SetAphiMaterialPropertiesCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *sigma_;
        Real *nu_;
        Real sigma_value_;
        Real nu_value_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    Real sigma_value_;
    Real nu_value_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_VARIABLES_CK_H
