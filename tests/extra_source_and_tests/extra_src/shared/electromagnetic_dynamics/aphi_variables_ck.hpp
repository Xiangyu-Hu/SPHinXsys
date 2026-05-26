#ifndef APHI_VARIABLES_CK_HPP
#define APHI_VARIABLES_CK_HPP

#include "electromagnetic_dynamics/aphi_variables_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline InitializeAphiVariablesCK::InitializeAphiVariablesCK(SPHBody &sph_body, Real default_sigma, Real default_nu,
                                                            const AphiVariableNames &variable_names)
    : LocalDynamics(sph_body), default_sigma_(default_sigma), default_nu_(default_nu), variable_names_(variable_names),
      dv_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.solution.a_real)),
      dv_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.solution.a_imag)),
      dv_phi_real_(particles_->registerStateVariable<Real>(variable_names_.solution.phi_real)),
      dv_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.solution.phi_imag)),
      dv_rhs_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.rhs.a_real)),
      dv_rhs_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.rhs.a_imag)),
      dv_rhs_phi_real_(particles_->registerStateVariable<Real>(variable_names_.rhs.phi_real)),
      dv_rhs_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.rhs.phi_imag)),
      dv_lhs_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.lhs.a_real)),
      dv_lhs_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.lhs.a_imag)),
      dv_lhs_phi_real_(particles_->registerStateVariable<Real>(variable_names_.lhs.phi_real)),
      dv_lhs_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.lhs.phi_imag)),
      dv_residual_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.residual.a_real)),
      dv_residual_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.residual.a_imag)),
      dv_residual_phi_real_(particles_->registerStateVariable<Real>(variable_names_.residual.phi_real)),
      dv_residual_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.residual.phi_imag)),
      dv_r_hat_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.r_hat.a_real)),
      dv_r_hat_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.r_hat.a_imag)),
      dv_r_hat_phi_real_(particles_->registerStateVariable<Real>(variable_names_.r_hat.phi_real)),
      dv_r_hat_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.r_hat.phi_imag)),
      dv_search_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.search.a_real)),
      dv_search_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.search.a_imag)),
      dv_search_phi_real_(particles_->registerStateVariable<Real>(variable_names_.search.phi_real)),
      dv_search_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.search.phi_imag)),
      dv_v_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.v.a_real)),
      dv_v_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.v.a_imag)),
      dv_v_phi_real_(particles_->registerStateVariable<Real>(variable_names_.v.phi_real)),
      dv_v_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.v.phi_imag)),
      dv_s_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.s.a_real)),
      dv_s_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.s.a_imag)),
      dv_s_phi_real_(particles_->registerStateVariable<Real>(variable_names_.s.phi_real)),
      dv_s_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.s.phi_imag)),
      dv_t_a_real_(particles_->registerStateVariable<Vecd>(variable_names_.t.a_real)),
      dv_t_a_imag_(particles_->registerStateVariable<Vecd>(variable_names_.t.a_imag)),
      dv_t_phi_real_(particles_->registerStateVariable<Real>(variable_names_.t.phi_real)),
      dv_t_phi_imag_(particles_->registerStateVariable<Real>(variable_names_.t.phi_imag)),
      dv_sigma_(particles_->registerStateVariable<Real>(variable_names_.material.sigma, default_sigma_)),
      dv_nu_(particles_->registerStateVariable<Real>(variable_names_.material.nu, default_nu_))
{
    addBlockVariablesToWrite(variable_names_.solution);
    addBlockVariablesToWrite(variable_names_.rhs);
    addBlockVariablesToWrite(variable_names_.lhs);
    addBlockVariablesToWrite(variable_names_.residual);
    addBlockVariablesToWrite(variable_names_.r_hat);
    addBlockVariablesToWrite(variable_names_.search);
    addBlockVariablesToWrite(variable_names_.v);
    addBlockVariablesToWrite(variable_names_.s);
    addBlockVariablesToWrite(variable_names_.t);
    particles_->addVariableToWrite<Real>(variable_names_.material.sigma);
    particles_->addVariableToWrite<Real>(variable_names_.material.nu);
}

inline void InitializeAphiVariablesCK::addBlockVariablesToWrite(const AphiBlockNames &block_names)
{
    particles_->addVariableToWrite<Vecd>(block_names.a_real);
    particles_->addVariableToWrite<Vecd>(block_names.a_imag);
    particles_->addVariableToWrite<Real>(block_names.phi_real);
    particles_->addVariableToWrite<Real>(block_names.phi_imag);
}

template <class ExecutionPolicy, class EncloserType>
inline InitializeAphiVariablesCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
      a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)),
      rhs_phi_real_(encloser.dv_rhs_phi_real_->DelegatedData(ex_policy)),
      rhs_phi_imag_(encloser.dv_rhs_phi_imag_->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy)),
      lhs_phi_real_(encloser.dv_lhs_phi_real_->DelegatedData(ex_policy)),
      lhs_phi_imag_(encloser.dv_lhs_phi_imag_->DelegatedData(ex_policy)),
      residual_a_real_(encloser.dv_residual_a_real_->DelegatedData(ex_policy)),
      residual_a_imag_(encloser.dv_residual_a_imag_->DelegatedData(ex_policy)),
      residual_phi_real_(encloser.dv_residual_phi_real_->DelegatedData(ex_policy)),
      residual_phi_imag_(encloser.dv_residual_phi_imag_->DelegatedData(ex_policy)),
      r_hat_a_real_(encloser.dv_r_hat_a_real_->DelegatedData(ex_policy)),
      r_hat_a_imag_(encloser.dv_r_hat_a_imag_->DelegatedData(ex_policy)),
      r_hat_phi_real_(encloser.dv_r_hat_phi_real_->DelegatedData(ex_policy)),
      r_hat_phi_imag_(encloser.dv_r_hat_phi_imag_->DelegatedData(ex_policy)),
      search_a_real_(encloser.dv_search_a_real_->DelegatedData(ex_policy)),
      search_a_imag_(encloser.dv_search_a_imag_->DelegatedData(ex_policy)),
      search_phi_real_(encloser.dv_search_phi_real_->DelegatedData(ex_policy)),
      search_phi_imag_(encloser.dv_search_phi_imag_->DelegatedData(ex_policy)),
      v_a_real_(encloser.dv_v_a_real_->DelegatedData(ex_policy)),
      v_a_imag_(encloser.dv_v_a_imag_->DelegatedData(ex_policy)),
      v_phi_real_(encloser.dv_v_phi_real_->DelegatedData(ex_policy)),
      v_phi_imag_(encloser.dv_v_phi_imag_->DelegatedData(ex_policy)),
      s_a_real_(encloser.dv_s_a_real_->DelegatedData(ex_policy)),
      s_a_imag_(encloser.dv_s_a_imag_->DelegatedData(ex_policy)),
      s_phi_real_(encloser.dv_s_phi_real_->DelegatedData(ex_policy)),
      s_phi_imag_(encloser.dv_s_phi_imag_->DelegatedData(ex_policy)),
      t_a_real_(encloser.dv_t_a_real_->DelegatedData(ex_policy)),
      t_a_imag_(encloser.dv_t_a_imag_->DelegatedData(ex_policy)),
      t_phi_real_(encloser.dv_t_phi_real_->DelegatedData(ex_policy)),
      t_phi_imag_(encloser.dv_t_phi_imag_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      default_sigma_(encloser.default_sigma_),
      default_nu_(encloser.default_nu_)
{
}

inline void InitializeAphiVariablesCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd zero_vec = Vecd::Zero();
    a_real_[index_i] = zero_vec;
    a_imag_[index_i] = zero_vec;
    phi_real_[index_i] = 0.0;
    phi_imag_[index_i] = 0.0;

    rhs_a_real_[index_i] = zero_vec;
    rhs_a_imag_[index_i] = zero_vec;
    rhs_phi_real_[index_i] = 0.0;
    rhs_phi_imag_[index_i] = 0.0;

    lhs_a_real_[index_i] = zero_vec;
    lhs_a_imag_[index_i] = zero_vec;
    lhs_phi_real_[index_i] = 0.0;
    lhs_phi_imag_[index_i] = 0.0;

    residual_a_real_[index_i] = zero_vec;
    residual_a_imag_[index_i] = zero_vec;
    residual_phi_real_[index_i] = 0.0;
    residual_phi_imag_[index_i] = 0.0;

    r_hat_a_real_[index_i] = zero_vec;
    r_hat_a_imag_[index_i] = zero_vec;
    r_hat_phi_real_[index_i] = 0.0;
    r_hat_phi_imag_[index_i] = 0.0;

    search_a_real_[index_i] = zero_vec;
    search_a_imag_[index_i] = zero_vec;
    search_phi_real_[index_i] = 0.0;
    search_phi_imag_[index_i] = 0.0;

    v_a_real_[index_i] = zero_vec;
    v_a_imag_[index_i] = zero_vec;
    v_phi_real_[index_i] = 0.0;
    v_phi_imag_[index_i] = 0.0;

    s_a_real_[index_i] = zero_vec;
    s_a_imag_[index_i] = zero_vec;
    s_phi_real_[index_i] = 0.0;
    s_phi_imag_[index_i] = 0.0;

    t_a_real_[index_i] = zero_vec;
    t_a_imag_[index_i] = zero_vec;
    t_phi_real_[index_i] = 0.0;
    t_phi_imag_[index_i] = 0.0;

    sigma_[index_i] = default_sigma_;
    nu_[index_i] = default_nu_;
}

inline SetAphiMaterialPropertiesCK::SetAphiMaterialPropertiesCK(SPHBody &sph_body, Real sigma, Real nu,
                                                                const AphiMaterialNames &material_names)
    : LocalDynamics(sph_body),
      dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(particles_->template getVariableByName<Real>(material_names.nu)),
      sigma_value_(sigma),
      nu_value_(nu)
{
}

template <class ExecutionPolicy, class EncloserType>
inline SetAphiMaterialPropertiesCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      sigma_value_(encloser.sigma_value_),
      nu_value_(encloser.nu_value_)
{
}

inline void SetAphiMaterialPropertiesCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    sigma_[index_i] = sigma_value_;
    nu_[index_i] = nu_value_;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_VARIABLES_CK_HPP
