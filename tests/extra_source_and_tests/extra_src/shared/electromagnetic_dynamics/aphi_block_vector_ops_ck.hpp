#ifndef APHI_BLOCK_VECTOR_OPS_CK_HPP
#define APHI_BLOCK_VECTOR_OPS_CK_HPP

#include "electromagnetic_dynamics/aphi_block_vector_ops_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiBlockAXPYCK::AphiBlockAXPYCK(SPHBody &sph_body, const AphiBlockNames &dst_names, Real alpha,
                                        const AphiBlockNames &src_names)
    : LocalDynamics(sph_body), alpha_(alpha),
      dv_dst_a_real_(particles_->template getVariableByName<Vecd>(dst_names.a_real)),
      dv_dst_a_imag_(particles_->template getVariableByName<Vecd>(dst_names.a_imag)),
      dv_dst_phi_real_(particles_->template getVariableByName<Real>(dst_names.phi_real)),
      dv_dst_phi_imag_(particles_->template getVariableByName<Real>(dst_names.phi_imag)),
      dv_src_a_real_(particles_->template getVariableByName<Vecd>(src_names.a_real)),
      dv_src_a_imag_(particles_->template getVariableByName<Vecd>(src_names.a_imag)),
      dv_src_phi_real_(particles_->template getVariableByName<Real>(src_names.phi_real)),
      dv_src_phi_imag_(particles_->template getVariableByName<Real>(src_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockAXPYCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : alpha_(encloser.alpha_),
      dst_a_real_(encloser.dv_dst_a_real_->DelegatedData(ex_policy)),
      dst_a_imag_(encloser.dv_dst_a_imag_->DelegatedData(ex_policy)),
      dst_phi_real_(encloser.dv_dst_phi_real_->DelegatedData(ex_policy)),
      dst_phi_imag_(encloser.dv_dst_phi_imag_->DelegatedData(ex_policy)),
      src_a_real_(encloser.dv_src_a_real_->DelegatedData(ex_policy)),
      src_a_imag_(encloser.dv_src_a_imag_->DelegatedData(ex_policy)),
      src_phi_real_(encloser.dv_src_phi_real_->DelegatedData(ex_policy)),
      src_phi_imag_(encloser.dv_src_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockAXPYCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    dst_a_real_[index_i] += alpha_ * src_a_real_[index_i];
    dst_a_imag_[index_i] += alpha_ * src_a_imag_[index_i];
    dst_phi_real_[index_i] += alpha_ * src_phi_real_[index_i];
    dst_phi_imag_[index_i] += alpha_ * src_phi_imag_[index_i];
}

inline AphiBlockLinearCombinationCK::AphiBlockLinearCombinationCK(
    SPHBody &sph_body, const AphiBlockNames &dst_names, Real coeff_x, Real coeff_y, const AphiBlockNames &block_x,
    const AphiBlockNames &block_y)
    : LocalDynamics(sph_body), coeff_x_(coeff_x), coeff_y_(coeff_y),
      dv_dst_a_real_(particles_->template getVariableByName<Vecd>(dst_names.a_real)),
      dv_dst_a_imag_(particles_->template getVariableByName<Vecd>(dst_names.a_imag)),
      dv_dst_phi_real_(particles_->template getVariableByName<Real>(dst_names.phi_real)),
      dv_dst_phi_imag_(particles_->template getVariableByName<Real>(dst_names.phi_imag)),
      dv_x_a_real_(particles_->template getVariableByName<Vecd>(block_x.a_real)),
      dv_x_a_imag_(particles_->template getVariableByName<Vecd>(block_x.a_imag)),
      dv_x_phi_real_(particles_->template getVariableByName<Real>(block_x.phi_real)),
      dv_x_phi_imag_(particles_->template getVariableByName<Real>(block_x.phi_imag)),
      dv_y_a_real_(particles_->template getVariableByName<Vecd>(block_y.a_real)),
      dv_y_a_imag_(particles_->template getVariableByName<Vecd>(block_y.a_imag)),
      dv_y_phi_real_(particles_->template getVariableByName<Real>(block_y.phi_real)),
      dv_y_phi_imag_(particles_->template getVariableByName<Real>(block_y.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockLinearCombinationCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                EncloserType &encloser)
    : coeff_x_(encloser.coeff_x_), coeff_y_(encloser.coeff_y_),
      dst_a_real_(encloser.dv_dst_a_real_->DelegatedData(ex_policy)),
      dst_a_imag_(encloser.dv_dst_a_imag_->DelegatedData(ex_policy)),
      dst_phi_real_(encloser.dv_dst_phi_real_->DelegatedData(ex_policy)),
      dst_phi_imag_(encloser.dv_dst_phi_imag_->DelegatedData(ex_policy)),
      x_a_real_(encloser.dv_x_a_real_->DelegatedData(ex_policy)),
      x_a_imag_(encloser.dv_x_a_imag_->DelegatedData(ex_policy)),
      x_phi_real_(encloser.dv_x_phi_real_->DelegatedData(ex_policy)),
      x_phi_imag_(encloser.dv_x_phi_imag_->DelegatedData(ex_policy)),
      y_a_real_(encloser.dv_y_a_real_->DelegatedData(ex_policy)),
      y_a_imag_(encloser.dv_y_a_imag_->DelegatedData(ex_policy)),
      y_phi_real_(encloser.dv_y_phi_real_->DelegatedData(ex_policy)),
      y_phi_imag_(encloser.dv_y_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockLinearCombinationCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    dst_a_real_[index_i] = coeff_x_ * x_a_real_[index_i] + coeff_y_ * y_a_real_[index_i];
    dst_a_imag_[index_i] = coeff_x_ * x_a_imag_[index_i] + coeff_y_ * y_a_imag_[index_i];
    dst_phi_real_[index_i] = coeff_x_ * x_phi_real_[index_i] + coeff_y_ * y_phi_real_[index_i];
    dst_phi_imag_[index_i] = coeff_x_ * x_phi_imag_[index_i] + coeff_y_ * y_phi_imag_[index_i];
}

inline AphiBlockBiCGStabUpdateSolutionCK::AphiBlockBiCGStabUpdateSolutionCK(
    SPHBody &sph_body, const AphiBlockNames &solution_names, Real alpha, Real omega,
    const AphiBlockNames &search_names, const AphiBlockNames &s_names)
    : LocalDynamics(sph_body), alpha_(alpha), omega_(omega),
      dv_x_a_real_(particles_->template getVariableByName<Vecd>(solution_names.a_real)),
      dv_x_a_imag_(particles_->template getVariableByName<Vecd>(solution_names.a_imag)),
      dv_x_phi_real_(particles_->template getVariableByName<Real>(solution_names.phi_real)),
      dv_x_phi_imag_(particles_->template getVariableByName<Real>(solution_names.phi_imag)),
      dv_p_a_real_(particles_->template getVariableByName<Vecd>(search_names.a_real)),
      dv_p_a_imag_(particles_->template getVariableByName<Vecd>(search_names.a_imag)),
      dv_p_phi_real_(particles_->template getVariableByName<Real>(search_names.phi_real)),
      dv_p_phi_imag_(particles_->template getVariableByName<Real>(search_names.phi_imag)),
      dv_s_a_real_(particles_->template getVariableByName<Vecd>(s_names.a_real)),
      dv_s_a_imag_(particles_->template getVariableByName<Vecd>(s_names.a_imag)),
      dv_s_phi_real_(particles_->template getVariableByName<Real>(s_names.phi_real)),
      dv_s_phi_imag_(particles_->template getVariableByName<Real>(s_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockBiCGStabUpdateSolutionCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                     EncloserType &encloser)
    : alpha_(encloser.alpha_), omega_(encloser.omega_),
      x_a_real_(encloser.dv_x_a_real_->DelegatedData(ex_policy)),
      x_a_imag_(encloser.dv_x_a_imag_->DelegatedData(ex_policy)),
      x_phi_real_(encloser.dv_x_phi_real_->DelegatedData(ex_policy)),
      x_phi_imag_(encloser.dv_x_phi_imag_->DelegatedData(ex_policy)),
      p_a_real_(encloser.dv_p_a_real_->DelegatedData(ex_policy)),
      p_a_imag_(encloser.dv_p_a_imag_->DelegatedData(ex_policy)),
      p_phi_real_(encloser.dv_p_phi_real_->DelegatedData(ex_policy)),
      p_phi_imag_(encloser.dv_p_phi_imag_->DelegatedData(ex_policy)),
      s_a_real_(encloser.dv_s_a_real_->DelegatedData(ex_policy)),
      s_a_imag_(encloser.dv_s_a_imag_->DelegatedData(ex_policy)),
      s_phi_real_(encloser.dv_s_phi_real_->DelegatedData(ex_policy)),
      s_phi_imag_(encloser.dv_s_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockBiCGStabUpdateSolutionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    x_a_real_[index_i] += alpha_ * p_a_real_[index_i] + omega_ * s_a_real_[index_i];
    x_a_imag_[index_i] += alpha_ * p_a_imag_[index_i] + omega_ * s_a_imag_[index_i];
    x_phi_real_[index_i] += alpha_ * p_phi_real_[index_i] + omega_ * s_phi_real_[index_i];
    x_phi_imag_[index_i] += alpha_ * p_phi_imag_[index_i] + omega_ * s_phi_imag_[index_i];
}

inline AphiBlockBiCGStabUpdateResidualCK::AphiBlockBiCGStabUpdateResidualCK(
    SPHBody &sph_body, const AphiBlockNames &residual_names, Real omega, const AphiBlockNames &s_names,
    const AphiBlockNames &t_names)
    : LocalDynamics(sph_body), omega_(omega),
      dv_r_a_real_(particles_->template getVariableByName<Vecd>(residual_names.a_real)),
      dv_r_a_imag_(particles_->template getVariableByName<Vecd>(residual_names.a_imag)),
      dv_r_phi_real_(particles_->template getVariableByName<Real>(residual_names.phi_real)),
      dv_r_phi_imag_(particles_->template getVariableByName<Real>(residual_names.phi_imag)),
      dv_s_a_real_(particles_->template getVariableByName<Vecd>(s_names.a_real)),
      dv_s_a_imag_(particles_->template getVariableByName<Vecd>(s_names.a_imag)),
      dv_s_phi_real_(particles_->template getVariableByName<Real>(s_names.phi_real)),
      dv_s_phi_imag_(particles_->template getVariableByName<Real>(s_names.phi_imag)),
      dv_t_a_real_(particles_->template getVariableByName<Vecd>(t_names.a_real)),
      dv_t_a_imag_(particles_->template getVariableByName<Vecd>(t_names.a_imag)),
      dv_t_phi_real_(particles_->template getVariableByName<Real>(t_names.phi_real)),
      dv_t_phi_imag_(particles_->template getVariableByName<Real>(t_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockBiCGStabUpdateResidualCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                     EncloserType &encloser)
    : omega_(encloser.omega_),
      r_a_real_(encloser.dv_r_a_real_->DelegatedData(ex_policy)),
      r_a_imag_(encloser.dv_r_a_imag_->DelegatedData(ex_policy)),
      r_phi_real_(encloser.dv_r_phi_real_->DelegatedData(ex_policy)),
      r_phi_imag_(encloser.dv_r_phi_imag_->DelegatedData(ex_policy)),
      s_a_real_(encloser.dv_s_a_real_->DelegatedData(ex_policy)),
      s_a_imag_(encloser.dv_s_a_imag_->DelegatedData(ex_policy)),
      s_phi_real_(encloser.dv_s_phi_real_->DelegatedData(ex_policy)),
      s_phi_imag_(encloser.dv_s_phi_imag_->DelegatedData(ex_policy)),
      t_a_real_(encloser.dv_t_a_real_->DelegatedData(ex_policy)),
      t_a_imag_(encloser.dv_t_a_imag_->DelegatedData(ex_policy)),
      t_phi_real_(encloser.dv_t_phi_real_->DelegatedData(ex_policy)),
      t_phi_imag_(encloser.dv_t_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockBiCGStabUpdateResidualCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    r_a_real_[index_i] = s_a_real_[index_i] - omega_ * t_a_real_[index_i];
    r_a_imag_[index_i] = s_a_imag_[index_i] - omega_ * t_a_imag_[index_i];
    r_phi_real_[index_i] = s_phi_real_[index_i] - omega_ * t_phi_real_[index_i];
    r_phi_imag_[index_i] = s_phi_imag_[index_i] - omega_ * t_phi_imag_[index_i];
}

inline AphiBlockBiCGStabUpdateSearchCK::AphiBlockBiCGStabUpdateSearchCK(
    SPHBody &sph_body, const AphiBlockNames &search_names, Real beta, Real omega,
    const AphiBlockNames &residual_names, const AphiBlockNames &v_names)
    : LocalDynamics(sph_body), beta_(beta), omega_(omega),
      dv_p_a_real_(particles_->template getVariableByName<Vecd>(search_names.a_real)),
      dv_p_a_imag_(particles_->template getVariableByName<Vecd>(search_names.a_imag)),
      dv_p_phi_real_(particles_->template getVariableByName<Real>(search_names.phi_real)),
      dv_p_phi_imag_(particles_->template getVariableByName<Real>(search_names.phi_imag)),
      dv_r_a_real_(particles_->template getVariableByName<Vecd>(residual_names.a_real)),
      dv_r_a_imag_(particles_->template getVariableByName<Vecd>(residual_names.a_imag)),
      dv_r_phi_real_(particles_->template getVariableByName<Real>(residual_names.phi_real)),
      dv_r_phi_imag_(particles_->template getVariableByName<Real>(residual_names.phi_imag)),
      dv_v_a_real_(particles_->template getVariableByName<Vecd>(v_names.a_real)),
      dv_v_a_imag_(particles_->template getVariableByName<Vecd>(v_names.a_imag)),
      dv_v_phi_real_(particles_->template getVariableByName<Real>(v_names.phi_real)),
      dv_v_phi_imag_(particles_->template getVariableByName<Real>(v_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockBiCGStabUpdateSearchCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                   EncloserType &encloser)
    : beta_(encloser.beta_), omega_(encloser.omega_),
      p_a_real_(encloser.dv_p_a_real_->DelegatedData(ex_policy)),
      p_a_imag_(encloser.dv_p_a_imag_->DelegatedData(ex_policy)),
      p_phi_real_(encloser.dv_p_phi_real_->DelegatedData(ex_policy)),
      p_phi_imag_(encloser.dv_p_phi_imag_->DelegatedData(ex_policy)),
      r_a_real_(encloser.dv_r_a_real_->DelegatedData(ex_policy)),
      r_a_imag_(encloser.dv_r_a_imag_->DelegatedData(ex_policy)),
      r_phi_real_(encloser.dv_r_phi_real_->DelegatedData(ex_policy)),
      r_phi_imag_(encloser.dv_r_phi_imag_->DelegatedData(ex_policy)),
      v_a_real_(encloser.dv_v_a_real_->DelegatedData(ex_policy)),
      v_a_imag_(encloser.dv_v_a_imag_->DelegatedData(ex_policy)),
      v_phi_real_(encloser.dv_v_phi_real_->DelegatedData(ex_policy)),
      v_phi_imag_(encloser.dv_v_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockBiCGStabUpdateSearchCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    p_a_real_[index_i] = r_a_real_[index_i] + beta_ * (p_a_real_[index_i] - omega_ * v_a_real_[index_i]);
    p_a_imag_[index_i] = r_a_imag_[index_i] + beta_ * (p_a_imag_[index_i] - omega_ * v_a_imag_[index_i]);
    p_phi_real_[index_i] = r_phi_real_[index_i] + beta_ * (p_phi_real_[index_i] - omega_ * v_phi_real_[index_i]);
    p_phi_imag_[index_i] = r_phi_imag_[index_i] + beta_ * (p_phi_imag_[index_i] - omega_ * v_phi_imag_[index_i]);
}

inline AphiBlockCGUpdateSearchCK::AphiBlockCGUpdateSearchCK(SPHBody &sph_body, const AphiBlockNames &search_names,
                                                            Real beta, const AphiBlockNames &z_names)
    : LocalDynamics(sph_body), beta_(beta),
      dv_p_a_real_(particles_->template getVariableByName<Vecd>(search_names.a_real)),
      dv_p_a_imag_(particles_->template getVariableByName<Vecd>(search_names.a_imag)),
      dv_p_phi_real_(particles_->template getVariableByName<Real>(search_names.phi_real)),
      dv_p_phi_imag_(particles_->template getVariableByName<Real>(search_names.phi_imag)),
      dv_z_a_real_(particles_->template getVariableByName<Vecd>(z_names.a_real)),
      dv_z_a_imag_(particles_->template getVariableByName<Vecd>(z_names.a_imag)),
      dv_z_phi_real_(particles_->template getVariableByName<Real>(z_names.phi_real)),
      dv_z_phi_imag_(particles_->template getVariableByName<Real>(z_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockCGUpdateSearchCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : beta_(encloser.beta_), p_a_real_(encloser.dv_p_a_real_->DelegatedData(ex_policy)),
      p_a_imag_(encloser.dv_p_a_imag_->DelegatedData(ex_policy)),
      p_phi_real_(encloser.dv_p_phi_real_->DelegatedData(ex_policy)),
      p_phi_imag_(encloser.dv_p_phi_imag_->DelegatedData(ex_policy)),
      z_a_real_(encloser.dv_z_a_real_->DelegatedData(ex_policy)),
      z_a_imag_(encloser.dv_z_a_imag_->DelegatedData(ex_policy)),
      z_phi_real_(encloser.dv_z_phi_real_->DelegatedData(ex_policy)),
      z_phi_imag_(encloser.dv_z_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockCGUpdateSearchCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    p_a_real_[index_i] = z_a_real_[index_i] + beta_ * p_a_real_[index_i];
    p_a_imag_[index_i] = z_a_imag_[index_i] + beta_ * p_a_imag_[index_i];
    p_phi_real_[index_i] = z_phi_real_[index_i] + beta_ * p_phi_real_[index_i];
    p_phi_imag_[index_i] = z_phi_imag_[index_i] + beta_ * p_phi_imag_[index_i];
}

inline AphiBlockScaleCopyCK::AphiBlockScaleCopyCK(SPHBody &sph_body, const AphiBlockNames &dst_names, Real alpha,
                                                  const AphiBlockNames &src_names)
    : LocalDynamics(sph_body), alpha_(alpha),
      dv_dst_a_real_(particles_->template getVariableByName<Vecd>(dst_names.a_real)),
      dv_dst_a_imag_(particles_->template getVariableByName<Vecd>(dst_names.a_imag)),
      dv_dst_phi_real_(particles_->template getVariableByName<Real>(dst_names.phi_real)),
      dv_dst_phi_imag_(particles_->template getVariableByName<Real>(dst_names.phi_imag)),
      dv_src_a_real_(particles_->template getVariableByName<Vecd>(src_names.a_real)),
      dv_src_a_imag_(particles_->template getVariableByName<Vecd>(src_names.a_imag)),
      dv_src_phi_real_(particles_->template getVariableByName<Real>(src_names.phi_real)),
      dv_src_phi_imag_(particles_->template getVariableByName<Real>(src_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockScaleCopyCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : alpha_(encloser.alpha_),
      dst_a_real_(encloser.dv_dst_a_real_->DelegatedData(ex_policy)),
      dst_a_imag_(encloser.dv_dst_a_imag_->DelegatedData(ex_policy)),
      dst_phi_real_(encloser.dv_dst_phi_real_->DelegatedData(ex_policy)),
      dst_phi_imag_(encloser.dv_dst_phi_imag_->DelegatedData(ex_policy)),
      src_a_real_(encloser.dv_src_a_real_->DelegatedData(ex_policy)),
      src_a_imag_(encloser.dv_src_a_imag_->DelegatedData(ex_policy)),
      src_phi_real_(encloser.dv_src_phi_real_->DelegatedData(ex_policy)),
      src_phi_imag_(encloser.dv_src_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiBlockScaleCopyCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    dst_a_real_[index_i] = alpha_ * src_a_real_[index_i];
    dst_a_imag_[index_i] = alpha_ * src_a_imag_[index_i];
    dst_phi_real_[index_i] = alpha_ * src_phi_real_[index_i];
    dst_phi_imag_[index_i] = alpha_ * src_phi_imag_[index_i];
}

inline AphiBlockDotProductCK::AphiBlockDotProductCK(SPHBody &sph_body, const AphiBlockNames &block_x,
                                                    const AphiBlockNames &block_y)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_x_a_real_(particles_->template getVariableByName<Vecd>(block_x.a_real)),
      dv_x_a_imag_(particles_->template getVariableByName<Vecd>(block_x.a_imag)),
      dv_x_phi_real_(particles_->template getVariableByName<Real>(block_x.phi_real)),
      dv_x_phi_imag_(particles_->template getVariableByName<Real>(block_x.phi_imag)),
      dv_y_a_real_(particles_->template getVariableByName<Vecd>(block_y.a_real)),
      dv_y_a_imag_(particles_->template getVariableByName<Vecd>(block_y.a_imag)),
      dv_y_phi_real_(particles_->template getVariableByName<Real>(block_y.phi_real)),
      dv_y_phi_imag_(particles_->template getVariableByName<Real>(block_y.phi_imag))
{
    quantity_name_ = "AphiBlockDotProduct";
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockDotProductCK::ReduceKernel::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vol_(encloser.dv_vol_->DelegatedData(ex_policy)),
      x_a_real_(encloser.dv_x_a_real_->DelegatedData(ex_policy)),
      x_a_imag_(encloser.dv_x_a_imag_->DelegatedData(ex_policy)),
      x_phi_real_(encloser.dv_x_phi_real_->DelegatedData(ex_policy)),
      x_phi_imag_(encloser.dv_x_phi_imag_->DelegatedData(ex_policy)),
      y_a_real_(encloser.dv_y_a_real_->DelegatedData(ex_policy)),
      y_a_imag_(encloser.dv_y_a_imag_->DelegatedData(ex_policy)),
      y_phi_real_(encloser.dv_y_phi_real_->DelegatedData(ex_policy)),
      y_phi_imag_(encloser.dv_y_phi_imag_->DelegatedData(ex_policy))
{
}

inline Real AphiBlockDotProductCK::ReduceKernel::reduce(size_t index_i, Real dt)
{
    (void)dt;
    const Real vol_i = vol_[index_i];
    return vol_i * (x_a_real_[index_i].dot(y_a_real_[index_i]) + x_a_imag_[index_i].dot(y_a_imag_[index_i]) +
                    x_phi_real_[index_i] * y_phi_real_[index_i] + x_phi_imag_[index_i] * y_phi_imag_[index_i]);
}

inline AphiBlockNormSquaredCK::AphiBlockNormSquaredCK(SPHBody &sph_body, const AphiBlockNames &block_names)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
      dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
      dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
      dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
{
    quantity_name_ = "AphiBlockNormSquared";
}

template <class ExecutionPolicy, class EncloserType>
inline AphiBlockNormSquaredCK::ReduceKernel::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vol_(encloser.dv_vol_->DelegatedData(ex_policy)),
      a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
      a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
{
}

inline Real AphiBlockNormSquaredCK::ReduceKernel::reduce(size_t index_i, Real dt)
{
    (void)dt;
    const Real vol_i = vol_[index_i];
    return vol_i * (a_real_[index_i].squaredNorm() + a_imag_[index_i].squaredNorm() +
                    phi_real_[index_i] * phi_real_[index_i] + phi_imag_[index_i] * phi_imag_[index_i]);
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_VECTOR_OPS_CK_HPP
