#ifndef ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_HPP
#define ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_HPP

#include "electromagnetic_ophelie_biot_savart.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline ComputeOphelieCoilToGlassBiotSavartCK::ComputeOphelieCoilToGlassBiotSavartCK(
    SPHBody &glass_body, SPHBody &coil_body, const OphelieGlassFieldNames &glass_names,
    const OphelieCoilFieldNames &coil_names, const OphelieParameters &params)
    : LocalDynamics(glass_body), coeff_(params.mu0_ / (4.0 * Pi)),
      eps2_(params.softening_length_ * params.softening_length_), n_coil_(coil_body.getBaseParticles().TotalRealParticles()),
      dv_glass_pos_(particles_->template getVariableByName<Vecd>("Position")),
      dv_a_coil_real_(particles_->template getVariableByName<Vecd>(glass_names.a_coil_real)),
      dv_a_coil_imag_(particles_->template getVariableByName<Vecd>(glass_names.a_coil_imag)),
      dv_b_coil_real_(particles_->template getVariableByName<Vecd>(glass_names.b_coil_real)),
      dv_b_coil_imag_(particles_->template getVariableByName<Vecd>(glass_names.b_coil_imag)),
      dv_coil_pos_(coil_body.getBaseParticles().template getVariableByName<Vecd>("Position")),
      dv_j_src_real_(coil_body.getBaseParticles().template getVariableByName<Vecd>(coil_names.j_src_real)),
      dv_coil_vol_(coil_body.getBaseParticles().template getVariableByName<Real>("VolumetricMeasure"))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ComputeOphelieCoilToGlassBiotSavartCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                         EncloserType &encloser)
    : coeff_(encloser.coeff_), eps2_(encloser.eps2_), n_coil_(encloser.n_coil_),
      glass_pos_(encloser.dv_glass_pos_->DelegatedData(ex_policy)),
      a_coil_real_(encloser.dv_a_coil_real_->DelegatedData(ex_policy)),
      a_coil_imag_(encloser.dv_a_coil_imag_->DelegatedData(ex_policy)),
      b_coil_real_(encloser.dv_b_coil_real_->DelegatedData(ex_policy)),
      b_coil_imag_(encloser.dv_b_coil_imag_->DelegatedData(ex_policy)),
      coil_pos_(encloser.dv_coil_pos_->DelegatedData(ex_policy)),
      j_src_real_(encloser.dv_j_src_real_->DelegatedData(ex_policy)),
      coil_vol_(encloser.dv_coil_vol_->DelegatedData(ex_policy))
{
}

inline void ComputeOphelieCoilToGlassBiotSavartCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    Vecd a_sum = Vecd::Zero();
    Vecd b_sum = Vecd::Zero();
    const Vecd xi = glass_pos_[index_i];

    for (size_t j = 0; j < n_coil_; ++j)
    {
        const Vecd r = xi - coil_pos_[j];
        const Real r2 = r.squaredNorm() + eps2_;
        const Real inv_r = 1.0 / std::sqrt(r2);
        const Real inv_r3 = inv_r / r2;
        const Vecd jv = j_src_real_[j] * coil_vol_[j];
        a_sum += coeff_ * jv * inv_r;
        b_sum += coeff_ * jv.cross(r) * inv_r3;
    }

    a_coil_real_[index_i] = a_sum;
    a_coil_imag_[index_i] = Vecd::Zero();
    b_coil_real_[index_i] = b_sum;
    b_coil_imag_[index_i] = Vecd::Zero();
}

inline ComputeOphelieGlassSelfInducedBiotSavartCK::ComputeOphelieGlassSelfInducedBiotSavartCK(
    SPHBody &glass_body, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const std::string &j_real_field, const std::string &j_imag_field)
    : LocalDynamics(glass_body), coeff_(params.mu0_ / (4.0 * Pi)),
      eps2_(params.softening_length_ * params.softening_length_),
      n_glass_(particles_->TotalRealParticles()),
      dv_glass_pos_(particles_->template getVariableByName<Vecd>("Position")),
      dv_j_real_(particles_->template getVariableByName<Vecd>(j_real_field.empty() ? names.j_real : j_real_field)),
      dv_j_imag_(particles_->template getVariableByName<Vecd>(j_imag_field.empty() ? names.j_imag : j_imag_field)),
      dv_glass_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_a_ind_real_(particles_->template getVariableByName<Vecd>(names.a_ind_real)),
      dv_a_ind_imag_(particles_->template getVariableByName<Vecd>(names.a_ind_imag)),
      dv_b_ind_real_(particles_->template getVariableByName<Vecd>(names.b_ind_real)),
      dv_b_ind_imag_(particles_->template getVariableByName<Vecd>(names.b_ind_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ComputeOphelieGlassSelfInducedBiotSavartCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                               EncloserType &encloser)
    : coeff_(encloser.coeff_), eps2_(encloser.eps2_), n_glass_(encloser.n_glass_),
      glass_pos_(encloser.dv_glass_pos_->DelegatedData(ex_policy)),
      j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)),
      j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
      glass_vol_(encloser.dv_glass_vol_->DelegatedData(ex_policy)),
      a_ind_real_(encloser.dv_a_ind_real_->DelegatedData(ex_policy)),
      a_ind_imag_(encloser.dv_a_ind_imag_->DelegatedData(ex_policy)),
      b_ind_real_(encloser.dv_b_ind_real_->DelegatedData(ex_policy)),
      b_ind_imag_(encloser.dv_b_ind_imag_->DelegatedData(ex_policy))
{
}

inline void ComputeOphelieGlassSelfInducedBiotSavartCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    Vecd a_sum_r = Vecd::Zero();
    Vecd a_sum_i = Vecd::Zero();
    Vecd b_sum_r = Vecd::Zero();
    Vecd b_sum_i = Vecd::Zero();
    const Vecd xi = glass_pos_[index_i];

    for (size_t j = 0; j < n_glass_; ++j)
    {
        if (j == index_i)
        {
            continue;
        }
        const Vecd r = xi - glass_pos_[j];
        const Real r2 = r.squaredNorm() + eps2_;
        const Real inv_r = 1.0 / std::sqrt(r2);
        const Real inv_r3 = inv_r / r2;
        const Vecd moment_r = j_real_[j] * glass_vol_[j];
        const Vecd moment_i = j_imag_[j] * glass_vol_[j];
        a_sum_r += coeff_ * moment_r * inv_r;
        a_sum_i += coeff_ * moment_i * inv_r;
        b_sum_r += coeff_ * moment_r.cross(r) * inv_r3;
        b_sum_i += coeff_ * moment_i.cross(r) * inv_r3;
    }

    a_ind_real_[index_i] = a_sum_r;
    a_ind_imag_[index_i] = a_sum_i;
    b_ind_real_[index_i] = b_sum_r;
    b_ind_imag_[index_i] = b_sum_i;
}

inline CombineOphelieCoilAndInducedVectorPotentialCK::CombineOphelieCoilAndInducedVectorPotentialCK(
    SPHBody &sph_body, const OphelieGlassFieldNames &names)
    : LocalDynamics(sph_body),
      dv_a_coil_real_(particles_->template getVariableByName<Vecd>(names.a_coil_real)),
      dv_a_coil_imag_(particles_->template getVariableByName<Vecd>(names.a_coil_imag)),
      dv_b_coil_real_(particles_->template getVariableByName<Vecd>(names.b_coil_real)),
      dv_b_coil_imag_(particles_->template getVariableByName<Vecd>(names.b_coil_imag)),
      dv_a_ind_real_(particles_->template getVariableByName<Vecd>(names.a_ind_real)),
      dv_a_ind_imag_(particles_->template getVariableByName<Vecd>(names.a_ind_imag)),
      dv_b_ind_real_(particles_->template getVariableByName<Vecd>(names.b_ind_real)),
      dv_b_ind_imag_(particles_->template getVariableByName<Vecd>(names.b_ind_imag)),
      dv_a_src_real_(particles_->template getVariableByName<Vecd>(names.a_src_real)),
      dv_a_src_imag_(particles_->template getVariableByName<Vecd>(names.a_src_imag)),
      dv_b_src_real_(particles_->template getVariableByName<Vecd>(names.b_src_real)),
      dv_b_src_imag_(particles_->template getVariableByName<Vecd>(names.b_src_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline CombineOphelieCoilAndInducedVectorPotentialCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                                 EncloserType &encloser)
    : a_coil_real_(encloser.dv_a_coil_real_->DelegatedData(ex_policy)),
      a_coil_imag_(encloser.dv_a_coil_imag_->DelegatedData(ex_policy)),
      b_coil_real_(encloser.dv_b_coil_real_->DelegatedData(ex_policy)),
      b_coil_imag_(encloser.dv_b_coil_imag_->DelegatedData(ex_policy)),
      a_ind_real_(encloser.dv_a_ind_real_->DelegatedData(ex_policy)),
      a_ind_imag_(encloser.dv_a_ind_imag_->DelegatedData(ex_policy)),
      b_ind_real_(encloser.dv_b_ind_real_->DelegatedData(ex_policy)),
      b_ind_imag_(encloser.dv_b_ind_imag_->DelegatedData(ex_policy)),
      a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
      a_src_imag_(encloser.dv_a_src_imag_->DelegatedData(ex_policy)),
      b_src_real_(encloser.dv_b_src_real_->DelegatedData(ex_policy)),
      b_src_imag_(encloser.dv_b_src_imag_->DelegatedData(ex_policy))
{
}

inline void CombineOphelieCoilAndInducedVectorPotentialCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    a_src_real_[index_i] = a_coil_real_[index_i] + a_ind_real_[index_i];
    a_src_imag_[index_i] = a_coil_imag_[index_i] + a_ind_imag_[index_i];
    b_src_real_[index_i] = b_coil_real_[index_i] + b_ind_real_[index_i];
    b_src_imag_[index_i] = b_coil_imag_[index_i] + b_ind_imag_[index_i];
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_HPP
