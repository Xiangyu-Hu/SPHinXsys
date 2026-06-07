#ifndef APHI_BLOCK_JACOBI_PRECONDITIONER_CK_HPP
#define APHI_BLOCK_JACOBI_PRECONDITIONER_CK_HPP

#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

namespace detail
{

inline void aphiZeroMatrix8(Real matrix[8][8])
{
    for (UnsignedInt row = 0; row < 8; ++row)
    {
        for (UnsignedInt col = 0; col < 8; ++col)
        {
            matrix[row][col] = Real(0);
        }
    }
}

inline bool aphiSolve8x8(Real matrix[8][8], const Real rhs[8], Real solution[8])
{
    Real augmented[8][9];
    for (UnsignedInt row = 0; row < 8; ++row)
    {
        for (UnsignedInt col = 0; col < 8; ++col)
        {
            augmented[row][col] = matrix[row][col];
        }
        augmented[row][8] = rhs[row];
    }

    for (UnsignedInt pivot_col = 0; pivot_col < 8; ++pivot_col)
    {
        UnsignedInt pivot_row = pivot_col;
        Real pivot_abs_max = std::abs(augmented[pivot_col][pivot_col]);
        for (UnsignedInt row = pivot_col + 1; row < 8; ++row)
        {
            const Real candidate = std::abs(augmented[row][pivot_col]);
            if (candidate > pivot_abs_max)
            {
                pivot_abs_max = candidate;
                pivot_row = row;
            }
        }

        if (pivot_abs_max <= TinyReal)
        {
            return false;
        }

        if (pivot_row != pivot_col)
        {
            for (UnsignedInt col = pivot_col; col < 9; ++col)
            {
                const Real temp = augmented[pivot_col][col];
                augmented[pivot_col][col] = augmented[pivot_row][col];
                augmented[pivot_row][col] = temp;
            }
        }

        const Real pivot_value = augmented[pivot_col][pivot_col];
        const Real inv_pivot = Real(1) / pivot_value;
        for (UnsignedInt col = pivot_col; col < 9; ++col)
        {
            augmented[pivot_col][col] *= inv_pivot;
        }

        for (UnsignedInt row = 0; row < 8; ++row)
        {
            if (row == pivot_col)
            {
                continue;
            }
            const Real factor = augmented[row][pivot_col];
            if (std::abs(factor) <= TinyReal)
            {
                continue;
            }
            for (UnsignedInt col = pivot_col; col < 9; ++col)
            {
                augmented[row][col] -= factor * augmented[pivot_col][col];
            }
        }
    }

    for (UnsignedInt row = 0; row < 8; ++row)
    {
        solution[row] = augmented[row][8];
    }
    return true;
}

inline Real aphiMatrixMinAbsDiagonal(const Real matrix[8][8])
{
    Real min_abs = std::abs(matrix[0][0]);
    for (UnsignedInt row = 1; row < 8; ++row)
    {
        min_abs = std::min(min_abs, std::abs(matrix[row][row]));
    }
    return min_abs;
}

inline void aphiZeroMatrix3(Real matrix[3][3])
{
    for (UnsignedInt row = 0; row < 3; ++row)
    {
        for (UnsignedInt col = 0; col < 3; ++col)
        {
            matrix[row][col] = Real(0);
        }
    }
}

inline bool aphiSolve3x3(Real matrix[3][3], const Real rhs[3], Real solution[3])
{
    Real augmented[3][4];
    for (UnsignedInt row = 0; row < 3; ++row)
    {
        for (UnsignedInt col = 0; col < 3; ++col)
        {
            augmented[row][col] = matrix[row][col];
        }
        augmented[row][3] = rhs[row];
    }

    for (UnsignedInt pivot_col = 0; pivot_col < 3; ++pivot_col)
    {
        UnsignedInt pivot_row = pivot_col;
        Real pivot_abs_max = std::abs(augmented[pivot_col][pivot_col]);
        for (UnsignedInt row = pivot_col + 1; row < 3; ++row)
        {
            const Real candidate = std::abs(augmented[row][pivot_col]);
            if (candidate > pivot_abs_max)
            {
                pivot_abs_max = candidate;
                pivot_row = row;
            }
        }

        if (pivot_abs_max <= TinyReal)
        {
            return false;
        }

        if (pivot_row != pivot_col)
        {
            for (UnsignedInt col = pivot_col; col < 4; ++col)
            {
                const Real temp = augmented[pivot_col][col];
                augmented[pivot_col][col] = augmented[pivot_row][col];
                augmented[pivot_row][col] = temp;
            }
        }

        const Real inv_pivot = Real(1) / augmented[pivot_col][pivot_col];
        for (UnsignedInt col = pivot_col; col < 4; ++col)
        {
            augmented[pivot_col][col] *= inv_pivot;
        }

        for (UnsignedInt row = 0; row < 3; ++row)
        {
            if (row == pivot_col)
            {
                continue;
            }
            const Real factor = augmented[row][pivot_col];
            if (std::abs(factor) <= TinyReal)
            {
                continue;
            }
            for (UnsignedInt col = pivot_col; col < 4; ++col)
            {
                augmented[row][col] -= factor * augmented[pivot_col][col];
            }
        }
    }

    for (UnsignedInt row = 0; row < 3; ++row)
    {
        solution[row] = augmented[row][3];
    }
    return true;
}

/** Local uncorrected pairwise block for -grad(div A): B[c,alpha] = sum_j g_ij,c (g_ji,alpha + S_alpha). */
inline void aphiAccumulateGradDivABlockNeighbor(Matd &block, const Vecd &g_ij, const Vecd &g_ji, const Vecd &grad_sum)
{
    for (UnsignedInt alpha = 0; alpha < 3; ++alpha)
    {
        const Real grad_div_coeff = g_ji[alpha] + grad_sum[alpha];
        for (UnsignedInt component = 0; component < 3; ++component)
        {
            block(component, alpha) += g_ij[component] * grad_div_coeff;
        }
    }
}

inline void aphiAssembleCoupledPointBlock8x8(Real matrix[8][8], Real d_a, Real d_phi, Real reaction_coeff,
                                             Real phi_penalty, const Vecd &grad_phi_coupling,
                                             const Vecd &div_a_coupling, const Matd &graddiv_a_block,
                                             const AphiLhsTermFlags &terms, bool use_phi_gauge_penalty,
                                             bool use_a_divergence_penalty, Real a_divergence_penalty)
{
    aphiZeroMatrix8(matrix);

    if (terms.laplace_a || terms.reaction || use_a_divergence_penalty)
    {
        for (UnsignedInt row_component = 0; row_component < 3; ++row_component)
        {
            for (UnsignedInt col_component = 0; col_component < 3; ++col_component)
            {
                Real block_entry = Real(0);
                if (terms.laplace_a && row_component == col_component)
                {
                    block_entry += d_a;
                }
                if (use_a_divergence_penalty && terms.laplace_a)
                {
                    block_entry += a_divergence_penalty * graddiv_a_block(row_component, col_component);
                }
                matrix[row_component][col_component] += block_entry;
                matrix[row_component + 3][col_component + 3] += block_entry;
            }
            if (terms.reaction)
            {
                matrix[row_component][row_component + 3] -= reaction_coeff;
                matrix[row_component + 3][row_component] += reaction_coeff;
            }
        }
    }

    if (terms.grad_phi_coupling)
    {
        for (UnsignedInt component = 0; component < 3; ++component)
        {
            matrix[component][6] += grad_phi_coupling[component];
            matrix[component + 3][7] += grad_phi_coupling[component];
        }
    }

    if (terms.div_sigma_a_coupling)
    {
        for (UnsignedInt component = 0; component < 3; ++component)
        {
            matrix[6][component + 3] += div_a_coupling[component];
            matrix[7][component] -= div_a_coupling[component];
        }
    }

    if (terms.laplace_phi || use_phi_gauge_penalty)
    {
        const Real d_phi_total = d_phi + (use_phi_gauge_penalty ? phi_penalty : Real(0));
        matrix[6][6] += d_phi_total;
        matrix[7][7] += d_phi_total;
    }
}

} // namespace detail

template <typename... Parameters>
inline AphiComputeBlockJacobiDiagonalCK<Inner<Parameters...>>::AphiComputeBlockJacobiDiagonalCK(
    Inner<Parameters...> &inner_relation, const AphiMaterialNames &material_names, Real omega,
    const AphiLhsAssemblyOptions &options, const AphiBlockJacobiDiagonalNames &diag_names,
    Real pair_weight_regularization)
    : BaseInteraction(inner_relation), omega_(omega), options_(options), diag_names_(diag_names),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength())
{
    options_.omega = omega;
    auto &particles = this->particles_;
    particles->template registerStateVariable<Real>(diag_names_.laplace_a_diag, Real(0));
    particles->template registerStateVariable<Real>(diag_names_.laplace_phi_diag, Real(0));
    particles->template registerStateVariable<Vecd>(diag_names_.grad_phi_coupling, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(diag_names_.div_a_coupling, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Matd>(diag_names_.graddiv_a_block, Matd::Zero().eval());
    dv_sigma_ = particles->template getVariableByName<Real>(material_names.sigma);
    dv_nu_ = particles->template getVariableByName<Real>(material_names.nu);
    dv_laplace_a_diag_ = particles->template getVariableByName<Real>(diag_names_.laplace_a_diag);
    dv_laplace_phi_diag_ = particles->template getVariableByName<Real>(diag_names_.laplace_phi_diag);
    dv_grad_phi_coupling_ = particles->template getVariableByName<Vecd>(diag_names_.grad_phi_coupling);
    dv_div_a_coupling_ = particles->template getVariableByName<Vecd>(diag_names_.div_a_coupling);
    dv_graddiv_a_block_ = particles->template getVariableByName<Matd>(diag_names_.graddiv_a_block);
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiComputeBlockJacobiDiagonalCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      laplace_a_diag_(encloser.dv_laplace_a_diag_->DelegatedData(ex_policy)),
      laplace_phi_diag_(encloser.dv_laplace_phi_diag_->DelegatedData(ex_policy)),
      grad_phi_coupling_(encloser.dv_grad_phi_coupling_->DelegatedData(ex_policy)),
      div_a_coupling_(encloser.dv_div_a_coupling_->DelegatedData(ex_policy)),
      graddiv_a_block_(encloser.dv_graddiv_a_block_->DelegatedData(ex_policy)),
      use_a_divergence_penalty_(encloser.options_.use_a_divergence_penalty), omega_(encloser.omega_),
      terms_(encloser.options_.terms), use_phi_gauge_penalty_(encloser.options_.use_phi_gauge_penalty),
      phi_gauge_penalty_(encloser.options_.phi_gauge_penalty),
      pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename... Parameters>
inline void AphiComputeBlockJacobiDiagonalCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    Real laplace_a = Real(0);
    Real laplace_phi = Real(0);
    Vecd grad_phi_coupling(Real(0), Real(0), Real(0));
    Vecd div_a_coupling(Real(0), Real(0), Real(0));
    Vecd grad_sum(Real(0), Real(0), Real(0));
    Matd graddiv_block = Matd::Zero();
    const Real sigma_i = sigma_[index_i];
    const Real nu_i = nu_[index_i];
    const Real vol_i = Vol_[index_i];

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ij = this->dW_ij(index_i, index_j);
        const Real dW_ijV_j = dW_ij * Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, e_ij);

        if (terms_.laplace_phi)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, sigma_[index_j]);
            laplace_phi += sigma_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq,
                                                                          pair_weight_regularization_,
                                                                          reference_smoothing_length_);
        }

        if (terms_.laplace_a)
        {
            const Real nu_ij = AphiHarmonicMean(nu_i, nu_[index_j]);
            laplace_a += nu_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq,
                                                                     pair_weight_regularization_,
                                                                     reference_smoothing_length_);
        }

        if (terms_.grad_phi_coupling)
        {
            grad_phi_coupling += sigma_i * g_ij;
        }

        if (terms_.div_sigma_a_coupling)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, sigma_[index_j]);
            div_a_coupling += omega_ * sigma_ij * g_ij;
        }

        if (use_a_divergence_penalty_ && terms_.laplace_a)
        {
            grad_sum += g_ij;
        }
    }

    if (use_a_divergence_penalty_ && terms_.laplace_a)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Real dW_ij = this->dW_ij(index_i, index_j);
            const Real dW_ijV_j = dW_ij * Vol_[index_j];
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            const Vecd g_ji =
                -AphiPairwiseGradientWeightUncorrected(dW_ij * vol_i, this->e_ij(index_i, index_j));
            detail::aphiAccumulateGradDivABlockNeighbor(graddiv_block, g_ij, g_ji, grad_sum);
        }
    }

    laplace_a_diag_[index_i] = laplace_a;
    laplace_phi_diag_[index_i] = laplace_phi;
    grad_phi_coupling_[index_i] = grad_phi_coupling;
    div_a_coupling_[index_i] = div_a_coupling;
    graddiv_a_block_[index_i] = graddiv_block;
    (void)use_phi_gauge_penalty_;
    (void)phi_gauge_penalty_;
}

template <typename... Parameters>
inline AphiComputeBlockJacobiDiagonalCK<Contact<Parameters...>>::AphiComputeBlockJacobiDiagonalCK(
    Contact<Parameters...> &contact_relation, const AphiMaterialNames &material_names, Real omega,
    const AphiLhsAssemblyOptions &options, const AphiBlockJacobiDiagonalNames &diag_names,
    Real pair_weight_regularization)
    : BaseInteraction(contact_relation), omega_(omega), options_(options), diag_names_(diag_names),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength())
{
    options_.omega = omega;
    auto &particles = this->particles_;
    dv_sigma_ = particles->template getVariableByName<Real>(material_names.sigma);
    dv_nu_ = particles->template getVariableByName<Real>(material_names.nu);
    dv_laplace_a_diag_ = particles->template getVariableByName<Real>(diag_names_.laplace_a_diag);
    dv_laplace_phi_diag_ = particles->template getVariableByName<Real>(diag_names_.laplace_phi_diag);
    dv_grad_phi_coupling_ = particles->template getVariableByName<Vecd>(diag_names_.grad_phi_coupling);
    dv_div_a_coupling_ = particles->template getVariableByName<Vecd>(diag_names_.div_a_coupling);
    dv_graddiv_a_block_ = particles->template getVariableByName<Matd>(diag_names_.graddiv_a_block);
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_sigma_.push_back(contact_particles->template getVariableByName<Real>(material_names.sigma));
        dv_contact_nu_.push_back(contact_particles->template getVariableByName<Real>(material_names.nu));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiComputeBlockJacobiDiagonalCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      contact_sigma_(encloser.dv_contact_sigma_[contact_index]->DelegatedData(ex_policy)),
      contact_nu_(encloser.dv_contact_nu_[contact_index]->DelegatedData(ex_policy)),
      laplace_a_diag_(encloser.dv_laplace_a_diag_->DelegatedData(ex_policy)),
      laplace_phi_diag_(encloser.dv_laplace_phi_diag_->DelegatedData(ex_policy)),
      grad_phi_coupling_(encloser.dv_grad_phi_coupling_->DelegatedData(ex_policy)),
      div_a_coupling_(encloser.dv_div_a_coupling_->DelegatedData(ex_policy)),
      graddiv_a_block_(encloser.dv_graddiv_a_block_->DelegatedData(ex_policy)),
      use_a_divergence_penalty_(encloser.options_.use_a_divergence_penalty),
      contact_a_divergence_penalty_stencil_(encloser.options_.contact_a_divergence_penalty_stencil),
      omega_(encloser.omega_), terms_(encloser.options_.terms),
      pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename... Parameters>
inline void AphiComputeBlockJacobiDiagonalCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    Real laplace_a = Real(0);
    Real laplace_phi = Real(0);
    Vecd grad_phi_coupling(Real(0), Real(0), Real(0));
    Vecd div_a_coupling(Real(0), Real(0), Real(0));
    Vecd grad_sum(Real(0), Real(0), Real(0));
    Matd graddiv_block = Matd::Zero();
    const bool accumulate_contact_graddiv =
        use_a_divergence_penalty_ && terms_.laplace_a &&
        contact_a_divergence_penalty_stencil_ != AphiContactADivergencePenaltyStencilMode::InnerOnly;
    const Real sigma_i = sigma_[index_i];
    const Real nu_i = nu_[index_i];
    const Real vol_i = Vol_[index_i];

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ij = this->dW_ij(index_i, index_j);
        const Real dW_ijV_j = dW_ij * contact_Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, e_ij);

        if (terms_.laplace_phi)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
            laplace_phi += sigma_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq,
                                                                          pair_weight_regularization_,
                                                                          reference_smoothing_length_);
        }

        if (terms_.laplace_a)
        {
            const Real nu_ij = AphiHarmonicMean(nu_i, contact_nu_[index_j]);
            laplace_a += nu_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq,
                                                                     pair_weight_regularization_,
                                                                     reference_smoothing_length_);
        }

        if (terms_.grad_phi_coupling)
        {
            grad_phi_coupling += sigma_i * g_ij;
        }

        if (terms_.div_sigma_a_coupling)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
            div_a_coupling += omega_ * sigma_ij * g_ij;
        }

        if (accumulate_contact_graddiv)
        {
            grad_sum += g_ij;
        }
    }

    if (accumulate_contact_graddiv)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Real dW_ij = this->dW_ij(index_i, index_j);
            const Real dW_ijV_j = dW_ij * contact_Vol_[index_j];
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            const Vecd g_ji =
                -AphiPairwiseGradientWeightUncorrected(dW_ij * vol_i, this->e_ij(index_i, index_j));
            detail::aphiAccumulateGradDivABlockNeighbor(graddiv_block, g_ij, g_ji, grad_sum);
        }
    }

    laplace_a_diag_[index_i] += laplace_a;
    laplace_phi_diag_[index_i] += laplace_phi;
    grad_phi_coupling_[index_i] += grad_phi_coupling;
    div_a_coupling_[index_i] += div_a_coupling;
    graddiv_a_block_[index_i] += graddiv_block;
}

inline RegisterAphiBlockJacobiDiagnosticFieldsCK::RegisterAphiBlockJacobiDiagnosticFieldsCK(
    SPHBody &sph_body, const AphiBlockJacobiDiagnosticNames &diag_names)
    : LocalDynamics(sph_body)
{
    auto &particles = particles_;
    particles->template registerStateVariable<Real>(diag_names.fallback_flag, Real(0));
    particles->template registerStateVariable<Real>(diag_names.min_pivot, Real(0));
}

inline AphiApplyBlockJacobiInverseCK::AphiApplyBlockJacobiInverseCK(
    SPHBody &sph_body, const AphiBlockNames &src_names, const AphiBlockNames &dst_names,
    const AphiMaterialNames &material_names, Real omega, const AphiLhsAssemblyOptions &options,
    AphiBlockJacobiPreconditionerKind preconditioner_kind, const AphiBlockJacobiDiagonalNames &diag_names,
    const AphiBlockJacobiDiagnosticNames *diagnostic_names)
    : LocalDynamics(sph_body), omega_(omega), options_(options), preconditioner_kind_(preconditioner_kind),
      record_diagnostics_(diagnostic_names != nullptr)
{
    auto &particles = particles_;
    dv_sigma_ = particles->template getVariableByName<Real>(material_names.sigma);
    dv_laplace_a_diag_ = particles->template getVariableByName<Real>(diag_names.laplace_a_diag);
    dv_laplace_phi_diag_ = particles->template getVariableByName<Real>(diag_names.laplace_phi_diag);
    dv_grad_phi_coupling_ = particles->template getVariableByName<Vecd>(diag_names.grad_phi_coupling);
    dv_div_a_coupling_ = particles->template getVariableByName<Vecd>(diag_names.div_a_coupling);
    dv_graddiv_a_block_ = particles->template getVariableByName<Matd>(diag_names.graddiv_a_block);
    dv_src_a_real_ = particles->template getVariableByName<Vecd>(src_names.a_real);
    dv_src_a_imag_ = particles->template getVariableByName<Vecd>(src_names.a_imag);
    dv_src_phi_real_ = particles->template getVariableByName<Real>(src_names.phi_real);
    dv_src_phi_imag_ = particles->template getVariableByName<Real>(src_names.phi_imag);
    dv_dst_a_real_ = particles->template getVariableByName<Vecd>(dst_names.a_real);
    dv_dst_a_imag_ = particles->template getVariableByName<Vecd>(dst_names.a_imag);
    dv_dst_phi_real_ = particles->template getVariableByName<Real>(dst_names.phi_real);
    dv_dst_phi_imag_ = particles->template getVariableByName<Real>(dst_names.phi_imag);
    if (record_diagnostics_)
    {
        dv_fallback_flag_ = particles->template getVariableByName<Real>(diagnostic_names->fallback_flag);
        dv_min_pivot_ = particles->template getVariableByName<Real>(diagnostic_names->min_pivot);
    }
    else
    {
        dv_fallback_flag_ = nullptr;
        dv_min_pivot_ = nullptr;
    }
}

template <class ExecutionPolicy, class EncloserType>
inline AphiApplyBlockJacobiInverseCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                 EncloserType &encloser)
    : omega_(encloser.omega_), terms_(encloser.options_.terms),
      use_phi_gauge_penalty_(encloser.options_.use_phi_gauge_penalty),
      phi_gauge_penalty_(encloser.options_.phi_gauge_penalty),
      use_a_divergence_penalty_(encloser.options_.use_a_divergence_penalty),
      a_divergence_penalty_(encloser.options_.a_divergence_penalty),
      preconditioner_kind_(encloser.preconditioner_kind_), sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      laplace_a_diag_(encloser.dv_laplace_a_diag_->DelegatedData(ex_policy)),
      laplace_phi_diag_(encloser.dv_laplace_phi_diag_->DelegatedData(ex_policy)),
      grad_phi_coupling_(encloser.dv_grad_phi_coupling_->DelegatedData(ex_policy)),
      div_a_coupling_(encloser.dv_div_a_coupling_->DelegatedData(ex_policy)),
      graddiv_a_block_(encloser.dv_graddiv_a_block_->DelegatedData(ex_policy)),
      src_a_real_(encloser.dv_src_a_real_->DelegatedData(ex_policy)),
      src_a_imag_(encloser.dv_src_a_imag_->DelegatedData(ex_policy)),
      src_phi_real_(encloser.dv_src_phi_real_->DelegatedData(ex_policy)),
      src_phi_imag_(encloser.dv_src_phi_imag_->DelegatedData(ex_policy)),
      dst_a_real_(encloser.dv_dst_a_real_->DelegatedData(ex_policy)),
      dst_a_imag_(encloser.dv_dst_a_imag_->DelegatedData(ex_policy)),
      dst_phi_real_(encloser.dv_dst_phi_real_->DelegatedData(ex_policy)),
      dst_phi_imag_(encloser.dv_dst_phi_imag_->DelegatedData(ex_policy)),
      fallback_flag_(encloser.record_diagnostics_ ? encloser.dv_fallback_flag_->DelegatedData(ex_policy) : nullptr),
      min_pivot_(encloser.record_diagnostics_ ? encloser.dv_min_pivot_->DelegatedData(ex_policy) : nullptr)
{
}

inline void AphiApplyBlockJacobiInverseCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    if (preconditioner_kind_ == AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8)
    {
        Real matrix[8][8];
        Real rhs[8];
        Real solution[8];
        const Vecd a_re_in = src_a_real_[index_i];
        const Vecd a_im_in = src_a_imag_[index_i];
        rhs[0] = a_re_in[0];
        rhs[1] = a_re_in[1];
        rhs[2] = a_re_in[2];
        rhs[3] = a_im_in[0];
        rhs[4] = a_im_in[1];
        rhs[5] = a_im_in[2];
        rhs[6] = src_phi_real_[index_i];
        rhs[7] = src_phi_imag_[index_i];

        detail::aphiAssembleCoupledPointBlock8x8(
            matrix, laplace_a_diag_[index_i], laplace_phi_diag_[index_i], omega_ * sigma_[index_i], phi_gauge_penalty_,
            grad_phi_coupling_[index_i], div_a_coupling_[index_i], graddiv_a_block_[index_i], terms_,
            use_phi_gauge_penalty_, use_a_divergence_penalty_, a_divergence_penalty_);

        const Real min_abs_diag = detail::aphiMatrixMinAbsDiagonal(matrix);
        if (min_pivot_ != nullptr)
        {
            min_pivot_[index_i] = min_abs_diag;
        }

        if (detail::aphiSolve8x8(matrix, rhs, solution))
        {
            if (fallback_flag_ != nullptr)
            {
                fallback_flag_[index_i] = Real(0);
            }
            dst_a_real_[index_i] = Vecd(solution[0], solution[1], solution[2]);
            dst_a_imag_[index_i] = Vecd(solution[3], solution[4], solution[5]);
            dst_phi_real_[index_i] = solution[6];
            dst_phi_imag_[index_i] = solution[7];
            return;
        }

        if (fallback_flag_ != nullptr)
        {
            fallback_flag_[index_i] = Real(1);
        }
    }
    else if (fallback_flag_ != nullptr)
    {
        fallback_flag_[index_i] = Real(0);
        if (min_pivot_ != nullptr)
        {
            min_pivot_[index_i] = Real(0);
        }
    }

    if (terms_.laplace_a || terms_.reaction)
    {
        const Real d = laplace_a_diag_[index_i];
        const Vecd a_re_in = src_a_real_[index_i];
        const Vecd a_im_in = src_a_imag_[index_i];
        if (use_a_divergence_penalty_ && terms_.laplace_a && !terms_.reaction)
        {
            Real matrix_re[3][3];
            Real matrix_im[3][3];
            Real rhs_re[3];
            Real rhs_im[3];
            Real sol_re[3];
            Real sol_im[3];
            detail::aphiZeroMatrix3(matrix_re);
            detail::aphiZeroMatrix3(matrix_im);
            for (UnsignedInt row = 0; row < 3; ++row)
            {
                for (UnsignedInt col = 0; col < 3; ++col)
                {
                    const Real penalty_entry =
                        a_divergence_penalty_ * graddiv_a_block_[index_i](row, col);
                    const Real laplace_entry = (row == col ? d : Real(0));
                    matrix_re[row][col] = laplace_entry + penalty_entry;
                    matrix_im[row][col] = laplace_entry + penalty_entry;
                }
            }
            rhs_re[0] = a_re_in[0];
            rhs_re[1] = a_re_in[1];
            rhs_re[2] = a_re_in[2];
            rhs_im[0] = a_im_in[0];
            rhs_im[1] = a_im_in[1];
            rhs_im[2] = a_im_in[2];
            if (detail::aphiSolve3x3(matrix_re, rhs_re, sol_re) && detail::aphiSolve3x3(matrix_im, rhs_im, sol_im))
            {
                dst_a_real_[index_i] = Vecd(sol_re[0], sol_re[1], sol_re[2]);
                dst_a_imag_[index_i] = Vecd(sol_im[0], sol_im[1], sol_im[2]);
            }
            else
            {
                dst_a_real_[index_i] = a_re_in;
                dst_a_imag_[index_i] = a_im_in;
            }
        }
        else if (terms_.reaction)
        {
            const Real c = omega_ * sigma_[index_i];
            const Real det = d * d + c * c + TinyReal;
            const Real inv_det = Real(1) / det;
            dst_a_real_[index_i] = inv_det * (d * a_re_in + c * a_im_in);
            dst_a_imag_[index_i] = inv_det * (-c * a_re_in + d * a_im_in);
        }
        else
        {
            const Real inv_d = Real(1) / (d + TinyReal);
            dst_a_real_[index_i] = inv_d * a_re_in;
            dst_a_imag_[index_i] = inv_d * a_im_in;
        }
    }
    else
    {
        dst_a_real_[index_i] = src_a_real_[index_i];
        dst_a_imag_[index_i] = src_a_imag_[index_i];
    }

    if (terms_.laplace_phi || use_phi_gauge_penalty_)
    {
        const Real d_phi = laplace_phi_diag_[index_i] +
                           (use_phi_gauge_penalty_ ? phi_gauge_penalty_ : Real(0));
        const Real inv_d_phi = Real(1) / (d_phi + TinyReal);
        dst_phi_real_[index_i] = inv_d_phi * src_phi_real_[index_i];
        dst_phi_imag_[index_i] = inv_d_phi * src_phi_imag_[index_i];
    }
    else
    {
        dst_phi_real_[index_i] = src_phi_real_[index_i];
        dst_phi_imag_[index_i] = src_phi_imag_[index_i];
    }
}

template <class ExecutionPolicy>
inline AphiComputeBlockJacobiContactDynamicsBundle<ExecutionPolicy>::AphiComputeBlockJacobiContactDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiMaterialNames &material_names,
    Real omega, const AphiLhsAssemblyOptions &options, const AphiBlockJacobiDiagonalNames &diag_names,
    Real pair_weight_regularization)
    : compute_inner_(DynamicsArgs(inner_relation, material_names, omega, options)),
      compute_contact_(DynamicsArgs(contact_relation, material_names, omega, options))
{
    (void)sph_body;
    (void)diag_names;
    (void)pair_weight_regularization;
}

template <class ExecutionPolicy>
inline void AphiComputeBlockJacobiContactDynamicsBundle<ExecutionPolicy>::exec()
{
    compute_inner_.exec();
    compute_contact_.exec();
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_JACOBI_PRECONDITIONER_CK_HPP
