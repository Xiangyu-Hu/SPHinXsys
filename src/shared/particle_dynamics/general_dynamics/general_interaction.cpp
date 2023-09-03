#include "general_interaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
KernelCorrectionMatrixInner::
    KernelCorrectionMatrixInner(BaseInnerRelation &inner_relation, int beta, Real alpha)
    : LocalDynamics(inner_relation.getSPHBody()),
      GeneralDataDelegateInner(inner_relation),
      beta_(beta), alpha_(alpha),
      B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")) {}
//=================================================================================================//
void KernelCorrectionMatrixInner::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity();

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ij.transpose();
    }

    B_[index_i] = local_configuration;
}
//=================================================================================================//
void KernelCorrectionMatrixInner::update(size_t index_i, Real dt)
{
    Real det_sqr = pow(B_[index_i].determinant(), beta_);
    Matd inverse = B_[index_i].inverse();
    B_[index_i] = (det_sqr * inverse + alpha_ * Matd::Identity()) / (alpha_ + det_sqr);
}
//=================================================================================================//
KernelCorrectionMatrixComplex::
    KernelCorrectionMatrixComplex(ComplexRelation &complex_relation, int beta, Real alpha)
    : KernelCorrectionMatrixInner(complex_relation.getInnerRelation(), beta, alpha),
      GeneralDataDelegateContactOnly(complex_relation.getContactRelation())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
void KernelCorrectionMatrixComplex::interaction(size_t index_i, Real dt)
{
    KernelCorrectionMatrixInner::interaction(index_i, dt);

    Matd local_configuration = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Vecd gradW_ij = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
    }
    B_[index_i] += local_configuration;
}
//=================================================================================================//
KernelGradientCorrectionInner::
    KernelGradientCorrectionInner(KernelCorrectionMatrixInner &kernel_correction_inner)
    : LocalDynamics(kernel_correction_inner.getSPHBody()),
      GeneralDataDelegateInner(kernel_correction_inner.getInnerRelation()),
      B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")){};
//=================================================================================================//
void KernelGradientCorrectionInner::interaction(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    correctKernelGradient(inner_neighborhood, index_i);
}
//=================================================================================================//
void KernelGradientCorrectionInner::correctKernelGradient(Neighborhood &neighborhood, size_t index_i)
{
    for (size_t n = 0; n != neighborhood.current_size_; ++n)
    {
        size_t index_j = neighborhood.j_[n];

        Vecd corrected_direction = 0.5 * (B_[index_i] + B_[index_j]) * neighborhood.e_ij_[n];
        Real direction_norm = corrected_direction.norm();
        neighborhood.dW_ijV_j_[n] *= direction_norm;
        neighborhood.e_ij_[n] = corrected_direction / (direction_norm + Eps);
        neighborhood.r_ij_[n] *= direction_norm;
    }
}
//=================================================================================================//
KernelGradientCorrectionComplex::
    KernelGradientCorrectionComplex(KernelCorrectionMatrixComplex &kernel_correction_complex)
    : KernelGradientCorrectionInner(kernel_correction_complex),
      GeneralDataDelegateContactOnly(kernel_correction_complex.getContactRelation()) {}
//=================================================================================================//
void KernelGradientCorrectionComplex::interaction(size_t index_i, Real dt)
{
    KernelGradientCorrectionInner::interaction(index_i, dt);

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        KernelGradientCorrectionInner::correctKernelGradient(contact_neighborhood, index_i);
    }
}
//=================================================================================================//
} // namespace SPH
