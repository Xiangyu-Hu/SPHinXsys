#include "general_interaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
KernelCorrectionMatrixInner::
    KernelCorrectionMatrixInner(BaseInnerRelation &inner_relation, Real alpha)
    : LocalDynamics(inner_relation.getSPHBody()),
      GeneralDataDelegateInner(inner_relation),
      alpha_(alpha), B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")) {}
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
    Real det_sqr = alpha_;
    Matd inverse = B_[index_i].inverse();
    Real weight1_ = B_[index_i].determinant() / (B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (B_[index_i].determinant() + det_sqr);
    B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
}
//=================================================================================================//
KernelCorrectionMatrixComplex::
    KernelCorrectionMatrixComplex(ComplexRelation &complex_relation, Real alpha)
    : KernelCorrectionMatrixInner(complex_relation.getInnerRelation(), alpha),
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
      average_correction_matrix_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")){};
//=================================================================================================//
void KernelGradientCorrectionInner::interaction(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    correctKernelGradient(average_correction_matrix_, inner_neighborhood, index_i);
}
//=================================================================================================//
KernelGradientCorrectionComplex::
    KernelGradientCorrectionComplex(KernelCorrectionMatrixComplex &kernel_correction_complex)
    : KernelGradientCorrectionInner(kernel_correction_complex),
      GeneralDataDelegateContactOnly(kernel_correction_complex.getContactRelation())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_average_correction_matrix_.push_back(
            ParticlesPairAverageContact<Matd>(
                *particles_->getVariableByName<Matd>("KernelCorrectionMatrix"),
                *contact_particles_[k]->getVariableByName<Matd>("KernelCorrectionMatrix")));
    }
}
//=================================================================================================//
void KernelGradientCorrectionComplex::interaction(size_t index_i, Real dt)
{
    KernelGradientCorrectionInner::interaction(index_i, dt);

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        correctKernelGradient(contact_average_correction_matrix_[k], contact_neighborhood, index_i);
    }
}
//=================================================================================================//
} // namespace SPH
