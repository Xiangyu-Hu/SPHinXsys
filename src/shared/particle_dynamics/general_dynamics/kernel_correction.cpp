#include "kernel_correction.hpp"

namespace SPH
{
//=================================================================================================//
void KernelCorrectionMatrix<Inner<>>::interaction(size_t index_i, Real dt)
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
void KernelCorrectionMatrix<Inner<>>::update(size_t index_i, Real dt)
{
    Matd inverse = B_[index_i].inverse();
    Real weight = alpha_ / (B_[index_i].determinant() + alpha_);
    B_[index_i] = weight * Matd::Identity() + (1.0 - weight) * inverse;
}
//=================================================================================================//
KernelCorrectionMatrix<Contact<>>::
    KernelCorrectionMatrix(BaseContactRelation &contact_relation)
    : KernelCorrectionMatrix<GeneralDataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
void KernelCorrectionMatrix<Contact<>>::interaction(size_t index_i, Real dt)
{
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
KernelGradientCorrection<Inner<>>::
    KernelGradientCorrection(BaseInnerRelation &inner_relation)
    : KernelGradientCorrection<GeneralDataDelegateInner>(inner_relation),
      average_correction_matrix_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")){};
//=================================================================================================//
void KernelGradientCorrection<Inner<>>::interaction(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    correctKernelGradient(average_correction_matrix_, inner_neighborhood, index_i);
}
//=================================================================================================//
KernelGradientCorrection<Contact<>>::
    KernelGradientCorrection(BaseContactRelation &contact_relation)
    : KernelGradientCorrection<GeneralDataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_average_correction_matrix_.push_back(
            PairAverageVariable<Matd>(
                *particles_->getVariableByName<Matd>("KernelCorrectionMatrix"),
                *contact_particles_[k]->getVariableByName<Matd>("KernelCorrectionMatrix")));
    }
}
//=================================================================================================//
void KernelGradientCorrection<Contact<>>::interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        correctKernelGradient(contact_average_correction_matrix_[k], contact_neighborhood, index_i);
    }
}
//=================================================================================================//
} // namespace SPH
