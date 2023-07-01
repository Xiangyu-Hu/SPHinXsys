
#include "common_shared_eulerian_classes.h"

namespace SPH
{
//=================================================================================================//
KernelGradientWithCorrectionInner::KernelGradientWithCorrectionInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation)
{
    particles_->registerVariable(B_, "CorrectionMatrix");
    particles_->registerVariable(local_configuration_inner_, "LocalConfigurationInner");
};
//=================================================================================================//
void KernelGradientWithCorrectionInner::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ijV_j.transpose();
    }
    B_[index_i] = local_configuration.inverse();
    local_configuration_inner_[index_i] = local_configuration;
}
//=================================================================================================//
void KernelGradientWithCorrectionInner::update(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd r_ji = r_ij * e_ij;

        Matd B_average = 0.5 * (B_[index_i] + B_[index_j]);
        Vecd kernel_gradient_with_B = B_average * dW_ijV_j * e_ij;
        if (dW_ijV_j < 0)
        {
            inner_neighborhood.dW_ijV_j_[n] = -kernel_gradient_with_B.norm();
        }
        else
        {
            inner_neighborhood.dW_ijV_j_[n] = kernel_gradient_with_B.norm();
        }
        inner_neighborhood.e_ij_[n] = kernel_gradient_with_B / inner_neighborhood.dW_ijV_j_[n];
        inner_neighborhood.r_ij_[n] = r_ji.dot(inner_neighborhood.e_ij_[n]);
    }
}
//=================================================================================================//
void KernelGradientWithCorrectionComplex::interaction(size_t index_i, Real dt)
{
    KernelGradientWithCorrectionInner::interaction(index_i, dt);

    Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd r_ji = r_ij * e_ij;
            Vecd gradW_ijV_j = dW_ijV_j * e_ij;

            local_configuration -= r_ji * gradW_ijV_j.transpose();
        }
    }
    B_[index_i] = (local_configuration + local_configuration_inner_[index_i]).inverse();
}
//=================================================================================================//
void KernelGradientWithCorrectionComplex::update(size_t index_i, Real dt)
{
    KernelGradientWithCorrectionInner::update(index_i, dt);

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd r_ji = r_ij * e_ij;

            Vecd kernel_gradient_with_B = B_[index_i] * dW_ijV_j * e_ij;
            if (dW_ijV_j < 0)
            {
                contact_neighborhood.dW_ijV_j_[n] = -kernel_gradient_with_B.norm();
            }
            else
            {
                contact_neighborhood.dW_ijV_j_[n] = kernel_gradient_with_B.norm();
            }
            contact_neighborhood.e_ij_[n] = kernel_gradient_with_B / contact_neighborhood.dW_ijV_j_[n];
            contact_neighborhood.r_ij_[n] = r_ji.dot(contact_neighborhood.e_ij_[n]);
        }
    }
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//