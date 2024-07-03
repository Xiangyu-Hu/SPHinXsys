#include "kernel_correction.hpp"

namespace SPH
{
//=================================================================================================//
void LinearGradientCorrectionMatrix<Inner<>>::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity();

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ij.transpose();
    }
    B_[index_i] = local_configuration;
}
//=================================================================================================//
void LinearGradientCorrectionMatrix<Inner<>>::update(size_t index_i, Real dt)
{
    Real det_sqr = SMAX(alpha_ - B_[index_i].determinant(), Real(0));
    Matd inverse = B_[index_i].inverse();
    Real weight1_ = B_[index_i].determinant() / (B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (B_[index_i].determinant() + det_sqr);
    B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
}
//=================================================================================================//
LinearGradientCorrectionMatrix<Contact<>>::
    LinearGradientCorrectionMatrix(BaseContactRelation &contact_relation)
    : LinearGradientCorrectionMatrix<DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(contact_particles_[k]->getVariableByName<Real>("Mass"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void LinearGradientCorrectionMatrix<Contact<>>::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
    }
    B_[index_i] += local_configuration;
}
//=================================================================================================//
KernelGradientCorrection<Inner<>>::
    KernelGradientCorrection(BaseInnerRelation &inner_relation)
    : KernelGradientCorrection<DataDelegateInner>(inner_relation),
      average_correction_matrix_(*particles_->getVariableByName<Matd>("LinearGradientCorrectionMatrix")){};
//=================================================================================================//
void KernelGradientCorrection<Inner<>>::interaction(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    correctKernelGradient(average_correction_matrix_, inner_neighborhood, index_i);
}
//=================================================================================================//
KernelGradientCorrection<Contact<>>::
    KernelGradientCorrection(BaseContactRelation &contact_relation)
    : KernelGradientCorrection<DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_average_correction_matrix_.push_back(
            PairAverageVariable<Matd>(
                *particles_->getVariableByName<Matd>("LinearGradientCorrectionMatrix"),
                *contact_particles_[k]->getVariableByName<Matd>("LinearGradientCorrectionMatrix")));
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
