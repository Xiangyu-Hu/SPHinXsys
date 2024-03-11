#include "general_interaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
KernelCorrectionMatrixInner::
    KernelCorrectionMatrixInner(BaseInnerRelation &inner_relation, Real alpha)
    : LocalDynamics(inner_relation.getSPHBody()),
      GeneralDataDelegateInner(inner_relation),
      alpha_(alpha), B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")),
      Vol_(this->particles_->Vol_) {}
//=================================================================================================//
void KernelCorrectionMatrixInner::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity();

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ij.transpose();
    }
    B_[index_i] = local_configuration;
}
//=================================================================================================//
void KernelCorrectionMatrixInner::update(size_t index_i, Real dt)
{
    ///* WKGC1 */
    Real det_sqr = SMAX(alpha_ - B_[index_i].determinant(), 0.0);
    Matd inverse = B_[index_i].inverse();
    Real weight1_ = B_[index_i].determinant() / (B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (B_[index_i].determinant() + det_sqr);
    B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
    
    /* WKGC2 */
    /* Real det_sqr = alpha_;
    Matd inverse = B_[index_i].inverse();
    Real weight1_ = B_[index_i].determinant() / (B_[index_i].determinant() + det_sqr);
    Real weight2_ = det_sqr / (B_[index_i].determinant() + det_sqr);
    B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();*/
}
//=================================================================================================//
KernelCorrectionMatrixInnerWithLevelSet::KernelCorrectionMatrixInnerWithLevelSet(BaseInnerRelation& inner_relation)
    : KernelCorrectionMatrixInner(inner_relation), pos_(this->particles_->pos_),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void KernelCorrectionMatrixInnerWithLevelSet::interaction(size_t index_i, Real dt)
{
    KernelCorrectionMatrixInner::interaction(index_i, dt);
    Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
    B_[index_i] -= level_set_shape_->computeDisplacementKernelGradientIntegral(pos_[index_i],
        sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
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
        StdLargeVec<Real>& vol_k = *(this->contact_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd gradW_ij = contact_neighborhood.dW_ijV_j_[n] * vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
    }
    B_[index_i] += local_configuration;
}
//=================================================================================================//    
KernelCorrectionMatrixComplexWithLevelSet::
KernelCorrectionMatrixComplexWithLevelSet(ComplexRelation& complex_relation, const std::string& shape_name)
    : KernelCorrectionMatrixComplex(complex_relation), pos_(this->particles_->pos_),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    ComplexShape& complex_shape = DynamicCast<ComplexShape>(this, *sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
}
//=================================================================================================//   
void KernelCorrectionMatrixComplexWithLevelSet::interaction(size_t index_i, Real dt)
{
    KernelCorrectionMatrixComplex::interaction(index_i, dt);
    Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
    B_[index_i] -= level_set_shape_->computeDisplacementKernelGradientIntegral(pos_[index_i],
        sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
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
