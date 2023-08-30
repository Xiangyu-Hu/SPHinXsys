#include "general_interaction.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
CorrectedConfigurationInner::
    CorrectedConfigurationInner(BaseInnerRelation &inner_relation, int beta, Real alpha)
    : LocalDynamics(inner_relation.getSPHBody()),
      GeneralDataDelegateInner(inner_relation),
      beta_(beta), alpha_(alpha),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      B_(*particles_->getVariableByName<Matd>("CorrectionMatrix")) {}
//=================================================================================================//
void CorrectedConfigurationInner::interaction(size_t index_i, Real dt)
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
void CorrectedConfigurationInner::update(size_t index_i, Real dt)
{
    Real det_sqr = pow(B_[index_i].determinant(), beta_);
    Matd inverse = B_[index_i].inverse();
    Real weight1_ = det_sqr / (alpha_ * pow(smoothing_length_, beta_) + det_sqr);
    Real weight2_ = alpha_ * pow(smoothing_length_, beta_) / (alpha_ * pow(smoothing_length_, beta_) + det_sqr);
    B_[index_i] = weight1_ * inverse + weight2_ * Matd::Identity();
}
//=================================================================================================//
CorrectedConfigurationComplex::
    CorrectedConfigurationComplex(ComplexRelation &complex_relation, int beta, Real alpha)
    : CorrectedConfigurationInner(complex_relation.getInnerRelation(), beta, alpha),
      GeneralDataDelegateContactOnly(complex_relation.getContactRelation())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
void CorrectedConfigurationComplex::interaction(size_t index_i, Real dt)
{
    CorrectedConfigurationInner::interaction(index_i, dt);

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
} // namespace SPH
  //=================================================================================================//