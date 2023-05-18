#include "general_dynamics_wkgc.h"

namespace SPH
{
//=================================================================================================//
void CorrectionMatrixInner::interaction(size_t index_i, Real dt)
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
void CorrectionMatrixInner::update(size_t index_i, Real dt)
{
    Real det_sqr = pow(B_[index_i].determinant(), beta_);
    Matd inverse = B_[index_i].inverse();
    B_[index_i] = (det_sqr * inverse + alpha_ * Matd::Identity()) / (alpha_ + det_sqr);
}
//=================================================================================================//
void CorrectionMatrixContact::interaction(size_t index_i, Real dt)
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
} // namespace SPH
  //=================================================================================================//