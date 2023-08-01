#include "transport_velocity_correction.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
TransportVelocityCorrectionInnerAdaptive::
    TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation &inner_relation, Real coefficient)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      sph_adaptation_(*sph_body_.sph_adaptation_),
      pos_(particles_->pos_), indicator_(*particles_->getVariableByName<int>("Indicator")),
      smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
      coefficient_(coefficient) {}

//=================================================================================================//
void TransportVelocityCorrectionInnerAdaptive::
    interaction(size_t index_i, Real dt)
{
    Vecd acceleration_trans = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        // acceleration for transport velocity
        acceleration_trans -= 2.0 * nablaW_ijV_j;
    }

    if (indicator_[index_i] == 0)
    {
        Real inv_h_ratio = 1.0 / sph_adaptation_.SmoothingLengthRatio(index_i);
        pos_[index_i] += coefficient_ * smoothing_length_sqr_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
void TransportVelocityCorrectionComplexAdaptive::interaction(size_t index_i, Real dt)
{
    TransportVelocityCorrectionInnerAdaptive::interaction(index_i, dt);

    Vecd acceleration_trans = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

            // acceleration for transport velocity
            acceleration_trans -= 2.0 * nablaW_ijV_j;
        }
    }

    /** correcting particle position */
    if (indicator_[index_i] == 0)
    {
        Real inv_h_ratio = 1.0 / sph_adaptation_.SmoothingLengthRatio(index_i);
        pos_[index_i] += coefficient_ * smoothing_length_sqr_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
