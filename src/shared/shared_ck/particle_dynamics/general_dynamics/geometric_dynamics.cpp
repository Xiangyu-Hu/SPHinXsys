#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
NormalFromBodyShapeCK::NormalFromBodyShapeCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      initial_shape_(&sph_body.getInitialShape()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_n_(particles_->registerStateVariableOnly<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariableOnly<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(particles_->registerStateVariableOnly<Real>("SignedDistance")),
      dv_phi0_(particles_->registerStateVariableOnly<Real>("InitialSignedDistance", dv_phi_)) {}
//=============================================================================================//
void NormalFromBodyShapeCK::UpdateKernel::update(size_t index_i, Real dt)
{
    Vecd normal_direction = initial_shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = initial_shape_->findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=================================================================================================//
} // namespace SPH
