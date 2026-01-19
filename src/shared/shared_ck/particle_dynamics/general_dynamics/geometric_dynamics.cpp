#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
NormalFromBodyShapeCK::NormalFromBodyShapeCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      initial_shape_(&sph_body.getInitialShape()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_n_(particles_->registerStateVariable<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariable<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(particles_->registerStateVariable<Real>("SignedDistance")),
      dv_phi0_(particles_->registerStateVariable<Real>("InitialSignedDistance", dv_phi_)) {}
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
NormalFromSubShapeAndOpCK::NormalFromSubShapeAndOpCK(
    SPHBody &sph_body, ComplexShape &complex_shape, const std::string &shape_name)
    : LocalDynamics(sph_body),
      shape_and_op_(complex_shape.getSubShapeAndOpByName(shape_name)),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_n_(particles_->registerStateVariable<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariable<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(particles_->registerStateVariable<Real>("SignedDistance")),
      dv_phi0_(particles_->registerStateVariable<Real>("InitialSignedDistance", dv_phi_)) {}
//=================================================================================================//
NormalFromSubShapeAndOpCK::NormalFromSubShapeAndOpCK(SPHBody &sph_body, const std::string &shape_name)
    : NormalFromSubShapeAndOpCK(
          sph_body, DynamicCast<ComplexShape>(this, sph_body.getInitialShape()), shape_name) {}
//=================================================================================================//
void NormalFromSubShapeAndOpCK::UpdateKernel::update(size_t index_i, Real /*dt*/)
{
    Vecd normal_direction = switch_sign_ * shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = switch_sign_ * shape_->findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=============================================================================================//
SurfaceIndicationFromBodyShape::SurfaceIndicationFromBodyShape(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      initial_shape_(&sph_body.getInitialShape()),
      spacing_ref_(sph_body.getSPHAdaptation().ReferenceSpacing()),
      dv_indicator_(particles_->registerStateVariable<int>("SurfaceIndicator")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")) {}
//=============================================================================================//
void SurfaceIndicationFromBodyShape::UpdateKernel::update(size_t index_i, Real dt)
{
    Real signed_distance = initial_shape_->findSignedDistance(pos_[index_i]);
    indicator_[index_i] = signed_distance > -spacing_ref_ ? 1 : 0;
}
//=============================================================================================//
RandomizeParticlePositionCK::RandomizeParticlePositionCK(SPHBody &sph_body, Real randomize_factor)
    : LocalDynamics(sph_body), dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      randomize_scale_(sph_body.getSPHAdaptation().MinimumSpacing() * randomize_factor) {}
//=============================================================================================//
void RandomizeParticlePositionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    for (int k = 0; k != Dimensions; ++k)
    {
        pos_[index_i][k] += rand_uniform(-1.0, 1.0) * randomize_scale_;
    }
}
//=================================================================================================//
} // namespace SPH
