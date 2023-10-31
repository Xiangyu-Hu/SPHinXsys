#include "general_geometric.h"
#include "base_particles.hpp"

namespace SPH
{
//=============================================================================================//
NormalDirectionFromBodyShape::NormalDirectionFromBodyShape(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      body_shape_(*sph_body.body_shape_), pos_(particles_->pos_),
      n_(*particles_->registerSharedVariable<Vecd>("NormalDirection")),
      n0_(*particles_->registerSharedVariable<Vecd>("InitialNormalDirection")) {}
//=============================================================================================//
void NormalDirectionFromBodyShape::update(size_t index_i, Real dt)
{
    Vecd normal_direction = body_shape_.findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
}
//=============================================================================================//
NormalDirectionFromShapeAndOp::
    NormalDirectionFromShapeAndOp(SPHBody &sph_body, const std::string &shape_name)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      shape_and_op_(DynamicCast<ComplexShape>(this, sph_body.body_shape_)->getShapeAndOpByName(shape_name)),
      shape_(shape_and_op_->first),
      switch_sign_(shape_and_op_->second == ShapeBooleanOps::add ? 1.0 : -1.0),
      pos_(particles_->pos_),
      n_(*particles_->registerSharedVariable<Vecd>("NormalDirection")),
      n0_(*particles_->registerSharedVariable<Vecd>("InitialNormalDirection")) {}
//=============================================================================================//
void NormalDirectionFromShapeAndOp::update(size_t index_i, Real dt)
{
    Vecd normal_direction = switch_sign_ * shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
}
//=================================================================================================//
} // namespace SPH
