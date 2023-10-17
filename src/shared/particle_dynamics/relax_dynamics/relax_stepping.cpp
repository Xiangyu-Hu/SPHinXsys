#include "relax_stepping.hpp"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
void RelaxationResidue<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        residue -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
    }
    residue_[index_i] = residue;
};
//=================================================================================================//
RelaxationResidue<Inner<LevelSetCorrection>>::
    RelaxationResidue(BaseInnerRelation &inner_relation)
    : RelaxationResidue<Inner<>>(inner_relation), pos_(particles_->pos_),
      level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body_.body_shape_)) {}
//=================================================================================================//
RelaxationResidue<Inner<LevelSetCorrection>>::
    RelaxationResidue(BaseInnerRelation &inner_relation, std::string shape_name)
    : RelaxationResidue<Inner<>>(inner_relation), pos_(particles_->pos_)
{
    ComplexShape &complex_shape = DynamicCast<ComplexShape>(this, *sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
}
//=================================================================================================//
void RelaxationResidue<Inner<LevelSetCorrection>>::interaction(size_t index_i, Real dt)
{
    RelaxationResidue<Inner<>>::interaction(index_i, dt);
    residue_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
                                   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
}
//=================================================================================================//
void RelaxationResidue<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            residue -= 2.0 * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
    }
    residue_[index_i] += residue;
}
//=================================================================================================//
RelaxationScaling::RelaxationScaling(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      RelaxDataDelegateSimple(sph_body),
      residue_(*particles_->getVariableByName<Vecd>("ZeroOrderResidue")),
      h_ref_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
Real RelaxationScaling::reduce(size_t index_i, Real dt)
{
    return residue_[index_i].norm();
}
//=================================================================================================//
Real RelaxationScaling::outputResult(Real reduced_value)
{
    return 0.0625 * h_ref_ / (reduced_value + TinyReal);
}
//=================================================================================================//
PositionRelaxation::PositionRelaxation(SPHBody &sph_body)
    : LocalDynamics(sph_body), RelaxDataDelegateSimple(sph_body),
      sph_adaptation_(sph_body.sph_adaptation_), pos_(particles_->pos_),
      residue_(*particles_->getVariableByName<Vecd>("ZeroOrderResidue")) {}
//=================================================================================================//
void PositionRelaxation::update(size_t index_i, Real dt_square)
{
    pos_[index_i] += residue_[index_i] * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape)
    : LocalDynamics(sph_body), RelaxDataDelegateSimple(sph_body),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")),
      Vol_(particles_->Vol_), pos_(particles_->pos_), target_shape_(target_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, sph_body.sph_adaptation_)),
      reference_spacing_(particle_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::UpdateSmoothingLengthRatioByShape(SPHBody &sph_body)
    : UpdateSmoothingLengthRatioByShape(sph_body, *sph_body.body_shape_) {}
//=================================================================================================//
void UpdateSmoothingLengthRatioByShape::update(size_t index_i, Real dt_square)
{
    Real local_spacing = particle_adaptation_->getLocalSpacing(target_shape_, pos_[index_i]);
    h_ratio_[index_i] = reference_spacing_ / local_spacing;
    Vol_[index_i] = pow(local_spacing, Dimensions);
}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
