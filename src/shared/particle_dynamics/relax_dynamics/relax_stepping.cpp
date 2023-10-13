#include "relax_stepping.hpp"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      RelaxDataDelegateSimple(sph_body), acc_(particles_->acc_),
      h_ref_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
Real GetTimeStepSizeSquare::reduce(size_t index_i, Real dt)
{
    return acc_[index_i].norm();
}
//=================================================================================================//
Real GetTimeStepSizeSquare::outputResult(Real reduced_value)
{
    return 0.0625 * h_ref_ / (reduced_value + TinyReal);
}
//=================================================================================================//
void ParticleRelaxation<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        acceleration -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
    }
    pos_[index_i] += acceleration * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
};
//=================================================================================================//
ParticleRelaxation<Inner<LevelSetCorrection>>::
    ParticleRelaxation(BaseInnerRelation &inner_relation)
    : ParticleRelaxation<Inner<>>(inner_relation),
      level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body_.body_shape_)) {}
//=================================================================================================//
void ParticleRelaxation<Inner<LevelSetCorrection>>::interaction(size_t index_i, Real dt)
{
    ParticleRelaxation<Inner<>>::interaction(index_i, dt);
    Vecd acceleration = -2.0 * level_set_shape_->computeKernelGradientIntegral(
                                   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
    pos_[index_i] += acceleration * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
}
//=================================================================================================//
void interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    /** Contact interaction. */
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            acceleration -= 2.0 * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
    }

    acc_[index_i] = acceleration;
};

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
