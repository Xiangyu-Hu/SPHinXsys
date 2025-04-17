#include "relax_stepping.hpp"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
RelaxationResidue<Inner<>>::RelaxationResidue(BaseInnerRelation &inner_relation)
    : RelaxationResidue<Base, DataDelegateInner>(inner_relation),
      relax_shape_(sph_body_.getInitialShape()){};
//=================================================================================================//
RelaxationResidue<Inner<>>::
    RelaxationResidue(BaseInnerRelation &inner_relation, const std::string &sub_shape_name)
    : RelaxationResidue<Base, DataDelegateInner>(inner_relation),
      relax_shape_(*DynamicCast<ComplexShape>(this, sph_body_.getInitialShape()).getSubShapeByName(sub_shape_name)) {}
//=================================================================================================//
void RelaxationResidue<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        residue -= 2.0 * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    }
    residue_[index_i] = residue;
};
//=================================================================================================//
void RelaxationResidue<Inner<LevelSetCorrection>>::interaction(size_t index_i, Real dt)
{
    RelaxationResidue<Inner<>>::interaction(index_i, dt);
    residue_[index_i] -= 2.0 * level_set_shape_.computeKernelGradientIntegral(
                                   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
}
//=================================================================================================//
void RelaxationResidue<Inner<Implicit>>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
    residue_[index_i] = -error_and_parameters.error_ / dt / dt;
    kinetic_energy_[index_i] = residue_[index_i].norm();
}
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationResidue<Inner<Implicit>>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = 2.0 * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
            kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) *
            Vol_[index_j] * dt * dt;

        error_and_parameters.error_ += 2.0 * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * 
            inner_neighborhood.e_ij_[n] * dt * dt;
        error_and_parameters.a_ -= parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }

    Matd evolution = Matd::Identity();
    error_and_parameters.a_ -= evolution;
    return error_and_parameters;
}
//=================================================================================================//
void RelaxationResidue<Inner<Implicit>>::updateStates(size_t index_i, Real dt,
    const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
    Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

    pos_[index_i] += error_and_parameters.a_ * parameter_k;

    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = 2.0 * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
            kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;
        pos_[index_j] -= parameter_b * parameter_k;
    }
}
//=================================================================================================//
void RelaxationResidue<Inner<LevelSetCorrection, Implicit>>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
    residue_[index_i] = -error_and_parameters.error_ / dt / dt;
    kinetic_energy_[index_i] = residue_[index_i].norm();
}
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationResidue<Inner<LevelSetCorrection, Implicit>>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters =
        RelaxationResidue<Inner<Implicit>>::computeErrorAndParameters(index_i, dt);
    Real overlap = level_set_shape_.computeKernelIntegral(pos_[index_i],
        sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;

    error_and_parameters.error_ += 2.0 * level_set_shape_.computeKernelGradientIntegral(pos_[index_i],
            sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);
    error_and_parameters.a_ -= 2.0 * level_set_shape_.computeKernelSecondGradientIntegral(pos_[index_i],
            sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);

    return error_and_parameters;
}
//=================================================================================================//
void RelaxationResidue<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            residue -= 2.0 * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
        }
    }
    residue_[index_i] += residue;
}
//=================================================================================================//
RelaxationScaling::RelaxationScaling(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      residue_(particles_->getVariableDataByName<Vecd>("ZeroOrderResidue")),
      h_ref_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()) {}
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
    : LocalDynamics(sph_body),
      sph_adaptation_(&sph_body.getSPHAdaptation()),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      residue_(particles_->getVariableDataByName<Vecd>("ZeroOrderResidue")) {}
//=================================================================================================//
void PositionRelaxation::update(size_t index_i, Real dt_square)
{
    pos_[index_i] += residue_[index_i] * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape)
    : LocalDynamics(sph_body),
      h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      target_shape_(target_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, &sph_body.getSPHAdaptation())),
      reference_spacing_(particle_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::UpdateSmoothingLengthRatioByShape(SPHBody &sph_body)
    : UpdateSmoothingLengthRatioByShape(sph_body, sph_body.getInitialShape()) {}
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
