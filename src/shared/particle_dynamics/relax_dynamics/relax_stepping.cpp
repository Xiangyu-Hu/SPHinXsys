#include "relax_stepping.hpp"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
RelaxationResidue<Inner<>>::RelaxationResidue(BaseInnerRelation &inner_relation)
    : RelaxationResidue<Base, RelaxDataDelegateInner>(inner_relation),
      relax_shape_(sph_body_.getInitialShape()){};
//=================================================================================================//
RelaxationResidue<Inner<>>::
    RelaxationResidue(BaseInnerRelation &inner_relation, const std::string &sub_shape_name)
    : RelaxationResidue<Base, RelaxDataDelegateInner>(inner_relation),
      relax_shape_(*DynamicCast<ComplexShape>(this, sph_body_.getInitialShape())
                        .getSubShapeByName(sub_shape_name)) {}
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
RelaxationResidue<Inner<ComplexShapeBounding>>::
    RelaxationResidue(BaseInnerRelation& inner_relation, ComplexShape &complex_bounding_shapes)
    :RelaxationResidue<Base, RelaxDataDelegateInner>(inner_relation),pos_(particles_->pos_),
    complex_bounding_shapes_(complex_bounding_shapes)
{};
//=================================================================================================//
void RelaxationResidue<Inner<ComplexShapeBounding>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        residue -= 2.0 * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    }
    residue_[index_i] = residue;
    for (size_t i = 0; i != complex_bounding_shapes_.getLevelSetShapes().size(); ++i)
    {
       residue_[index_i] -=  2.0 * complex_bounding_shapes_.getLevelSetShapes()[i]->computeKernelGradientIntegral(pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));   
    };
}
//=================================================================================================//
void RelaxationResidue<Inner<LevelSetCorrection>>::interaction(size_t index_i, Real dt)
{
    RelaxationResidue<Inner<>>::interaction(index_i, dt);
    residue_[index_i] -= 2.0 * level_set_shape_.computeKernelGradientIntegral(
                                   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
}
//=================================================================================================//
void RelaxationResidue<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
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
    : UpdateSmoothingLengthRatioByShape(sph_body, sph_body.getInitialShape()) {}
//=================================================================================================//
void UpdateSmoothingLengthRatioByShape::update(size_t index_i, Real dt_square)
{
    Real local_spacing = particle_adaptation_->getLocalSpacing(target_shape_, pos_[index_i]);
    h_ratio_[index_i] = reference_spacing_ / local_spacing;
    Vol_[index_i] = pow(local_spacing, Dimensions);
}
//=================================================================================================//
RelaxationStepWithComplexBounding::RelaxationStepWithComplexBounding(BaseInnerRelation& inner_relation, ComplexShape& bounding_shapes)
    :BaseDynamics<void>(inner_relation.getSPHBody()),real_body_(DynamicCast<RealBody>(this, inner_relation.getSPHBody())),
      body_relations_(real_body_.getBodyRelations()),
      relaxation_residue_(inner_relation, bounding_shapes),
      relaxation_scaling_(real_body_),
      position_relaxation_(real_body_),
      near_shape_surface_(bounding_shapes),
      surface_bounding_(real_body_, near_shape_surface_) 
{}
//=================================================================================================//
void RelaxationStepWithComplexBounding::exec(Real dt)
{
    real_body_.updateCellLinkedList();
    for (size_t k = 0; k != body_relations_.size(); ++k)
    {
        body_relations_[k]->updateConfiguration();
    }
    relaxation_residue_.exec();
    Real scaling = relaxation_scaling_.exec();
    position_relaxation_.exec(scaling);
    surface_bounding_.exec();
}
} // namespace relax_dynamics
} // namespace SPH
