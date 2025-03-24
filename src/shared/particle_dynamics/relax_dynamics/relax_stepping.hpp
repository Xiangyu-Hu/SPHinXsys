#ifndef RELAX_STEPPING_HPP
#define RELAX_STEPPING_HPP

#include "relax_stepping.h"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
RelaxationResidue<Base, DataDelegationType>::RelaxationResidue(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      sph_adaptation_(&this->sph_body_.getSPHAdaptation()),
      kernel_(base_relation.getSPHBody().getSPHAdaptation().getKernel()),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      kinetic_energy_(this->particles_->template registerStateVariable<Real>("ParticleKineticEnergy")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      residue_(this->particles_->template registerStateVariable<Vecd>("ZeroOrderResidue")){}
//=================================================================================================//
template <typename... Args>
RelaxationResidue<Inner<LevelSetCorrection>>::RelaxationResidue(Args &&...args)
    : RelaxationResidue<Inner<>>(std::forward<Args>(args)...),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      level_set_shape_(DynamicCast<LevelSetShape>(this, this->getRelaxShape())){};
//=================================================================================================//
template <typename... Args>
RelaxationResidue<Inner<Implicit>>::RelaxationResidue(Args &&...args)
    : RelaxationResidue<Inner<>>(std::forward<Args>(args)...) {};
//=================================================================================================//
template <typename... Args>
RelaxationResidue<Inner<LevelSetCorrection, Implicit>>::RelaxationResidue(Args &&...args)
    : RelaxationResidue<Inner<Implicit>>(std::forward<Args>(args)...),
    pos_(particles_->getVariableDataByName<Vecd>("Position")),
    level_set_shape_(DynamicCast<LevelSetShape>(this, this->getRelaxShape())) {};
//=================================================================================================//
template <class RelaxationResidueType>
template <typename FirstArg, typename... OtherArgs>
RelaxationStep<RelaxationResidueType>::
    RelaxationStep(FirstArg &&first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      real_body_(DynamicCast<RealBody>(this, first_arg.getSPHBody())),
      body_relations_(real_body_.getBodyRelations()),
      relaxation_residue_(first_arg, std::forward<OtherArgs>(other_args)...),
      relaxation_scaling_(real_body_), position_relaxation_(real_body_),
      near_shape_surface_(real_body_, DynamicCast<LevelSetShape>(this, relaxation_residue_.getRelaxShape())),
      surface_bounding_(near_shape_surface_) {}
//=================================================================================================//
template <class RelaxationResidueType>
void RelaxationStep<RelaxationResidueType>::exec(Real dt)
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
//=================================================================================================//
template <class RelaxationResidueType>
template <typename FirstArg, typename... OtherArgs>
RelaxationStepImplicit<RelaxationResidueType>::
    RelaxationStepImplicit(FirstArg&& first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
    real_body_(DynamicCast<RealBody>(this, first_arg.getSPHBody())),
    body_relations_(real_body_.getBodyRelations()),
    relaxation_residue_(first_arg, std::forward<OtherArgs>(other_args)...),
    relaxation_scaling_(real_body_), position_relaxation_(real_body_),
    near_shape_surface_(real_body_, DynamicCast<LevelSetShape>(this, relaxation_residue_.getRelaxShape())),
    surface_bounding_(near_shape_surface_) {
}
//=================================================================================================//
template <class RelaxationResidueType>
void RelaxationStepImplicit<RelaxationResidueType>::exec(Real dt)
{
    Real scaling = SMIN(Real(sqrt(relaxation_scaling_.exec())), Real(0.01));
    relaxation_residue_.exec(scaling);
    real_body_.updateCellLinkedList();
    for (size_t k = 0; k != body_relations_.size(); ++k)
    {
        body_relations_[k]->updateConfiguration();
    }
    surface_bounding_.exec();
}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_STEPPING_HPP
