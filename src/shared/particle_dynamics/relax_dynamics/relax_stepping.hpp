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
RelaxationResidual<Base, DataDelegationType>::RelaxationResidual(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      sph_adaptation_(&this->sph_body_.getSPHAdaptation()),
      kernel_(base_relation.getSPHBody().getSPHAdaptation().getKernel()),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      kinetic_energy_(this->particles_->template registerStateVariableData<Real>("ParticleKineticEnergy")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      residual_(this->particles_->template registerStateVariableData<Vecd>("ZeroOrderResidual")) {}
//=================================================================================================//
template <typename... Args>
RelaxationResidual<Inner<LevelSetCorrection>>::RelaxationResidual(Args &&...args)
    : RelaxationResidual<Inner<>>(std::forward<Args>(args)...),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      level_set_shape_(DynamicCast<LevelSetShape>(this, this->getRelaxShape())){};
//=================================================================================================//
template <typename... Args>
RelaxationResidual<Inner<Implicit>>::RelaxationResidual(Args &&...args)
    : RelaxationResidual<Inner<>>(std::forward<Args>(args)...){};
//=================================================================================================//
template <typename... Args>
RelaxationResidual<Inner<LevelSetCorrection, Implicit>>::RelaxationResidual(Args &&...args)
    : RelaxationResidual<Inner<Implicit>>(std::forward<Args>(args)...),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      level_set_shape_(DynamicCast<LevelSetShape>(this, this->getRelaxShape())){};
//=================================================================================================//
template <class RelaxationResidualType>
template <typename FirstArg, typename... OtherArgs>
RelaxationStep<RelaxationResidualType>::
    RelaxationStep(FirstArg &&first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      real_body_(DynamicCast<RealBody>(this, first_arg.getSPHBody())),
      body_relations_(real_body_.getBodyRelations()),
      relaxation_residual_(first_arg, std::forward<OtherArgs>(other_args)...),
      relaxation_scaling_(real_body_), position_relaxation_(real_body_),
      near_shape_surface_(real_body_, DynamicCast<LevelSetShape>(this, relaxation_residual_.getRelaxShape())),
      surface_bounding_(near_shape_surface_) {}
//=================================================================================================//
template <class RelaxationResidualType>
void RelaxationStep<RelaxationResidualType>::exec(Real dt)
{
    real_body_.updateCellLinkedList();
    for (size_t k = 0; k != body_relations_.size(); ++k)
    {
        body_relations_[k]->updateConfiguration();
    }
    relaxation_residual_.exec();
    Real scaling = relaxation_scaling_.exec();
    position_relaxation_.exec(scaling);
    surface_bounding_.exec();
}
//=================================================================================================//
template <class RelaxationResidualType>
template <typename FirstArg, typename... OtherArgs>
RelaxationStepImplicit<RelaxationResidualType>::
    RelaxationStepImplicit(FirstArg &&first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      real_body_(DynamicCast<RealBody>(this, first_arg.getSPHBody())),
      body_relations_(real_body_.getBodyRelations()),
      relaxation_residual_(first_arg, std::forward<OtherArgs>(other_args)...),
      relaxation_scaling_(real_body_), position_relaxation_(real_body_),
      near_shape_surface_(real_body_, DynamicCast<LevelSetShape>(this, relaxation_residual_.getRelaxShape())),
      surface_bounding_(near_shape_surface_)
{
}
//=================================================================================================//
template <class RelaxationResidualType>
void RelaxationStepImplicit<RelaxationResidualType>::exec(Real dt)
{
    Real scaling = SMIN(Real(sqrt(relaxation_scaling_.exec())), Real(0.01));
    relaxation_residual_.exec(scaling);
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
