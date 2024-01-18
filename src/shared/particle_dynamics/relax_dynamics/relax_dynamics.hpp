#include "relax_dynamics.h"
#include "base_particles.hpp"
#include "level_set_shape.h"
#include "particle_generator_lattice.h"
//========================================================================================================//
namespace SPH
{
//=====================================================================================================//
namespace relax_dynamics
{
//=================================================================================================//
template <class RelaxationType>
RelaxationAccelerationInner<RelaxationType>::RelaxationAccelerationInner(BaseInnerRelation& inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
    relaxation_type(), acc_(particles_->acc_), pos_(particles_->pos_),
    B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")) {}
//=================================================================================================//
template <class RelaxationType>
RelaxationAccelerationInnerWithLevelSetCorrection<RelaxationType>::
RelaxationAccelerationInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
    : RelaxationAccelerationInner<RelaxationType>(inner_relation), 
    sph_adaptation_(this->sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, this->sph_body_.body_shape_);
}
//=================================================================================================//
template <class RelaxationType>
RelaxationStepInner<RelaxationType>::
RelaxationStepInner(BaseInnerRelation& inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
    inner_relation_(inner_relation), near_shape_surface_(*real_body_),
    get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
    surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationInner<RelaxationType>>>(inner_relation);
    }
    else
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationInnerWithLevelSetCorrection<RelaxationType>>>(inner_relation);
    }
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationStepInner<RelaxationType>::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_bounding_.exec();
}
//=================================================================================================//
template <class RelaxationType>
RelaxationAccelerationComplex<RelaxationType>::
RelaxationAccelerationComplex(ComplexRelation& complex_relation)
    : LocalDynamics(complex_relation.getSPHBody()),
    RelaxDataDelegateComplex(complex_relation),
    relaxation_type(), acc_(particles_->acc_), pos_(particles_->pos_),
    B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix"))
{
    contact_B_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_B_[k] = this->contact_particles_[k]->template registerSharedVariable<Matd>("KernelCorrectionMatrix");
    }
}
//=================================================================================================//
template <class RelaxationType>
RelaxationAccelerationComplexWithLevelSetCorrection<RelaxationType>::
RelaxationAccelerationComplexWithLevelSetCorrection(ComplexRelation& complex_relation, const std::string& shape_name)
    : RelaxationAccelerationComplex<RelaxationType>(complex_relation), 
    sph_adaptation_(this->sph_body_.sph_adaptation_)
{
    ComplexShape& complex_shape = DynamicCast<ComplexShape>(this, *this->sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
}
//=================================================================================================//
template <class RelaxationType>
RelaxationStepComplex<RelaxationType>::
RelaxationStepComplex(ComplexRelation& complex_relation, const std::string& shape_name, bool level_set_correction)
    : BaseDynamics<void>(complex_relation.getSPHBody()), real_body_(complex_relation.getInnerRelation().real_body_),
    complex_relation_(complex_relation), near_shape_surface_(*real_body_, shape_name),
    get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
    surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationComplex<RelaxationType>>>(complex_relation);
    }
    else
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationComplexWithLevelSetCorrection<RelaxationType>>>(complex_relation, shape_name);
    }
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationStepComplex<RelaxationType>::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    complex_relation_.updateConfiguration();
    relaxation_acceleration_complex_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    //surface_bounding_.exec();
}
//=================================================================================================//
template <class RelaxationType>
RelaxationInnerImplicit<RelaxationType>::
    RelaxationInnerImplicit(BaseInnerRelation& inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
    kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
    Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
    B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_), relaxation_type()
{
    particles_->registerVariable(implicit_residual_, "ImplicitResidual");
    particles_->addVariableToWrite<Real>("ImplicitResidual");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
};
//=================================================================================================//
template <class RelaxationType>
ErrorAndParameters<Vecd, Matd, Matd> RelaxationInnerImplicit<RelaxationType>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) *
                           inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                           kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * 
                           Vol_[index_j] * dt * dt;

        error_and_parameters.error_ += relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) * 
                                       inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] * dt * dt;
        error_and_parameters.a_ -= parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }

    Matd evolution = Matd::Identity();
    error_and_parameters.a_ -= evolution;
    return error_and_parameters;
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationInnerImplicit<RelaxationType>::updateStates(size_t index_i, Real dt,
    const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
    Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

    pos_[index_i] += error_and_parameters.a_ * parameter_k;

    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) * 
                           inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                           kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * 
                           Vol_[index_j] * dt * dt;
        pos_[index_j] -= parameter_b * parameter_k;
    }
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationInnerImplicit<RelaxationType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
    acc_[index_i] = -error_and_parameters.error_ / dt / dt;
    implicit_residual_[index_i] = (error_and_parameters.error_ / dt / dt).norm();
}
//=================================================================================================//
template <class RelaxationType>
RelaxationInnerWithLevelSetCorrectionImplicit<RelaxationType>::
RelaxationInnerWithLevelSetCorrectionImplicit(BaseInnerRelation& inner_relation)
    : RelaxationInnerImplicit<RelaxationType>(inner_relation) {}
//=================================================================================================//
template <class RelaxationType>
ErrorAndParameters<Vecd, Matd, Matd> RelaxationInnerWithLevelSetCorrectionImplicit<RelaxationType>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = 
        RelaxationInnerImplicit<RelaxationType>::computeErrorAndParameters(index_i, dt);
    Real overlap = this->level_set_shape_->computeKernelIntegral(this->pos_[index_i], 
                   this->sph_adaptation_->SmoothingLengthRatio(index_i));

    error_and_parameters.error_ += this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) *
        this->level_set_shape_->computeKernelGradientIntegral(this->pos_[index_i],
            this->sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);
    error_and_parameters.a_ -= this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) *
        this->level_set_shape_->computeKernelSecondGradientIntegral(this->pos_[index_i],
            this->sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);

    return error_and_parameters;
}
//=================================================================================================//
template <class RelaxationType>
RelaxationStepInnerImplicit<RelaxationType>::
RelaxationStepInnerImplicit(BaseInnerRelation& inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()), time_step_size_(0.01),
    real_body_(inner_relation.real_body_), inner_relation_(inner_relation), 
    near_shape_surface_(*real_body_), get_time_step_(*real_body_),
    relaxation_evolution_inner_(inner_relation), surface_bounding_(near_shape_surface_) {}
//=================================================================================================//
template <class RelaxationType>
void RelaxationStepInnerImplicit<RelaxationType>::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    time_step_size_ =  sqrt(get_time_step_.exec());
    relaxation_evolution_inner_.exec(dt * time_step_size_);
}
//=================================================================================================//
template <class RelaxationType>
RelaxationComplexImplicit<RelaxationType>::
RelaxationComplexImplicit(ComplexRelation &complex_relation, const std::string& shape_name)
    :LocalDynamics(complex_relation.getSPHBody()), RelaxDataDelegateComplex(complex_relation),
    kernel_(complex_relation.getSPHBody().sph_adaptation_->getKernel()),
    Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
    B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")),
    implicit_residual_(*particles_->registerSharedVariable<Real>("ImplicitResidual")),
    sph_adaptation_(sph_body_.sph_adaptation_), relaxation_type()
{
    contact_B_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
        contact_B_[k] = this->contact_particles_[k]->template registerSharedVariable<Matd>("KernelCorrectionMatrix");
    }
    particles_->addVariableToWrite<Real>("ImplicitResidual");

    ComplexShape& complex_shape = DynamicCast<ComplexShape>(this, *this->sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
};
//=================================================================================================//
template <class RelaxationType>
ErrorAndParameters<Vecd, Matd, Matd> RelaxationComplexImplicit<RelaxationType>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;

    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) *
                           inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                           kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * 
                           Vol_[index_j] * dt * dt;

        error_and_parameters.error_ += relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) *
                                       inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] * dt * dt;
        error_and_parameters.a_ -= parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
        StdLargeVec<Matd>& B_k = *(contact_B_[k]);
        Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_k[index_j]) *
                               contact_neighborhood.e_ij_[n] * contact_neighborhood.e_ij_[n].transpose() *
                               kernel_->d2W(contact_neighborhood.r_ij_[n], contact_neighborhood.e_ij_[n]) * 
                               Vol_k[index_j] * dt * dt;

            error_and_parameters.error_ += relaxation_type.getBackgroundForce(B_[index_i], B_k[index_j]) *
                                           contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n] * dt * dt;
            error_and_parameters.a_ -= parameter_b;

            /* With the wall*/
           /* size_t index_j = contact_neighborhood.j_[n];
            Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_[index_i]) *
                contact_neighborhood.e_ij_[n] * contact_neighborhood.e_ij_[n].transpose() *
                kernel_->d2W(contact_neighborhood.r_ij_[n], contact_neighborhood.e_ij_[n]) *
                Vol_k[index_j] * dt * dt;

            error_and_parameters.error_ += relaxation_type.getBackgroundForce(B_[index_i], B_[index_i]) *
                contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n] * dt * dt;
            error_and_parameters.a_ -= parameter_b;*/
        }
    }

    Matd evolution = Matd::Identity();
    error_and_parameters.a_ -= evolution;
    return error_and_parameters;
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationComplexImplicit<RelaxationType>::updateStates(size_t index_i, Real dt,
    const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
    Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

    pos_[index_i] += error_and_parameters.a_ * parameter_k;

    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Matd parameter_b = relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) *
                           inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                           kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * 
                           Vol_[index_j] * dt * dt;
        pos_[index_j] -= parameter_b * parameter_k;
    }
}
//=================================================================================================//
template <class RelaxationType>
void RelaxationComplexImplicit<RelaxationType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
    acc_[index_i] = -error_and_parameters.error_ / dt / dt;
    implicit_residual_[index_i] = (error_and_parameters.error_ / dt / dt).norm();
}
//=================================================================================================//
template <class RelaxationType>
RelaxationComplexWithLevelSetCorrectionImplicit<RelaxationType>::
RelaxationComplexWithLevelSetCorrectionImplicit(ComplexRelation& complex_relation, const std::string& shape_name)
    : RelaxationComplexImplicit<RelaxationType>(complex_relation, shape_name) {}
//=================================================================================================//
template <class RelaxationType>
ErrorAndParameters<Vecd, Matd, Matd> RelaxationComplexWithLevelSetCorrectionImplicit<RelaxationType>::
computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters =
        RelaxationComplexImplicit<RelaxationType>::computeErrorAndParameters(index_i, dt);
    Real overlap = this->level_set_shape_->computeKernelIntegral(this->pos_[index_i], 
                   this->sph_adaptation_->SmoothingLengthRatio(index_i));

    error_and_parameters.error_ += this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) *
                                   this->level_set_shape_->computeKernelGradientIntegral(this->pos_[index_i],
                                   this->sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);
    error_and_parameters.a_ -= this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) *
                               this->level_set_shape_->computeKernelSecondGradientIntegral(this->pos_[index_i],
                               this->sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt * (1 + overlap);

    return error_and_parameters;
}
//=================================================================================================//
template <class RelaxationType>
RelaxationStepComplexImplicit<RelaxationType>::
RelaxationStepComplexImplicit(ComplexRelation& complex_relation, const std::string& shape_name, 
                              bool level_set_correction)
    : BaseDynamics<void>(complex_relation.getSPHBody()), time_step_size_(0.01),
    real_body_(complex_relation.getInnerRelation().real_body_), complex_relation_(complex_relation), 
    near_shape_surface_(*real_body_, shape_name), get_time_step_(*real_body_),
    relaxation_evolution_complex_(complex_relation, shape_name), surface_bounding_(near_shape_surface_) {}
//=================================================================================================//
template <class RelaxationType>
void RelaxationStepComplexImplicit<RelaxationType>::exec(Real dt)
{
    time_step_size_ = sqrt(get_time_step_.exec());
    relaxation_evolution_complex_.exec(dt * time_step_size_);
    real_body_->updateCellLinkedList();
    complex_relation_.updateConfiguration();
}
//=================================================================================================//
} // namespace relax_dynamics
  //=================================================================================================//

} // namespace SPH
//=================================================================================================//