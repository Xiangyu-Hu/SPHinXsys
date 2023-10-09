/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef RELAX_DYNAMICS_H
#define RELAX_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "cell_linked_list.h"
#include "general_dynamics.h"

namespace SPH
{
class GeometryShape;
class LevelSetShape;

namespace relax_dynamics
{
typedef DataDelegateSimple<BaseParticles> RelaxDataDelegateSimple;

typedef DataDelegateInner<BaseParticles> RelaxDataDelegateInner;

typedef DataDelegateComplex<BaseParticles, BaseParticles> RelaxDataDelegateComplex;

/**
 * @class GetTimeStepSizeSquare
 * @brief relaxation dynamics for particle initialization
 * computing the square of time step size
 */
class GetTimeStepSizeSquare : public LocalDynamicsReduce<Real, ReduceMax>,
                              public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &acc_;
    Real h_ref_;

  public:
    explicit GetTimeStepSizeSquare(SPHBody &sph_body);
    virtual ~GetTimeStepSizeSquare(){};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value);
};

/**
 * @class RelaxationAccelerationInner
 * @brief simple algorithm for physics relaxation
 * without considering contact interaction.
 * this is usually used for solid like bodies
 */
class RelaxationAccelerationInner : public LocalDynamics, 
                                    public RelaxDataDelegateInner
{
  public:
    explicit RelaxationAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~RelaxationAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            acceleration -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
        acc_[index_i] = acceleration;
    };

  protected:
    StdLargeVec<Vecd> &acc_, &pos_;
};

/**
 * @class RelaxationAccelerationInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationAccelerationInnerWithLevelSetCorrection : public RelaxationAccelerationInner
{
  public:
    explicit RelaxationAccelerationInnerWithLevelSetCorrection(BaseInnerRelation &inner_relation);
    virtual ~RelaxationAccelerationInnerWithLevelSetCorrection(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationInner::interaction(index_i, dt);
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));

        //if (phi > -constrained_distance_)
        //{
        //    acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i], 
        //        sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        //}
        //else
        //{
        //    acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i], 
        //        sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 - overlap);
        //};
        acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i],
                sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
    };

  protected:
    LevelSetShape *level_set_shape_;
    SPHAdaptation *sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class UpdateParticlePosition
 * @brief update the particle position for a time step
 */
class UpdateParticlePosition : public LocalDynamics,
                               public RelaxDataDelegateSimple
{
  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &pos_, &acc_;

  public:
    explicit UpdateParticlePosition(SPHBody &sph_body);
    virtual ~UpdateParticlePosition(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class UpdateSmoothingLengthRatioByShape
 * @brief update the particle smoothing length ratio
 */
class UpdateSmoothingLengthRatioByShape : public LocalDynamics,
                                          public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Real> &h_ratio_, &Vol_;
    StdLargeVec<Vecd> &pos_;
    Shape &target_shape_;
    ParticleRefinementByShape *particle_adaptation_;
    Real reference_spacing_;

  public:
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape);
    explicit UpdateSmoothingLengthRatioByShape(SPHBody &sph_body);
    virtual ~UpdateSmoothingLengthRatioByShape(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class RelaxationAccelerationComplex
 * @brief compute relaxation acceleration while consider the present of contact bodies
 * with considering contact interaction
 * this is usually used for fluid like bodies //TODO: seems better called as Contact
 */
class RelaxationAccelerationComplex : public LocalDynamics,
                                      public RelaxDataDelegateComplex
{
  public:
    explicit RelaxationAccelerationComplex(ComplexRelation &complex_relation);
    virtual ~RelaxationAccelerationComplex(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            acceleration -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }

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

  protected:
    StdLargeVec<Vecd> &acc_, &pos_;
};

/**
 * @class ShapeSurfaceBounding
 * @brief constrain surface particles by
 * map constrained particles to geometry face and
 * r = r + phi * norm (vector distance to face)
 */
class ShapeSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                             public RelaxDataDelegateSimple
{
  public:
    ShapeSurfaceBounding(NearShapeSurface &body_part);
    virtual ~ShapeSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    LevelSetShape *level_set_shape_;
    Real constrained_distance_;
};

/**
 * @class NearSurfaceVolumeCorrection
 * @brief
 */
class NearSurfaceVolumeCorrection : public BaseLocalDynamics<BodyPartByCell>,
    public RelaxDataDelegateSimple
{
public:
    NearSurfaceVolumeCorrection(NearShapeSurface& body_part);
    virtual ~NearSurfaceVolumeCorrection() {};
    void update(size_t index_i, Real dt = 0.0);

protected:
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Real>& Vol_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class RelaxationStepInner
 * @brief carry out particle relaxation step of particles within the body
 */
class RelaxationStepInner : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepInner(BaseInnerRelation &inner_relation,
                                 bool level_set_correction = false);
    virtual ~RelaxationStepInner(){};
    SimpleDynamics<ShapeSurfaceBounding> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    BaseInnerRelation &inner_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_inner_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
    SharedPtr<BaseDynamics<void>> surface_correction_;
};

/**
 * @class RelaxationAccelerationByStressInner
 * @brief simple algorithm for physics relaxation by stress
 * without considering contact interaction.
 * this is usually used for solid like bodies.
 */
class RelaxationAccelerationByStressInner : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationAccelerationByStressInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationAccelerationByStressInner() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            acceleration -= (B_[index_i] + B_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
        acc_[index_i] = acceleration;
    };

protected:
    StdLargeVec<Vecd>& acc_, & pos_;
    StdLargeVec<Matd>& B_;
};

/**
 * @class RelaxationAccelerationByStressInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationAccelerationByStressInnerWithLevelSetCorrection : public RelaxationAccelerationByStressInner
{
public:
    explicit RelaxationAccelerationByStressInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationAccelerationByStressInnerWithLevelSetCorrection() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationByStressInner::interaction(index_i, dt);
        acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
                          pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
    };

protected:
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class RelaxationStepByStressInner
 * @brief carry out particle relaxation step of particle by stress
 *  within the body with the first order consisitency.
 */
class RelaxationStepByStressInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepByStressInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepByStressInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_inner_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
    SharedPtr<BaseDynamics<void>> surface_correction_;
};

/**
 * @class RelaxationAccelerationComplexWithLevelSetCorrection
 * @brief compute relaxation acceleration while consider the present of contact bodies
 * with considering contact interaction
 * this is usually used for fluid like bodies
 * we constrain particles with a level-set correction function when the fluid boundary is not contacted with solid.
 */
class RelaxationAccelerationComplexWithLevelSetCorrection : public RelaxationAccelerationComplex
{
  public:
    RelaxationAccelerationComplexWithLevelSetCorrection(
        ComplexRelation &complex_relation, const std::string &shape_name);
    virtual ~RelaxationAccelerationComplexWithLevelSetCorrection(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationComplex::interaction(index_i, dt);
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));
        if (phi > -constrained_distance_)
        {
            acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
                         pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        }
        else
        {
            acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
                         pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 - overlap);
        };
    };

  protected:
    LevelSetShape *level_set_shape_;
    SPHAdaptation *sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class RelaxationStepComplex
 * @brief carry out particle relaxation step of particles within multi bodies
 */
class RelaxationStepComplex : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepComplex(ComplexRelation &complex_relation,
                                   const std::string &shape_name, bool level_set_correction = false);
    virtual ~RelaxationStepComplex(){};
    SimpleDynamics<ShapeSurfaceBounding> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    ComplexRelation &complex_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_complex_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

template <typename ErrorDataType, typename ParameterADataType, typename ParameterCDataType>
struct ErrorAndParameters
{
    ErrorDataType error_;
    ParameterADataType a_;
    ParameterCDataType c_;

    ErrorAndParameters<ErrorDataType, ParameterADataType, ParameterCDataType>() : 
        error_(ZeroData<ErrorDataType>::value),
        a_(ZeroData<ParameterADataType>::value),
        c_(ZeroData<ParameterCDataType>::value) {};
};

/**
 * @class RelaxationImplicitInner
 * @brief carry out particle relaxation by position with implicit evolution.
 */
class RelaxationImplicitInner : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationImplicitInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationImplicitInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Real target_error_p_; //It is generally to be the maximum error in the whole domain.
    StdLargeVec<Real> error_p_; //It contains the intermediate error.

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Real> implicit_residue_p_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;

public:
    inline void updateTargetError(Real target_error) { Real target_error_p_ = target_error; }
};

/**
 * @class RelaxationImplicitInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationImplicitInnerWithLevelSetCorrection : public RelaxationImplicitInner
{
public:
    explicit RelaxationImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationImplicitInnerWithLevelSetCorrection() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepImplicitInner
 * @brief carry out the particle relaxation evolution within the body
 */
class RelaxationStepImplicitInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepImplicitInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real target_error_p_;

    Real time_step_size_;
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationImplicitInner> relaxation_evolution_inner_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
    SimpleDynamics<NearSurfaceVolumeCorrection> surface_correction_;

    ReduceDynamics<QuantityMaximum<Real>> update_averaged_error_;
};

/**
 * @class RelaxationByStressImplicitInner
 * @brief carry out particle relaxation by first ordre consistency implicit evolution.
 */
class RelaxationByStressImplicitInner : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationByStressImplicitInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationByStressImplicitInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Real target_error_s_;
    StdLargeVec<Real> error_s_;

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> implicit_residue_s_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;

public:
    inline void updateTargetError(Real target_error) { target_error_s_ = target_error; }
};

/**
 * @class RelaxationByStressEvolutionInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationByStressImplicitInnerWithLevelSetCorrection : public RelaxationByStressImplicitInner
{
public:
    explicit RelaxationByStressImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationByStressImplicitInnerWithLevelSetCorrection() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepByStressInner
 * @brief carry out the particle relaxation evolution from first order consistency within the body
 */
class RelaxationStepByStressImplicitInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepByStressImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepByStressImplicitInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real target_error_s_;

    Real time_step_size_;
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationByStressImplicitInner> relaxation_evolution_inner_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
    SimpleDynamics<NearSurfaceVolumeCorrection> surface_correction_;

    ReduceDynamics<QuantityMaximum<Real>> update_averaged_error_;
};

/**
 * @class CalcualteParticleStress
 * @brief calculate the particle stress with first order consistency
 */
class CalculateParticleStress : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit CalculateParticleStress(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CalculateParticleStress() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    StdLargeVec<Vecd> pos_;
    StdLargeVec<Matd> B_;
    StdLargeVec<Matd> stress_;
    bool level_set_correction_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class UpdateParticleKineticEnergy
 * @brief calculate the particle kinetic energy
 */
class UpdateParticleKineticEnergy : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation);
    virtual ~UpdateParticleKineticEnergy() {};
    void interaction(size_t index_i, Real dt);

protected:
    StdLargeVec<Real>& mass_;
    StdLargeVec<Vecd>& acc_;
    StdLargeVec<Real> particle_kinetic_energy;
};


/**
 * @class CheckCorrectedZeroOrderConsistency
 * @brief calculate the corrected zero order consistency
 */
class CheckCorrectedZeroOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedZeroOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedZeroOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_zero_order_error_;
    StdLargeVec<Vecd> corrected_zero_order_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class CheckCorrectedFirstOrderConsistency
 * @brief calculate the corrected first order consistency
 */
class CheckCorrectedFirstOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedFirstOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_first_order_error_;
    StdLargeVec<Matd> corrected_first_order_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class CheckConsistencyRealization
 * @brief check the consistency of SPH conservative formulation.
 */
class CheckConsistencyRealization : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckConsistencyRealization() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& pressure_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> pressure_gradient_norm_;
    StdLargeVec<Vecd> pressure_gradient_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class RelaxationAccelerationByCMInner
 * @brief simple algorithm for physics relaxation by correction 
 *        matrix without considering contact interaction.
 */
class RelaxationAccelerationByCMInner : public LocalDynamics, 
                                        public RelaxDataDelegateInner
{
public:
    explicit RelaxationAccelerationByCMInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationAccelerationByCMInner() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            acceleration -= (B_[index_i] + B_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
        acc_[index_i] = acceleration;
    };

protected:
    StdLargeVec<Vecd>& acc_, & pos_;
    StdLargeVec<Matd>& B_;
};

/**
 * @class RelaxationAccelerationByCMInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationAccelerationByCMInnerWithLevelSetCorrection : 
    public RelaxationAccelerationByCMInner
{
public:
    explicit RelaxationAccelerationByCMInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationAccelerationByCMInnerWithLevelSetCorrection() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationByCMInner::interaction(index_i, dt);
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
        //if (phi > -constrained_distance_)
        //{
        //    acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
        //                      pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        //}
        //else
        //{
        //    acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
        //                      pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 - overlap);
        //};

        acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
                          pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
    };

protected:
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class RelaxationStepByCMInner
 * @brief carry out particle relaxation step of particle by stress
 *  within the body with the first order consisitency.
 */
class RelaxationStepByCMInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepByCMInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepByCMInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_inner_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
 * @class RelaxationAccelerationByCMComplex
 */
class RelaxationAccelerationByCMComplex : public LocalDynamics,
    public RelaxDataDelegateComplex
{
public:
    explicit RelaxationAccelerationByCMComplex(ComplexRelation& complex_relation);
    virtual ~RelaxationAccelerationByCMComplex() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            acceleration -= (B_[index_i] + B_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }

        /** Contact interaction. */
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Matd>& B_k = *(contact_B_[k]);
            Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                acceleration -= (B_[index_i] + B_k[index_j]) * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }

        acc_[index_i] = acceleration;
    };

protected:
    StdLargeVec<Vecd>& acc_, & pos_;
    StdLargeVec<Matd>& B_;
    StdVec<StdLargeVec<Matd>*> contact_B_;
};

/**
 * @class RelaxationAccelerationByCMComplexWithLevelSetCorrection
 */
class RelaxationAccelerationByCMComplexWithLevelSetCorrection : 
    public RelaxationAccelerationByCMComplex
{
public:
    RelaxationAccelerationByCMComplexWithLevelSetCorrection(
        ComplexRelation& complex_relation, const std::string& shape_name);
    virtual ~RelaxationAccelerationByCMComplexWithLevelSetCorrection() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationByCMComplex::interaction(index_i, dt);
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));
        if (phi > -constrained_distance_)
        {
            acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
                         pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        }
        else
        {
            acc_[index_i] -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
                         pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 - overlap);
        };
    };

protected:
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class RelaxationStepByCMComplex
 */
class RelaxationStepByCMComplex : public BaseDynamics<void>
{
public:
    explicit RelaxationStepByCMComplex(ComplexRelation& complex_relation,
        const std::string& shape_name, bool level_set_correction = false);
    virtual ~RelaxationStepByCMComplex() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    RealBody* real_body_;
    ComplexRelation& complex_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_complex_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

//****************************************IMPLICIT SCHEME******************************************//
template <typename ErrorDataType, typename ParameterADataType, typename ParameterCDataType>
struct ErrorAndParameters
{
    ErrorDataType error_;
    ParameterADataType a_;
    ParameterCDataType c_;

    ErrorAndParameters<ErrorDataType, ParameterADataType, ParameterCDataType>() : 
        error_(ZeroData<ErrorDataType>::value),
        a_(ZeroData<ParameterADataType>::value),
        c_(ZeroData<ParameterCDataType>::value) {};
};

/**
 * @class RelaxationImplicitInner
 * @brief carry out particle relaxation by position with implicit evolution.
 */
class RelaxationImplicitInner : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationImplicitInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationImplicitInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Real> implicit_residual_pressure_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class RelaxationImplicitInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationImplicitInnerWithLevelSetCorrection : public RelaxationImplicitInner
{
public:
    explicit RelaxationImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationImplicitInnerWithLevelSetCorrection() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepImplicitInner
 * @brief carry out the particle relaxation evolution within the body
 */
class RelaxationStepImplicitInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepImplicitInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real time_step_size_;
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationImplicitInnerWithLevelSetCorrection> relaxation_evolution_inner_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
 * @class RelaxationByCMImplicitInner
 * @brief carry out particle relaxation by first ordre consistency implicit evolution.
 */
class RelaxationByCMImplicitInner : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationByCMImplicitInner(BaseInnerRelation& inner_relation);
    virtual ~RelaxationByCMImplicitInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> implicit_residual_cm_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class RelaxationByCMImplicitInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
class RelaxationByCMImplicitInnerWithLevelSetCorrection : public RelaxationByCMImplicitInner
{
public:
    explicit RelaxationByCMImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation);
    virtual ~RelaxationByCMImplicitInnerWithLevelSetCorrection() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepByStressInner
 * @brief carry out the particle relaxation evolution from first order consistency within the body
 */
class RelaxationStepByCMImplicitInner : public BaseDynamics<void>
{
public:
    explicit RelaxationStepByCMImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~RelaxationStepByCMImplicitInner() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real time_step_size_;
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationByCMImplicitInnerWithLevelSetCorrection> relaxation_evolution_inner_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
  * @class CorrectedConfigurationInnerWithLevelSet
  * @brief calculate the correction matrix based on the level set
  */
class CorrectedConfigurationInnerWithLevelSet : public LocalDynamics, public RelaxDataDelegateInner
{
  public:
    explicit CorrectedConfigurationInnerWithLevelSet(BaseInnerRelation &inner_relation, bool level_set_correction = false);
    virtual ~CorrectedConfigurationInnerWithLevelSet(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Matd> &B_with_level_set_;
    bool level_set_correction_;
    LevelSetShape *level_set_shape_;
    SPHAdaptation *sph_adaptation_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class UpdateParticleKineticEnergy
 * @brief calculate the particle kinetic energy
 */
class UpdateParticleKineticEnergy : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation);
    virtual ~UpdateParticleKineticEnergy() {};
    void interaction(size_t index_i, Real dt);

protected:
    StdLargeVec<Real>& mass_;
    StdLargeVec<Vecd>& acc_;
    StdLargeVec<Real> particle_kinetic_energy;
};

/**
 * @class CheckCorrectedZeroOrderConsistency
 * @brief calculate the corrected zero order consistency
 */
class CheckCorrectedZeroOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedZeroOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedZeroOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_zero_order_error_norm_;
    StdLargeVec<Vecd> corrected_zero_order_error_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;

};

/**
 * @class CheckCorrectedFirstOrderConsistency
 * @brief calculate the corrected first order consistency
 */
class CheckCorrectedFirstOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedFirstOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_first_order_error_norm_;
    StdLargeVec<Matd> corrected_first_order_error_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class CheckConsistencyRealization
 * @brief check the consistency of SPH conservative formulation.
 */
class CheckConsistencyRealization : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckConsistencyRealization() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& pressure_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> relaxation_error_;
    StdLargeVec<Real> pressure_gradient_error_norm_;
    StdLargeVec<Vecd> pressure_gradient_;
    StdLargeVec<Real> zero_order_error_norm_;
    StdLargeVec<Vecd> zero_order_error_;
    StdLargeVec<Real> reproduce_gradient_error_norm_;
    StdLargeVec<Vecd> reproduce_gradient_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class CheckReverseConsistencyRealization
 */
class CheckReverseConsistencyRealization : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckReverseConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckReverseConsistencyRealization() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& pressure_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> pressure_gradient_error_norm_;
    StdLargeVec<Vecd> pressure_gradient_;
    StdLargeVec<Real> zero_order_error_norm_;
    StdLargeVec<Vecd> zero_order_error_;
    StdLargeVec<Real> reproduce_gradient_error_norm_;
    StdLargeVec<Vecd> reproduce_gradient_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class ShellMidSurfaceBounding
 * @brief constrain particles by constraining particles to mid-surface.
 * Note that level_set_refinement_ratio should be smaller than particle_spacing_ref_ / (0.05 * thickness_)
 * because if level_set_refinement_ratio > particle_spacing_ref_ / (0.05 * thickness_),
 * there will be no level set field.
 */
class ShellMidSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                                public RelaxDataDelegateInner
{
  public:
    ShellMidSurfaceBounding(NearShapeSurface &body_part, BaseInnerRelation &inner_relation);
    virtual ~ShellMidSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    Real constrained_distance_;
    Real particle_spacing_ref_;
    LevelSetShape *level_set_shape_;
};

/**
 * @class ShellNormalDirectionPrediction
 * @brief predict the normal direction of shell particles.
 */
class ShellNormalDirectionPrediction : public BaseDynamics<void>
{
    const Real convergence_criterion_;
    const Real consistency_criterion_;

    void predictNormalDirection();
    void correctNormalDirection();

  public:
    explicit ShellNormalDirectionPrediction(BaseInnerRelation &inner_relation,
                                            Real thickness, Real consistency_criterion = cos(Pi / 20.0));
    virtual ~ShellNormalDirectionPrediction(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class NormalPrediction : public RelaxDataDelegateSimple, public LocalDynamics
    {
        Real thickness_;
        LevelSetShape *level_set_shape_;
        StdLargeVec<Vecd> &pos_, &n_, n_temp_;

      public:
        NormalPrediction(SPHBody &sph_body, Real thickness);
        virtual ~NormalPrediction(){};
        void update(size_t index_i, Real dt = 0.0);
    };

    class PredictionConvergenceCheck : public LocalDynamicsReduce<bool, ReduceAND>,
                                       public RelaxDataDelegateSimple
    {
      protected:
        const Real convergence_criterion_;
        StdLargeVec<Vecd> &n_, &n_temp_;

      public:
        PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion);
        virtual ~PredictionConvergenceCheck(){};

        bool reduce(size_t index_i, Real dt = 0.0);
    };

    class ConsistencyCorrection : public LocalDynamics, public RelaxDataDelegateInner
    {
      public:
        explicit ConsistencyCorrection(BaseInnerRelation &inner_relation, Real consistency_criterion);
        virtual ~ConsistencyCorrection(){};

        inline void interaction(size_t index_i, Real dt = 0.0)
        {
            mutex_modify_neighbor_.lock();
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                if (updated_indicator_[index_i] == 1)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (updated_indicator_[index_j] == 0)
                    {
                        updated_indicator_[index_j] = 1;
                        if (n_[index_i].dot(n_[index_j]) < consistency_criterion_)
                        {
                            if (n_[index_i].dot(-n_[index_j]) < consistency_criterion_)
                            {
                                n_[index_j] = n_[index_i];
                                updated_indicator_[index_j] = 2;
                            }
                            else
                            {
                                n_[index_j] = -n_[index_j];
                                updated_indicator_[index_j] = 1;
                            }
                        }
                    }
                }
            }
            mutex_modify_neighbor_.unlock();
        };

      protected:
        std::mutex mutex_modify_neighbor_; /**< mutex exclusion for memory conflict */
        const Real consistency_criterion_;
        StdLargeVec<int> updated_indicator_; /**> 0 not updated, 1 updated with reliable prediction, 2 updated from a reliable neighbor */
        StdLargeVec<Vecd> &n_;
    };

    class ConsistencyUpdatedCheck : public LocalDynamicsReduce<bool, ReduceAND>,
                                    public RelaxDataDelegateSimple
    {
      protected:
        StdLargeVec<int> &updated_indicator_;

      public:
        explicit ConsistencyUpdatedCheck(SPHBody &sph_body);
        virtual ~ConsistencyUpdatedCheck(){};

        bool reduce(size_t index_i, Real dt = 0.0);
    };

    class SmoothingNormal : public ParticleSmoothing<Vecd>
    {
      public:
        explicit SmoothingNormal(BaseInnerRelation &inner_relation);
        virtual ~SmoothingNormal(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
    };

    SimpleDynamics<NormalPrediction> normal_prediction_;
    ReduceDynamics<PredictionConvergenceCheck> normal_prediction_convergence_check_;
    InteractionDynamics<ConsistencyCorrection, execution::SequencedPolicy> consistency_correction_;
    ReduceDynamics<ConsistencyUpdatedCheck> consistency_updated_check_;
    InteractionWithUpdate<SmoothingNormal> smoothing_normal_;
};

/**
 * @class ShellRelaxationStepInner
 * @brief carry out particle relaxation step of particles within the shell body
 */
class ShellRelaxationStepInner : public RelaxationStepInner
{
  public:
    explicit ShellRelaxationStepInner(BaseInnerRelation &inner_relation, bool level_set_correction = false);
    virtual ~ShellRelaxationStepInner(){};

    SimpleDynamics<UpdateParticlePosition> update_shell_particle_position_;
    SimpleDynamics<ShellMidSurfaceBounding> mid_surface_bounding_;

    virtual void exec(Real dt = 0.0) override;
};
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_DYNAMICS_H
