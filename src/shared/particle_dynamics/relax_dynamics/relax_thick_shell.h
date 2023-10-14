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
 * @file 	relax_bounding.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Xiangyu Hu
 */

#ifndef RELAX_BOUNDING_H
#define RELAX_BOUNDING_H

#include "base_relax_dynamics.h"

namespace SPH
{
class GeometryShape;
class LevelSetShape;

namespace relax_dynamics
{
/**
 * @class ShellMidSurfaceBounding
 * @brief constrain particles by constraining particles to mid-surface.
 * Note that level_set_refinement_ratio should be smaller than particle_spacing_ref_ / (0.05 * thickness_)
 * because if level_set_refinement_ratio > particle_spacing_ref_ / (0.05 * thickness_),
 * there will be no level set field.
 */
class ShellMidSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                                public RelaxDataDelegateSimple
{
  public:
    ShellMidSurfaceBounding(NearShapeSurface &body_part, Real thickness, Real level_set_refinement_ratio);
    virtual ~ShellMidSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    Real constrained_distance_;
    Real particle_spacing_ref_, thickness_, level_set_refinement_ratio_;
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

        void interaction(size_t index_i, Real dt = 0.0);

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
    explicit ShellRelaxationStepInner(BaseInnerRelation &inner_relation, Real thickness,
                                      Real level_set_refinement_ratio, bool level_set_correction = false);
    virtual ~ShellRelaxationStepInner(){};

    SimpleDynamics<UpdateParticlePosition> update_shell_particle_position_;
    SimpleDynamics<ShellMidSurfaceBounding> mid_surface_bounding_;

    virtual void exec(Real dt = 0.0) override;
};
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_BOUNDING_H
