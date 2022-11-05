/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	combined_particle_dynamics.h
 * @brief 	This is the classes for combined particle dynamics
 *           which share the particle loop but are independent from each other,
 *          aiming to increase computing intensity under the data caching environment
 * @author	Xiangyu Hu
 */

#ifndef COMBINED_PARTICLE_DYNAMICS_H
#define COMBINED_PARTICLE_DYNAMICS_H

#include "particle_dynamics_algorithms.h"

namespace SPH
{
    /**
     * @class CombinedLocalDynamicsSimple
     * @brief
     */
    template <typename... MultipleLocalDynamics>
    class CombinedLocalDynamicsSimple;

    template <>
    class CombinedLocalDynamicsSimple<> : public LocalDynamics
    {
    public:
        template <class BodyRelationType>
        CombinedLocalDynamicsSimple(BodyRelationType &body_relation)
            : LocalDynamics(body_relation.getDynamicsRange()){};

        void update(size_t index_i, Real dt = 0.0){};
    };

    template <class FirstLocalDynamics, class... OtherLocalDynamics>
    class CombinedLocalDynamicsSimple<FirstLocalDynamics, OtherLocalDynamics...> : public LocalDynamics
    {
    protected:
        FirstLocalDynamics first_local_dynamics_;
        CombinedLocalDynamicsSimple<OtherLocalDynamics...> other_local_dynamics_;

    public:
        template <typename BodyRelationType, typename... FirstArgs, typename... OtherArgs>
        CombinedLocalDynamicsSimple(BodyRelationType &body_relation, FirstArgs &&...first_args, OtherArgs &&...other_args)
            : LocalDynamics(body_relation.getDynamicsRange()),
              first_local_dynamics_(body_relation, std::forward<FirstArgs>(first_args)...),
              other_local_dynamics_(body_relation, std::forward<OtherArgs>(other_args)...){};

        virtual void setupDynamics(Real dt = 0.0) override
        {
            first_local_dynamics_.setupDynamics(dt);
            other_local_dynamics_.setupDynamics(dt);
        };

        void update(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.update(index_i, dt);
            other_local_dynamics_.update(index_i, dt);
        };
    };

    /**
     * @class CombinedLocalInteraction
     * @brief
     */
    template <typename... MultipleLocalDynamics>
    class CombinedLocalInteraction;

    template <>
    class CombinedLocalInteraction<> : public LocalDynamics
    {
    public:
        template <class BodyRelationType>
        CombinedLocalInteraction(BodyRelationType &body_relation)
            : LocalDynamics(body_relation.getDynamicsRange()){};

        void interaction(size_t index_i, Real dt = 0.0){};
    };

    template <class FirstLocalDynamics, class... OtherLocalDynamics>
    class CombinedLocalInteraction<FirstLocalDynamics, OtherLocalDynamics...> : public LocalDynamics
    {
    protected:
        FirstLocalDynamics first_local_dynamics_;
        CombinedLocalInteraction<OtherLocalDynamics...> other_local_dynamics_;

    public:
        template <typename BodyRelationType, typename... FirstArgs, typename... OtherArgs>
        CombinedLocalInteraction(BodyRelationType &body_relation, FirstArgs &&...first_args, OtherArgs &&...other_args)
            : LocalDynamics(body_relation.getDynamicsRange()),
              first_local_dynamics_(body_relation, std::forward<FirstArgs>(first_args)...),
              other_local_dynamics_(body_relation, std::forward<OtherArgs>(other_args)...){};

        virtual void setupDynamics(Real dt = 0.0) override
        {
            first_local_dynamics_.setupDynamics(dt);
            other_local_dynamics_.setupDynamics(dt);
        };

        void interaction(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.interaction(index_i, dt);
            other_local_dynamics_.interaction(index_i, dt);
        };
    };

    /**
     * @class CombinedLocalInteractionWithUpdate
     * @brief
     */
    template <typename... MultipleLocalDynamics>
    class CombinedLocalInteractionWithUpdate;

    template <>
    class CombinedLocalInteractionWithUpdate<> : public LocalDynamics
    {
    public:
        template <class BodyRelationType>
        CombinedLocalInteractionWithUpdate(BodyRelationType &body_relation)
            : LocalDynamics(body_relation.getDynamicsRange()){};

        void interaction(size_t index_i, Real dt = 0.0){};
        void update(size_t index_i, Real dt = 0.0){};
    };

    template <class FirstLocalDynamics, class... OtherLocalDynamics>
    class CombinedLocalInteractionWithUpdate<FirstLocalDynamics, OtherLocalDynamics...> : public LocalDynamics
    {
    protected:
        FirstLocalDynamics first_local_dynamics_;
        CombinedLocalInteractionWithUpdate<OtherLocalDynamics...> other_local_dynamics_;

    public:
        template <typename BodyRelationType, typename... FirstArgs, typename... OtherArgs>
        CombinedLocalInteractionWithUpdate(BodyRelationType &body_relation, FirstArgs &&...first_args, OtherArgs &&...other_args)
            : LocalDynamics(body_relation.getDynamicsRange()),
              first_local_dynamics_(body_relation, std::forward<FirstArgs>(first_args)...),
              other_local_dynamics_(body_relation, std::forward<OtherArgs>(other_args)...){};

        virtual void setupDynamics(Real dt = 0.0) override
        {
            first_local_dynamics_.setupDynamics(dt);
            other_local_dynamics_.setupDynamics(dt);
        };

        void interaction(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.interaction(index_i, dt);
            other_local_dynamics_.interaction(index_i, dt);
        };

        void update(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.update(index_i, dt);
            other_local_dynamics_.update(index_i, dt);
        };
    };

    /**
     * @class CombinedLocalDynamics1Level
     * @brief
     */
    template <typename... MultipleLocalDynamics>
    class CombinedLocalDynamics1Level;

    template <>
    class CombinedLocalDynamics1Level<> : public LocalDynamics
    {
    public:
        template <class BodyRelationType>
        CombinedLocalDynamics1Level(BodyRelationType &body_relation)
            : LocalDynamics(body_relation.getDynamicsRange()){};

        void initialize(size_t index_i, Real dt = 0.0){};
        void interaction(size_t index_i, Real dt = 0.0){};
        void update(size_t index_i, Real dt = 0.0){};
    };

    template <class FirstLocalDynamics, class... OtherLocalDynamics>
    class CombinedLocalDynamics1Level<FirstLocalDynamics, OtherLocalDynamics...> : public LocalDynamics
    {
    protected:
        FirstLocalDynamics first_local_dynamics_;
        CombinedLocalDynamics1Level<OtherLocalDynamics...> other_local_dynamics_;

    public:
        template <typename BodyRelationType, typename... FirstArgs, typename... OtherArgs>
        CombinedLocalDynamics1Level(BodyRelationType &body_relation, FirstArgs &&...first_args, OtherArgs &&...other_args)
            : LocalDynamics(body_relation.getDynamicsRange()),
              first_local_dynamics_(body_relation, std::forward<FirstArgs>(first_args)...),
              other_local_dynamics_(body_relation, std::forward<OtherArgs>(other_args)...){};

        virtual void setupDynamics(Real dt = 0.0) override
        {
            first_local_dynamics_.setupDynamics(dt);
            other_local_dynamics_.setupDynamics(dt);
        };

        void initialize(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.initialize(index_i, dt);
            other_local_dynamics_.initialize(index_i, dt);
        };

        void interaction(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.interaction(index_i, dt);
            other_local_dynamics_.interaction(index_i, dt);
        };

        void update(size_t index_i, Real dt = 0.0)
        {
            first_local_dynamics_.update(index_i, dt);
            other_local_dynamics_.update(index_i, dt);
        };
    };
}
#endif // COMBINED_PARTICLE_DYNAMICS_H