/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef RELAX_DYNAMICS_SURFACE_H
#define RELAX_DYNAMICS_SURFACE_H

#include "sphinxsys.h"
#include "vector.h"
#include "surface_shape.h"

namespace SPH
{
	
    class SurfaceShape;

	namespace relax_dynamics
    {
     class ShapeSurfaceBounding2 : public LocalDynamics,
                                  public RelaxDataDelegateSimple
    {
      public:
        ShapeSurfaceBounding2(RealBody &real_body_);
        virtual ~ShapeSurfaceBounding2(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        StdLargeVec<Vecd> &pos_;
        Shape *shape_;
    };


      class RelaxationStepInnerFirstHalf : public BaseDynamics<void>
    {
      public:
        explicit RelaxationStepInnerFirstHalf(BaseInnerRelation &inner_relation,
                                              bool level_set_correction = false);
        virtual ~RelaxationStepInnerFirstHalf(){};
        virtual void exec(Real dt = 0.0) override;
       

      protected:
        RealBody *real_body_;
        BaseInnerRelation &inner_relation_;
        UniquePtr<BaseDynamics<void>> relaxation_acceleration_inner_;
    };

    class RelaxationStepInnerSecondHalf : public BaseDynamics<void>
    {
      public:
        explicit RelaxationStepInnerSecondHalf(BaseInnerRelation &inner_relation,
                                               bool level_set_correction = false);
        virtual ~RelaxationStepInnerSecondHalf(){};
        SimpleDynamics<ShapeSurfaceBounding2> &SurfaceBounding() { return surface_bounding_; };
        virtual void exec(Real dt = 0.0) override;
        

      protected:
        RealBody *real_body_;
        ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
        SimpleDynamics<UpdateParticlePosition> update_particle_position_;
        SimpleDynamics<ShapeSurfaceBounding2> surface_bounding_;
    };

    /**
     * @class SurfaceNormalDirection
     * @brief get the normal direction of surface particles.
     */
    class SurfaceNormalDirection : public RelaxDataDelegateSimple, public LocalDynamics
    {
      public:
        explicit SurfaceNormalDirection(SPHBody &sph_body);
        virtual ~SurfaceNormalDirection(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        SurfaceShape *surface_shape_;
        StdLargeVec<Vecd> &pos_, &n_;
    };


    /**@class ConstrainSuefaceBodyRegion
     * @brief Fix the position surafce body part.
     */
    class ConstrainSurfaceBodyRegion : public BaseLocalDynamics<BodyPartByParticle>, public RelaxDataDelegateSimple
    {
      public:
        ConstrainSurfaceBodyRegion(BodyPartByParticle &body_part);
        virtual ~ConstrainSurfaceBodyRegion(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        StdLargeVec<Vecd> &acc_;
    };


    }
  }
#endif // RELAX_DYNAMICS_H
