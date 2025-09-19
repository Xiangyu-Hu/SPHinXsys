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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    geometric_dynamics.h
 * @brief   This is the particle dynamics applicable for all type bodies
 * @author	Xiangyu Hu
 */

#ifndef GEOMETRIC_DYNAMICS_H
#define GEOMETRIC_DYNAMICS_H

#include "base_general_dynamics.h"

namespace SPH
{
class HostKernel
{
  public:
    template <class ExecutionPolicy, class EncloserType>
    HostKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    {
        // not implemented for device policy due to virtual function call in inital_shape_,
        // which is not allowed in device code
        static_assert(!std::is_base_of<execution::DeviceExecution<>, ExecutionPolicy>::value,
                      "This compute kernel is not designed for execution on device!");
    }
};

class NormalFromBodyShapeCK : public LocalDynamics
{
  public:
    explicit NormalFromBodyShapeCK(SPHBody &sph_body);
    virtual ~NormalFromBodyShapeCK() {};

    class UpdateKernel : public HostKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     NormalFromBodyShapeCK &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Shape *initial_shape_;
        Vecd *pos_, *n_, *n0_;
        Real *phi_, *phi0_;
    };

  protected:
    Shape *initial_shape_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_n_, *dv_n0_;
    DiscreteVariable<Real> *dv_phi_, *dv_phi0_;
};

class SurfaceIndicationFromBodyShape : public LocalDynamics
{
  public:
    explicit SurfaceIndicationFromBodyShape(SPHBody &sph_body);
    virtual ~SurfaceIndicationFromBodyShape() {};

    class UpdateKernel : public HostKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     SurfaceIndicationFromBodyShape &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Shape *initial_shape_;
        Real spacing_ref_;
        int *indicator_;
        Vecd *pos_;
    };

  protected:
    Shape *initial_shape_;
    Real spacing_ref_;
    DiscreteVariable<int> *dv_indicator_;
    DiscreteVariable<Vecd> *dv_pos_;
};

class RandomizeParticlePositionCK : public LocalDynamics
{
  public:
    explicit RandomizeParticlePositionCK(SPHBody &sph_body, Real randomize_factor = 0.25);
    virtual ~RandomizeParticlePositionCK() {};

    class UpdateKernel : public HostKernel // due to random number generation
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *pos_;
        Real randomize_scale_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    Real randomize_scale_;
};
} // namespace SPH
#endif // GEOMETRIC_DYNAMICS_H
