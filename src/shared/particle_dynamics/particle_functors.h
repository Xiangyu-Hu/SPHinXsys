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
 * @file    particle_functors.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_FUNCTORS_H
#define PARTICLE_FUNCTORS_H

#include "base_particles.hpp"
#include "reduce_functors.h"

namespace SPH
{
/**
 * @class WithinScope
 * Base class introduce the concept of "within the scope".
 * The derived class should implement the operator bool(size_t...)
 * to indicate whether an indexed element is within the scope.
 * Generally, the object of the derived class
 * should be named as "within_scope" or "within_scope_" (class member)
 * so that the code can be more readable.
 */
class WithinScope
{
};

//----------------------------------------------------------------------
// Particle scope functors
//---------------------------------------------------------------------

class AllParticles : public WithinScope
{
  public:
    explicit AllParticles(BaseParticles *base_particles) : WithinScope() {};
    bool operator()(size_t index_i)
    {
        return true;
    };
};

template <int INDICATOR>
class IndicatedParticles : public WithinScope
{
    int *indicator_;

  public:
    explicit IndicatedParticles(BaseParticles *base_particles)
        : WithinScope(),
          indicator_(base_particles->getVariableDataByName<int>("Indicator")) {};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] == INDICATOR;
    };
};

using BulkParticles = IndicatedParticles<0>;

template <int INDICATOR>
class NotIndicatedParticles : public WithinScope
{
    int *indicator_;

  public:
    explicit NotIndicatedParticles(BaseParticles *base_particles)
        : WithinScope(),
          indicator_(base_particles->getVariableDataByName<int>("Indicator")) {};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] != INDICATOR;
    };
};

//----------------------------------------------------------------------
// Particle parameter functors
//----------------------------------------------------------------------
class ParticleParameter // base class to indicate the concept of particle parameter
{
};

template <typename T>
class ParameterFixed : public ParticleParameter
{
    T parameter_;

  public:
    explicit ParameterFixed(const T &c) : ParticleParameter(), parameter_(c) {};
    T &operator()(size_t index_i)
    {
        return parameter_;
    };
};

template <typename T>
class ParameterVariable : public ParticleParameter
{
    T *v_;

  public:
    explicit ParameterVariable(T *v) : ParticleParameter(), v_(v) {};
    T &operator()(size_t index_i)
    {
        return v_[index_i];
    };
};
//----------------------------------------------------------------------
// Particle average functors
//----------------------------------------------------------------------
class ParticleAverage // base class to indicate the concept of particle average
{
};

template <typename T>
class PairAverageFixed : public ParticleAverage
{
    T average_;

  public:
    PairAverageFixed(const T &c1, const T &c2)
        : ParticleAverage(), average_(0.5 * (c1 + c2)) {};
    explicit PairAverageFixed(const T &c)
        : PairAverageFixed(c, c) {};
    T &operator()(size_t index_i, size_t index_j)
    {
        return average_;
    };
};

class GeomAverage : public ParticleAverage
{
  public:
    GeomAverage() : ParticleAverage() {};

  protected:
    Real inverse(const Real &x) { return 1.0 / x; };

    template <typename MatrixType>
    MatrixType inverse(const MatrixType &x) { return x.inverse(); };
};

template <typename T>
class PairGeomAverageFixed : public GeomAverage
{
    T geom_average_;

  public:
    PairGeomAverageFixed(const T &c1, const T &c2)
        : GeomAverage(), geom_average_(2.0 * c1 * c2 * inverse(c1 + c2)) {};
    explicit PairGeomAverageFixed(const T &c)
        : PairGeomAverageFixed(c, c) {};
    T &operator()(size_t index_i, size_t index_j)
    {
        return geom_average_;
    };
};

template <typename T>
class PairAverageVariable : public ParticleAverage
{
    T *v1_, *v2_;

  public:
    PairAverageVariable(T *v1, T *v2)
        : ParticleAverage(), v1_(v1), v2_(v2) {};
    explicit PairAverageVariable(T *v)
        : PairAverageVariable(v, v) {};
    T operator()(size_t index_i, size_t index_j)
    {
        return 0.5 * (v1_[index_i] + v2_[index_j]);
    };
};

template <typename T>
class PairGeomAverageVariable : public GeomAverage
{
    T *v1_, *v2_;

  public:
    PairGeomAverageVariable(T *v1, T *v2)
        : GeomAverage(), v1_(v1), v2_(v2) {};
    explicit PairGeomAverageVariable(T *v)
        : PairGeomAverageVariable(v, v) {};

    T operator()(size_t index_i, size_t index_j)
    {
        return 2.0 * v1_[index_i] * v2_[index_j] * inverse(v1_[index_i] + v2_[index_j]);
    };
};
//----------------------------------------------------------------------
// Particle kernel correction functors
//----------------------------------------------------------------------
class KernelCorrection // base class to indicate the concept of kernel correction
{
};

class NoKernelCorrection : public KernelCorrection
{
  public:
    NoKernelCorrection(BaseParticles *particles) : KernelCorrection() {};
    Real operator()(size_t index_i, size_t index_j)
    {
        return 1.0;
    };

    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};

class LinearGradientCorrection : public KernelCorrection
{
  public:
    LinearGradientCorrection(BaseParticles *particles)
        : KernelCorrection(),
          B_(particles->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {};

    Matd operator()(size_t index_i, size_t index_j)
    {
        return B_[index_i];
    };

    Matd operator()(size_t index_i)
    {
        return B_[index_i];
    };

  protected:
    Matd *B_;
};

class LinearGradientCorrectionWithBulkScope : public KernelCorrection
{
  public:
    LinearGradientCorrectionWithBulkScope(BaseParticles *particles)
        : KernelCorrection(),
          B_(particles->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          within_scope_(particles) {};

    Matd operator()(size_t index_j, size_t index_i)
    {
        return within_scope_(index_j) ? B_[index_j] : B_[index_i];
    };

    Matd operator()(size_t index_i)
    {
        return B_[index_i];
    };

  protected:
    Matd *B_;
    BulkParticles within_scope_;
};

class SingleResolution
{
  public:
    SingleResolution(BaseParticles *particles) {};
    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};
//----------------------------------------------------------------------
// Particle adaptation functors
//----------------------------------------------------------------------
class AdaptiveResolution
{
  public:
    AdaptiveResolution(BaseParticles *particles)
        : h_ratio_(particles->getVariableDataByName<Real>("SmoothingLengthRatio")) {};

    Real operator()(size_t index_i)
    {
        return h_ratio_[index_i];
    };

  protected:
    Real *h_ratio_;
};
} // namespace SPH
#endif // PARTICLE_FUNCTORS_H
