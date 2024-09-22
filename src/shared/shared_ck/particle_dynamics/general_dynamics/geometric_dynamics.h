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
 * @file    geometric_dynamics.h
 * @brief   This is the particle dynamics applicable for all type bodies
 * @author	Xiangyu Hu
 */

#ifndef GEOMETRIC_DYNAMICS_H
#define GEOMETRIC_DYNAMICS_H

#include "base_general_dynamics.h"

namespace SPH
{
class NormalFromBodyShapeCK : public LocalDynamics
{
  public:
    explicit NormalFromBodyShapeCK(SPHBody &sph_body, Shape &shape);
    virtual ~NormalFromBodyShapeCK(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     NormalFromBodyShapeCK &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Shape *shape_;
        Vecd *pos_, *n_, *n0_;
        Real *phi_, *phi0_;
    };

  protected:
    Shape *shape_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_n_, *dv_n0_;
    DiscreteVariable<Real> *dv_phi_, *dv_phi0_;
};
} // namespace SPH
#endif // GEOMETRIC_DYNAMICS_H
