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
 * @file 	fluid_mixture_dynamics.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_MIXTURE_DYNAMICS_H
#define FLUID_MIXTURE_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class LocalDynamicsType>
class MultiPhaseAcousticStep : public AcousticStep<LocalDynamicsType>
{
  public:
    template <class DynamicsIdentifier>
    explicit MultiPhaseAcousticStep(DynamicsIdentifier &identifier);
    virtual ~MultiPhaseAcousticStep(){};

  protected:
    WeaklyCompressibleMultiPhase &multi_phase_fluid_;
    DiscreteVariable<Real> *dv_phi_list_, *dv_mass_list_, *dev_dmass_dt_list_;
    DiscreteVariable<Vecd> *dv_velocity_list_, dv_momentum_list_, *dv_dmomentum_dt_list_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_MIXTURE_DYNAMICS_H
