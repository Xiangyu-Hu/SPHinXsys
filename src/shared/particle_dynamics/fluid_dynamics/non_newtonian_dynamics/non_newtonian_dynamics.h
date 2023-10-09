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
 * @file 	non_newtonian_dynamics.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Xiangyu Hu
 */

#ifndef NON_NEWTONIAN_DYNAMICS_H
#define NON_NEWTONIAN_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "fluid_integration.hpp"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class Oldroyd_BIntegration1stHalf;

using Integration1stHalfInnerDissipative =
    Integration1stHalf<Inner<>, DissipativeRiemannSolver, NoKernelCorrection>;

template <>
class Oldroyd_BIntegration1stHalf<Inner<>>
    : public Integration1stHalfInnerDissipative
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> tau_, dtau_dt_;
};

using Integration1stHalfWithWallDissipative =
    Integration1stHalf<ContactWall<>, DissipativeRiemannSolver, NoKernelCorrection>;
template <>
class Oldroyd_BIntegration1stHalf<ContactWall<>>
    : public Integration1stHalfWithWallDissipative
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &tau_;
};

template <typename... InteractionTypes>
class Oldroyd_BIntegration2ndHalf;

using Integration2ndHalfInnerDissipative =
    Integration2ndHalf<Inner<>, DissipativeRiemannSolver>;
template <>
class Oldroyd_BIntegration2ndHalf<Inner<>>
    : public Integration2ndHalfInnerDissipative
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};

using Integration2ndHalfWithWallDissipative =
    Integration2ndHalf<ContactWall<>, DissipativeRiemannSolver>;
template <>
class Oldroyd_BIntegration2ndHalf<ContactWall<>>
    : public Integration2ndHalfWithWallDissipative
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // NON_NEWTONIAN_DYNAMICS_H