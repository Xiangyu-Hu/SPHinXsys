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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	fluid_dynamics_complex_correction.h
 * @brief Here, we define the algorithm classes for fluid dynamics,
 *        in which correction matrix is used to increase the approximation
 *        of pressure gradient.
 * @author Yaru Ren and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_COMPLEX_CORRECTION_H
#define FLUID_DYNAMICS_COMPLEX_CORRECTION_H

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_inner_correction.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class BaseIntegration1stHalfCorrectWithWall
 * @brief  template class pressure relaxation scheme together with wall boundary
 */
template <class BaseIntegration1stHalfCorrectType>
class BaseIntegration1stHalfCorrectWithWall : public InteractionWithWall<BaseIntegration1stHalfCorrectType>
{
  public:
    template <typename... Args>
    BaseIntegration1stHalfCorrectWithWall(Args &&...args)
        : InteractionWithWall<BaseIntegration1stHalfCorrectType>(std::forward<Args>(args)...){};
    virtual ~BaseIntegration1stHalfCorrectWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

using Integration1stHalfCorrectWithWall = BaseIntegration1stHalfCorrectWithWall<Integration1stHalfCorrect>;
using Integration1stHalfRiemannCorrectWithWall = BaseIntegration1stHalfCorrectWithWall<Integration1stHalfRiemannCorrect>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_COMPLEX_CORRECTION_H