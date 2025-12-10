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
 * @file 	base_dynamics.h
 * @brief 	This is for the base classes of all dynamics.
 * @author	Xiangyu Hu
 */
#ifndef BASE_DYNAMICS_H
#define BASE_DYNAMICS_H

#include "base_data_type_package.h"

namespace SPH
{
class AbstractDynamics
{
  public:
    AbstractDynamics() {};
    virtual ~AbstractDynamics() {};
};

template <class ReturnType = void>
class BaseDynamics : public AbstractDynamics
{
  public:
    BaseDynamics() : AbstractDynamics(), is_newly_updated_(false) {};
    virtual ~BaseDynamics() {};
    bool checkNewlyUpdated() { return is_newly_updated_; };
    void setNotNewlyUpdated() { is_newly_updated_ = false; };

    template <class IdentifierType>
    void setUpdated(IdentifierType &identifier)
    {
        identifier.setNewlyUpdated();
        is_newly_updated_ = true;
    };

    /** There is the interface functions for computing. */
    virtual ReturnType exec(Real dt = 0.0) = 0;

  private:
    bool is_newly_updated_;
};
} // namespace SPH
#endif // BASE_DYNAMICS_H