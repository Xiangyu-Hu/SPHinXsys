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
 * @file general_constraint_ck.h
 * @brief Particles are constrained on their position according to
 * different rules.
 * @author	Xiangyu Hu
 */

#ifndef GENERAL_CONSTRAINT_CK_H
#define GENERAL_CONSTRAINT_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class DynamicsIdentifier, typename DataType>
class ConstantConstraintCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    ConstantConstraintCK(DynamicsIdentifier &identifier, const std::string &variable_name,
                         const DataType &constrained_value)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
          constrained_value_(constrained_value){};
    virtual ~ConstantConstraintCK(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
              constrained_value_(encloser.constrained_value_){};

        void update(size_t index_i, Real dt = 0.0)
        {
            variable_[index_i] = constrained_value_;
        };

      protected:
        DataType *variable_;
        DataType constrained_value_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
    DataType constrained_value_;
};
} // namespace SPH
#endif // GENERAL_CONSTRAINT_CK_H
