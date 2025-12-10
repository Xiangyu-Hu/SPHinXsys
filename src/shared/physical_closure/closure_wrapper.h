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
 * @file 	closure_wrapper.h
 * @brief 	tbd.
 * @author	Xiangyu Hu
 */

#ifndef CLOSURE_WRAPPER_H
#define CLOSURE_WRAPPER_H

#include "base_data_type_package.h"
#include "sphinxsys_containers.h"

namespace SPH
{
class BaseParticles;

template <typename...>
class Closure;

template <>
class Closure<>
{
  public:
    Closure() {};
    virtual ~Closure() {};
    virtual void registerLocalParameters(BaseParticles *base_particles) {};
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) {};
    virtual void initializeLocalParameters(BaseParticles *base_particles) {};
};

template <class BaseModel, class... AuxiliaryModels>
class Closure<BaseModel, AuxiliaryModels...>
    : public BaseModel, public Closure<AuxiliaryModels...>
{
  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit Closure(FirstParameterSet &&first_parameter_set,
                     OtherParameterSets &&...other_parameter_sets)
        : BaseModel(first_parameter_set),
          Closure<AuxiliaryModels...>(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    virtual void registerLocalParameters(BaseParticles *base_particles) override
    {
        BaseModel::registerLocalParameters(base_particles);
        Closure<AuxiliaryModels...>::registerLocalParameters(base_particles);
    };

    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override
    {
        BaseModel::registerLocalParametersFromReload(base_particles);
        Closure<AuxiliaryModels...>::registerLocalParametersFromReload(base_particles);
    };

    virtual void initializeLocalParameters(BaseParticles *base_particles) override
    {
        BaseModel::initializeLocalParameters(base_particles);
        Closure<AuxiliaryModels...>::initializeLocalParameters(base_particles);
    };
};
} // namespace SPH
#endif // CLOSURE_WRAPPER_H