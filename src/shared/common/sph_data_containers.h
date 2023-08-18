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
 * @file 	sph_data_containers.h
 * @brief 	Set up of basic data structure.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPH_DATA_CONTAINERS_H
#define SPH_DATA_CONTAINERS_H

#include "base_data_package.h"
#include "base_data_type.h"
#include "base_variable.h"

namespace SPH
{
/**
 * @brief Pre-claimed classes.
 */
class BaseMaterial;
class SPHBody;
class RealBody;
class SolidBody;
class BodyPart;
class BaseParticles;

using MaterialVector = StdVec<BaseMaterial *>;
using SPHBodyVector = StdVec<SPHBody *>;
using SolidBodyVector = StdVec<SolidBody *>;
using RealBodyVector = StdVec<RealBody *>;
using BodyPartVector = StdVec<BodyPart *>;

using IndexVector = StdVec<size_t>;
using ConcurrentIndexVector = ConcurrentVec<size_t>;

/** List data pair: first for indexes, second for particle position. */
using ListData = std::tuple<size_t, Vecd, Real>;
using ListDataVector = StdLargeVec<ListData>;
using DataListsInCells = StdLargeVec<ListDataVector *>;
using ConcurrentCellLists = ConcurrentVec<ConcurrentIndexVector *>;
/** Cell list for splitting algorithms. */
using SplitCellLists = StdVec<ConcurrentCellLists>;
/** Cell list for periodic boundary condition algorithms. */
using CellLists = std::pair<ConcurrentCellLists, DataListsInCells>;

/** Generalized particle data type */
typedef DataContainerAddressAssemble<StdLargeVec> ParticleData;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<DiscreteVariable> ParticleVariables;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<GlobalVariable> GlobalVariables;

/** Generalized mesh data type */
template <typename DataType>
using MeshVariable = DiscreteVariable<DataType>;
typedef DataContainerAddressAssemble<MeshVariable> MeshVariableAssemble;

/** operation by looping or going through a particle variables */
template <typename DataType>
struct loopParticleVariables
{
    template <typename VariableOperation>
    void operator()(ParticleData &particle_data,
                    ParticleVariables &particle_variables, VariableOperation &variable_operation) const
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        for (DiscreteVariable<DataType> *variable : std::get<type_index>(particle_variables))
        {
            StdLargeVec<DataType> &variable_data = *(std::get<type_index>(particle_data)[variable->IndexInContainer()]);
            variable_operation(variable->Name(), variable_data);
        }
    };
};
} // namespace SPH
#endif // SPH_DATA_CONTAINERS_H
