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
 * @file 	sphinxsys_containers.h
 * @brief 	Set up of basic data structure.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPHINXSYS_CONTAINERS_H
#define SPHINXSYS_CONTAINERS_H

#include "base_data_type_package.h"
#include "base_data_type.h"
#include "sphinxsys_constant.h"
#include "sphinxsys_variable.h"

namespace SPH
{
class Base;  // Indicating base class
class SingleValued
{
};
class Continuous
{
};
class Lattice;          // Indicating with lattice points
class UnstructuredMesh; // Indicating with unstructured mesh
class BaseMaterial;
class SPHBody;
class RealBody;
class SolidBody;
class BodyPart;
class BaseParticles;
class UserDefined; // Indicating with user defined type in apps
//----------------------------------------------------------------------
// Interaction type identifies
//----------------------------------------------------------------------
template <typename... InnerParameters>
class Inner; /**< Inner interaction: interaction within a body*/

template <typename... ContactParameters>
class Contact; /**< Contact interaction: interaction between a body with one or several another bodies */

class Boundary; /**< Interaction with boundary */
class Wall;     /**< Interaction with wall boundary */
class Extended; /**< An extened method of an interaction type */

template <typename...>
class Dirichlet; /**< Contact interaction with Dirichlet boundary condition */
//----------------------------------------------------------------------
// Time stepping type identifies
//----------------------------------------------------------------------
class ForwardEuler;
class RungeKutta;
class RungeKutta1stStage;
class RungeKutta2ndStage;
template <typename... ControlTypes>
class Dirichlet; /**< Contact interaction with Dirichlet boundary condition */
template <typename... ControlTypes>
class Neumann; /**< Contact interaction with Neumann boundary condition */
template <typename... ControlTypes>
class Robin; /**< Contact interaction with Neumann boundary condition */
//----------------------------------------------------------------------
// Spatial temporal type identifies
//----------------------------------------------------------------------
class SpatialTemporal;
//----------------------------------------------------------------------
// Other type identifies
//----------------------------------------------------------------------
using MaterialVector = StdVec<BaseMaterial *>;
using SPHBodyVector = StdVec<SPHBody *>;
using SolidBodyVector = StdVec<SolidBody *>;
using RealBodyVector = StdVec<RealBody *>;
using BodyPartVector = StdVec<BodyPart *>;

using IndexVector = StdVec<size_t>;
using ConcurrentIndexVector = ConcurrentVec<size_t>;
using ParticlesBound = std::pair<size_t, size_t>;

/** List data pair: first for indexes, second for particle position. */
using ListData = std::pair<size_t, Vecd>;
using ListDataVector = StdVec<ListData>;
using DataListsInCells = StdVec<ListDataVector *>;
using ConcurrentCellLists = ConcurrentVec<ConcurrentIndexVector *>;
/** Cell list for periodic boundary condition algorithms. */
using CellLists = std::pair<ConcurrentCellLists, DataListsInCells>;

/** Generalized particle data type */
typedef DataContainerAssemble<AllocatedData> ParticleData;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<DiscreteVariable> ParticleVariables;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<SingularVariable> SingularVariables;
} // namespace SPH
#endif // SPHINXSYS_CONTAINERS_H
