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
 * @file 	mesh_dynamics.h
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_DYNAMICS_H
#define MESH_DYNAMICS_H

#include "mesh_local_dynamics.h"
#include "mesh_with_data_packages.hpp"

#include <functional>

using namespace std::placeholders;

namespace SPH
{
/**
 * @class BaseMeshDynamics
 * @brief The base class for all mesh dynamics
 * This class contains only the interface functions available
 * for all dynamics. An specific implementation should be realized.
 */
class BaseMeshDynamics
{
  public:
    BaseMeshDynamics(MeshWithGridDataPackages<4> &mesh_data)
        : mesh_data_(mesh_data){};
    virtual ~BaseMeshDynamics(){};

  protected:
    MeshWithGridDataPackages<4> &mesh_data_;
};

/**
 * @class MeshAllDynamics
 * @brief Mesh dynamics for all cell on the mesh
 */
template <class LocalDynamicsType>
class MeshAllDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
  public:
    template <typename... Args>
    MeshAllDynamics(MeshWithGridDataPackages<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data){};
    virtual ~MeshAllDynamics(){};

    void exec()
    {
        mesh_data_.grid_parallel_for(
            [&](Arrayi cell_index)
            {
              this->update(cell_index);
            }
        );
    };
};

/**
 * @class MeshInnerDynamics
 * @brief Mesh dynamics for only inner cells on the mesh
 */
template <class LocalDynamicsType>
class MeshInnerDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
  public:
    template <typename... Args>
    MeshInnerDynamics(MeshWithGridDataPackages<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data){};
    virtual ~MeshInnerDynamics(){};

    void exec()
    {
        mesh_data_.package_parallel_for(
            [&](size_t package_index)
            {
              this->update(package_index);
            }
        );
    };
};

/**
 * @class MeshCoreDynamics
 * @brief Mesh dynamics for only core cells on the mesh
 */
template <class LocalDynamicsType>
class MeshCoreDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    MeshCoreDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(identifier, std::forward<Args>(args)...),
          BaseMeshDynamics(identifier){};
    virtual ~MeshCoreDynamics(){};

    void exec()
    {
        mesh_data_.package_parallel_for(
          [&](size_t package_index)
          {
              std::pair<Arrayi, int> &metadata = mesh_data_.meta_data_cell_[package_index];
              if (metadata.second == 1)
              {
                  this->update(package_index);
              }
          }
        );
    };
};

/**
 * @class MeshSingleDynamics
 * @brief Mesh dynamics for only core cells on the mesh
 */
template <class LocalDynamicsType>
class MeshSingleDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    MeshSingleDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(identifier, std::forward<Args>(args)...),
          BaseMeshDynamics(identifier){};
    virtual ~MeshSingleDynamics(){};

    template <typename... Args>
    void exec(Args &&...args)
    {
        this->update(std::forward<Args>(args)...);
    };
};

template <typename ReturnType, class LocalDynamicsType>
class MeshCalculateDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    MeshCalculateDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(identifier, std::forward<Args>(args)...),
          BaseMeshDynamics(identifier){};
    virtual ~MeshCalculateDynamics(){};

    template <typename... Args>
    ReturnType exec(Args &&...args)
    {
        return this->update(std::forward<Args>(args)...);
    };
};
} // namespace SPH
#endif // MESH_DYNAMICS_H