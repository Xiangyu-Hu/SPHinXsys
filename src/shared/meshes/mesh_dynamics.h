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
#include "mesh_iterators.h"

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
        : mesh_data_(mesh_data),
          all_cells_(mesh_data.AllCells()),
          num_grid_pkgs_(mesh_data.num_grid_pkgs_),
          meta_data_cell_(mesh_data.meta_data_cell_){};
    virtual ~BaseMeshDynamics(){};

  protected:
    MeshWithGridDataPackages<4> &mesh_data_;
    Arrayi all_cells_;
    size_t &num_grid_pkgs_;
    std::pair<Arrayi, int>* &meta_data_cell_; 

    template <typename FunctionOnData>
    void grid_parallel_for(const FunctionOnData &function)
    {
        mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                          [&](Arrayi cell_index)
                          {
                              function(cell_index);
                          });
    }

    /** Iterator on a collection of mesh data packages. parallel computing. */
    template <typename FunctionOnData>
    void package_parallel_for(const FunctionOnData &function)
    {
        parallel_for(IndexRange(2, num_grid_pkgs_),
                    [&](const IndexRange &r)
                    {
                        for (size_t i = r.begin(); i != r.end(); ++i)
                        {
                            function(i);
                        }
                    },
                    ap);
    }
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
        grid_parallel_for(
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
        package_parallel_for(
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
        package_parallel_for(
          [&](size_t package_index)
          {
              std::pair<Arrayi, int> &metadata = meta_data_cell_[package_index];
              if (metadata.second == 1)
              {
                  this->update(package_index);
              }
          }
        );
    };
};
} // namespace SPH
#endif // MESH_DYNAMICS_H