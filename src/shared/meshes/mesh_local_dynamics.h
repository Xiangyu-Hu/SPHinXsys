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
 * @file 	base_relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_LOCAL_DYNAMICS_H
#define MESH_LOCAL_DYNAMICS_H

#include "base_variable.h"
#include "mesh_with_data_packages.hpp"
#include "base_geometry.h"

namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
template <typename IndexType>
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackages<4> &mesh_data)
        : mesh_data_(mesh_data){};
    virtual ~BaseMeshLocalDynamics(){};

    virtual void update(const IndexType &index) = 0;
  protected:
    MeshWithGridDataPackages<4> &mesh_data_;
};

class InitializeDataInACell : public BaseMeshLocalDynamics<Arrayi>
{
  public:
    explicit InitializeDataInACell(MeshWithGridDataPackages<4> &mesh_data, Shape &shape)
        : BaseMeshLocalDynamics(mesh_data),
          all_cells_(mesh_data.AllCells()),
          shape_(shape){};
    virtual ~InitializeDataInACell(){};

    void update(const Arrayi &index);
  
  private:
    Shape &shape_;
    Real grid_spacing_ = mesh_data_.GridSpacing();
    Arrayi all_cells_;

    size_t SortIndexFromCellIndex(const Arrayi &cell_index);
};

class TagACellIsInnerPackage : public BaseMeshLocalDynamics<Arrayi>
{
  public:
    explicit TagACellIsInnerPackage(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data),
          all_cells_(mesh_data.AllCells()){};
    virtual ~TagACellIsInnerPackage(){};

    void update(const Arrayi &index);

  private:
    Arrayi all_cells_;

    bool isInnerPackage(const Arrayi &cell_index);
};

class InitializeIndexMesh : public BaseMeshLocalDynamics<size_t>
{
  public:
    explicit InitializeIndexMesh(MeshWithGridDataPackagesType &mesh_data)
        : BaseMeshLocalDynamics(mesh_data),
          all_cells_(mesh_data.AllCells()){};
    virtual ~InitializeIndexMesh(){};

    void update(const size_t &index);

  private:
    Arrayi all_cells_;
};

// class InitializeCellNeighborhood
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit InitializeCellNeighborhood(MeshWithGridDataPackages &mesh_data){};
//     virtual ~InitializeCellNeighborhood(){};

//     void update();
// };

class UpdateLevelSetGradient : public BaseMeshLocalDynamics<size_t>
{
  public:
    explicit UpdateLevelSetGradient(MeshWithGridDataPackages<4> &mesh_data)
        : BaseMeshLocalDynamics(mesh_data),
          phi_(*mesh_data.getMeshVariable<Real>("Levelset")),
          phi_gradient_(*mesh_data.getMeshVariable<Vecd>("LevelsetGradient")){};
    virtual ~UpdateLevelSetGradient(){};

    void update(const size_t &index);
    
  private:
    MeshVariable<Real> &phi_;
    MeshVariable<Vecd> &phi_gradient_;
};

// class UpdateKernelIntegrals
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit UpdateKernelIntegrals(MeshWithGridDataPackages &mesh_data){};
//     virtual ~UpdateKernelIntegrals(){};

//     void update(size_t package_index);

//   private:
//     Real computeKernelIntegral(const Vecd &position);
//     Vecd computeKernelGradientIntegral(const Vecd &position);
// };
// class ReinitializeLevelSet
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit ReinitializeLevelSet(MeshWithGridDataPackages &mesh_data){};
//     virtual ~ReinitializeLevelSet(){};

//     void update(size_t package_index);
// };
// class MarkNearInterface
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit MarkNearInterface(MeshWithGridDataPackages &mesh_data){};
//     virtual ~MarkNearInterface(){};

//     void update(size_t package_index);
// };
// class RedistanceInterface
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit RedistanceInterface(MeshWithGridDataPackages &mesh_data){};
//     virtual ~RedistanceInterface(){};

//     void update(size_t package_index);
// };
// class DiffuseLevelSetSign
//     : public BaseMeshLocalDynamics
// {
//   public:
//     explicit DiffuseLevelSetSign(MeshWithGridDataPackages &mesh_data){};
//     virtual ~DiffuseLevelSetSign(){};

//     void update(size_t package_index);
// };

// } // namespace mesh_dynamics
} // namespace SPH
#endif // MESH_LOCAL_DYNAMICS_H
