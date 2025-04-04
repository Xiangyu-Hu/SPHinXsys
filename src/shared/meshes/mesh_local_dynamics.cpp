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

#include "mesh_local_dynamics.h"

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
void InitializeDataInACell::update(const Arrayi &cell_index)
{
    Vecd cell_position = mesh_data_.CellPositionFromIndex(cell_index);
    Real signed_distance = shape_.findSignedDistance(cell_position);
    Vecd normal_direction = shape_.findNormalDirection(cell_position);
    Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
    if (measure < grid_spacing_)
    {
        std::pair<size_t, int> occupied;
        occupied.first = SortIndexFromCellIndex(cell_index);
        occupied.second = 1;

        mesh_data_.assignDataPackageIndex(cell_index, 2);
        mesh_data_.registerOccupied(occupied);
    }
    else
    {
        size_t package_index = shape_.checkContain(cell_position) ? 0 : 1;
        mesh_data_.assignDataPackageIndex(cell_index, package_index);
    }
}
//=================================================================================================//
void TagACellIsInnerPackage::update(const Arrayi &cell_index)
{
    if (isInnerPackage(cell_index))
    {
        if (!mesh_data_.isInnerDataPackage(cell_index))
        {
            std::pair<size_t, int> occupied;
            occupied.first = SortIndexFromCellIndex(cell_index);
            occupied.second = 0;

            mesh_data_.registerOccupied(occupied);
        }
    }
}
//=================================================================================================//
void InitializeIndexMesh::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_.occupied_data_pkgs_[package_index-2].first;
    Arrayi cell_index = CellIndexFromSortIndex(sort_index);
    mesh_data_.assignDataPackageIndex(cell_index, package_index);
}
//=================================================================================================//
void UpdateLevelSetGradient::update(const size_t &index)
{
    mesh_data_.computeGradient(phi_, phi_gradient_, index);
}
//=================================================================================================//
void UpdateKernelIntegrals::update(const size_t &package_index)
{
    Arrayi cell_index = mesh_data_.meta_data_cell_[package_index].first;
    mesh_data_.assignByPosition(
        kernel_weight_, cell_index, [&](const Vecd &position) -> Real
        { return computeKernelIntegral(position); });
    mesh_data_.assignByPosition(
        kernel_gradient_, cell_index, [&](const Vecd &position) -> Vecd
        { return computeKernelGradientIntegral(position); });
    mesh_data_.assignByPosition(
        kernel_second_gradient_, cell_index, [&](const Vecd& position) -> Matd
        { return computeKernelSecondGradientIntegral(position); });
}
//=================================================================================================//
void InitializeDataInACellFromCoarse::update(const Arrayi &cell_index)
{
    Vecd cell_position = mesh_data_.CellPositionFromIndex(cell_index);
    MeshVariable<Real> &coarse_phi_ = *coarse_mesh_.getMeshVariable<Real>("Levelset");
    size_t package_index = coarse_mesh_.probeMesh(coarse_phi_, cell_position) < 0.0 ? 0 : 1;
    mesh_data_.assignDataPackageIndex(cell_index, package_index);
    if (coarse_mesh_.isWithinCorePackage(cell_position))
    {
        Real signed_distance = shape_.findSignedDistance(cell_position);
        Vecd normal_direction = shape_.findNormalDirection(cell_position);
        Real measure = (signed_distance * normal_direction).cwiseAbs().maxCoeff();
        if (measure < grid_spacing_)
        {
            std::pair<size_t, int> occupied;
            occupied.first = SortIndexFromCellIndex(cell_index);
            occupied.second = 1;

            mesh_data_.assignDataPackageIndex(cell_index, 2);
            mesh_data_.registerOccupied(occupied);
        }
    }
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
