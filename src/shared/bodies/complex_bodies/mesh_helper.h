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
 * @file mesh_helper.h
 * @brief Here, we define the functions used for reading the mesh data.
 * @author Yash Mandaokar Zhentong Wang and Xiangyu Hu
 */
#ifndef MESH_HELPER_H
#define MESH_HELPER_H

#include "unstructured_mesh.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace SPH
{

class MeshFileHelpers
{
  public:
    MeshFileHelpers() {};
    virtual ~MeshFileHelpers() {};

    static void meshDimension(std::ifstream &mesh_file, size_t &dimension, std::string &text_line);
    static void numberOfNodes(std::ifstream &mesh_file, size_t &number_of_points, std::string &text_line);
    static void nodeCoordinates(std::ifstream &mesh_file, StdVec<Vecd> &node_coordinates_, std::string &text_line, size_t &dimension);
    static void numberOfElements(std::ifstream &mesh_file, size_t &number_of_elements, std::string &text_line);
    static void dataStruct(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_, StdVec<StdVec<size_t>> &elements_nodes_connection_,
                           size_t number_of_elements, size_t mesh_type, size_t dimension);
    static size_t findBoundaryType(std::string &text_line, size_t boundary_type);
    static Vecd nodeIndex(std::string &text_line);
    static Vec2d cellIndex(std::string &text_line);
    static void updateElementsNodesConnection(StdVec<StdVec<size_t>> &elements_nodes_connection_, Vecd nodes, Vec2d cells,
                                              bool &check_neighbor_cell1, bool &check_neighbor_cell2);
    static void updateCellLists(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_, StdVec<StdVec<size_t>> &elements_nodes_connection_, Vecd nodes,
                                Vec2d cells, bool &check_neighbor_cell1, bool &check_neighbor_cell2, size_t boundary_type);
    static void updateBoundaryCellLists(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_, StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                        Vecd nodes, Vec2d cells, bool &check_neighbor_cell1, bool &check_neighbor_cell2, size_t boundary_type);
    static void cellCenterCoordinates(StdVec<StdVec<size_t>> &elements_nodes_connection_, std::size_t &element,
                                      StdVec<Vecd> &node_coordinates_, StdVec<Vecd> &elements_center_coordinates_, Vecd &center_coordinate);
    static void elementVolume(StdVec<StdVec<size_t>> &elements_nodes_connection_, std::size_t &element,
                              StdVec<Vecd> &node_coordinates_, StdVec<Real> &elements_volumes_);
    static void minimumDistance(StdVec<Real> &all_data_of_distance_between_nodes, StdVec<Real> &elements_volumes_,
                                StdVec<StdVec<StdVec<size_t>>> &mesh_topology_, StdVec<Vecd> &node_coordinates_);
    /*Functions for .msh file from FLUENT*/
    static void numberOfNodesFluent(std::ifstream &mesh_file, size_t &number_of_points, std::string &text_line);
    static void numberOfElementsFluent(std::ifstream &mesh_file, size_t &number_of_elements, std::string &text_line);
    static void nodeCoordinatesFluent(std::ifstream &mesh_file, StdVec<Vecd> &node_coordinates_, std::string &text_line,
                                      size_t &dimension);
    static void updateBoundaryCellListsFluent(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_, StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                              Vecd nodes, Vec2d cells, bool &check_neighbor_cell1, bool &check_neighbor_cell2, size_t boundary_type);
};

} // namespace SPH
#endif // MESH_HELPER_H