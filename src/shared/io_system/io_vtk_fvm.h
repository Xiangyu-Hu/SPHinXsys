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
 * @file 	io_vtk_fvm.h
 * @brief Classes for input and output with vtk (Paraview) for FVM unstructured mesh.
 * @author Xiangyu Hu
 */

#ifndef IO_VTK_FVM_H
#define IO_VTK_FVM_H

#include "io_vtk.h"
#include "unstructured_mesh.h"

namespace SPH
{
/**
 * @class BodyStatesRecordingInMeshToVtp
 * @brief  Write files for bodies
 * the output file is VTK XML format in FVMcan visualized by ParaView the data type vtkPolyData
 */
class BodyStatesRecordingInMeshToVtp : public BodyStatesRecordingToVtp
{
  public:
    BodyStatesRecordingInMeshToVtp(SPHBody &body, ANSYSMesh &ansys_mesh);
    virtual ~BodyStatesRecordingInMeshToVtp(){};

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    StdLargeVec<Vecd> &node_coordinates_;
    StdLargeVec<StdVec<size_t>> &elements_nodes_connection_;
};

class BodyStatesRecordingInMeshToVtu : public BodyStatesRecordingToVtp
{
  public:
    BodyStatesRecordingInMeshToVtu(SPHBody &body, ANSYSMesh &ansys_mesh);
    virtual ~BodyStatesRecordingInMeshToVtu(){};

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    StdLargeVec<Vecd> &node_coordinates_;
    StdLargeVec<StdVec<size_t>> &elements_nodes_connection_;
    SPHBody &bounds_;
};
} // namespace SPH
#endif // IO_VTK_FVM_H
