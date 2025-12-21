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
 * @file 	io_plt.h
 * @brief 	Classes for save data in Tecplot file format.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef IO_PLT_H
#define IO_PLT_H

#include "io_base.h"

namespace SPH
{
/**
 * @class PltEngine
 * @brief The base class which defines Tecplot file related operation.
 */
class PltEngine
{
  public:
    PltEngine() {};
    virtual ~PltEngine() {};

    void writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name);
    template <int N>
    void writeAQuantityHeader(std::ofstream &out_file, const Eigen::Matrix<Real, N, 1> &quantity, const std::string &quantity_name);
    template <int N, int M>
    void writeAQuantityHeader(std::ofstream &out_file, const Eigen::Matrix<Real, N, M> &quantity, const std::string &quantity_name);
    void writeAQuantityHeader(std::ofstream &out_file, const SimTK::SpatialVec &quantity, const std::string &quantity_name);
    void writeAQuantity(std::ofstream &out_file, const Real &quantity);
    template <int N>
    void writeAQuantity(std::ofstream &out_file, const Eigen::Matrix<Real, N, 1> &quantity);
    template <int N, int M>
    void writeAQuantity(std::ofstream &out_file, const Eigen::Matrix<Real, N, M> &quantity);
    void writeAQuantity(std::ofstream &out_file, const SimTK::SpatialVec &quantity);
};

/**
 * @class BodyStatesRecordingToPlt
 * @brief  Write files for bodies
 * the output file is dat format can visualized by TecPlot
 */
class BodyStatesRecordingToPlt : public BodyStatesRecording
{
  public:
    BodyStatesRecordingToPlt(SPHBody &body) : BodyStatesRecording(body) {};
    BodyStatesRecordingToPlt(SPHSystem &sph_system) : BodyStatesRecording(sph_system) {};
    virtual ~BodyStatesRecordingToPlt() {};

  protected:
    void writePltFileHeader(std::ofstream &output_file, ParticleVariables &variables_to_write);
    void writePltFileParticleData(std::ofstream &output_file, ParticleVariables &variables_to_write, Vecd *position, size_t index);
    virtual void writeWithFileName(const std::string &sequence) override;
};

/**
 * @class MeshRecordingToPlt
 * @brief  write the mesh data in Tecplot format
 */
class MeshRecordingToPlt : public BaseIO
{
  protected:
    BaseMeshField &mesh_field_;
    std::string partial_file_name_;

  public:
    MeshRecordingToPlt(SPHSystem &sph_system, BaseMeshField &mesh_field);
    virtual ~MeshRecordingToPlt() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};
} // namespace SPH
#endif // IO_PLT_H
