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
 * @file 	io_plt.h
 * @brief 	Classes for save data in tecplot file format.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once

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
    PltEngine(){};
    virtual ~PltEngine(){};

    void writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name);
    void writeAQuantityHeader(std::ofstream &out_file, const Vecd &quantity, const std::string &quantity_name);
    void writeAQuantity(std::ofstream &out_file, const Real &quantity);
    void writeAQuantity(std::ofstream &out_file, const Vecd &quantity);
};

/**
 * @class BodyStatesRecordingToPlt
 * @brief  Write files for bodies
 * the output file is dat format can visualized by TecPlot
 */
class BodyStatesRecordingToPlt : public BodyStatesRecording
{
  public:
    BodyStatesRecordingToPlt(IOEnvironment &io_environment, SPHBody &body)
        : BodyStatesRecording(io_environment, body){};
    BodyStatesRecordingToPlt(IOEnvironment &io_environment, SPHBodyVector bodies)
        : BodyStatesRecording(io_environment, bodies){};
    virtual ~BodyStatesRecordingToPlt(){};

  protected:
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
    std::string filefullpath_;

  public:
    MeshRecordingToPlt(IOEnvironment &io_environment, BaseMeshField &mesh_field);
    virtual ~MeshRecordingToPlt(){};
    virtual void writeToFile(size_t iteration_step = 0) override;
};
} // namespace SPH
