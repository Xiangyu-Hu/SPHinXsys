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
 * @file 	io_vtk_mesh_ck.h
 * @brief Classes for input and output with vtk (Paraview) for FVM unstructured mesh.
 * @author Xiangyu Hu
 */

#ifndef IO_VTK_MESH_CK_H
#define IO_VTK_MESH_CK_H

#include "execution_policy.h"
#include "io_vtk_mesh.h"

namespace SPH
{
template <class ExecutionPolicy>
class BodyStatesRecordingToTriangleMeshVtpCK : public BodyStatesRecordingToTriangleMeshVtp
{
  public:
    template <typename... Args>
    BodyStatesRecordingToTriangleMeshVtpCK(Args &&...args)
        : BodyStatesRecordingToTriangleMeshVtp(std::forward<Args>(args)...){};
    virtual ~BodyStatesRecordingToTriangleMeshVtpCK() {};

    virtual void writeToFile()
    {
        if (state_recording_)
        {
            for (size_t i = 0; i < bodies_.size(); ++i)
            {
                if (bodies_[i]->checkNewlyUpdated())
                {
                    BaseParticles &base_particles = bodies_[i]->getBaseParticles();
                    base_particles.dvParticlePosition()->prepareForOutput(ExecutionPolicy{});
                    prepare_variable_to_write_(base_particles.VariablesToWrite(), ExecutionPolicy{});
                }
            }
            BodyStatesRecordingToVtp::writeToFile();
        }
    };

  protected:
    OperationOnDataAssemble<ParticleVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_write_;
};
} // namespace SPH
#endif // IO_VTK_MESH_CK_H
