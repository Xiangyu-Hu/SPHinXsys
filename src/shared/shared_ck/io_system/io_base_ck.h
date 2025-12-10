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
 * @file 	io_base_ck.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef IO_BASE_CK_H
#define IO_BASE_CK_H

#include "execution_policy.h"
#include "io_base.h"
#include "io_vtk.h"

namespace SPH
{

template <class ExecutionPolicy>
class BodyStatesRecordingToVtpCK : public BodyStatesRecordingToVtp
{
  protected:
    OperationOnDataAssemble<ParticleVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_write_;

    void prepareToWrite()
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
    }

  public:
    template <typename... Args>
    BodyStatesRecordingToVtpCK(Args &&...args) : BodyStatesRecordingToVtp(std::forward<Args>(args)...){};
    virtual ~BodyStatesRecordingToVtpCK() {};

    virtual void writeToFile() override
    {
        if (state_recording_)
        {
            prepareToWrite();
            BodyStatesRecordingToVtp::writeToFile();
        }
    };

    virtual void writeToFile(size_t iteration_step) override
    {
        if (state_recording_)
        {
            prepareToWrite();
            BodyStatesRecordingToVtp::writeToFile(iteration_step);
        }
    }
};

template <class ExecutionPolicy>
class RestartIOCK : public RestartIO
{
  public:
    template <typename... Args>
    RestartIOCK(Args &&...args) : RestartIO(std::forward<Args>(args)...){};
    virtual ~RestartIOCK() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        for (size_t i = 0; i < bodies_.size(); ++i)
        {
            BaseParticles &base_particles = bodies_[i]->getBaseParticles();
            prepare_variable_to_write_(base_particles.EvolvingVariables(), ExecutionPolicy{});
        }
        RestartIO::writeToFile(iteration_step);
    };

    virtual void readFromFile(size_t iteration_step = 0) override
    {
        RestartIO::readFromFile(iteration_step);

        for (size_t i = 0; i < bodies_.size(); ++i)
        {
            BaseParticles &base_particles = bodies_[i]->getBaseParticles();
            finalize_variables_after_read_(base_particles.EvolvingVariables(), ExecutionPolicy{});
        }
    };

  protected:
    OperationOnDataAssemble<ParticleVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_write_;
    OperationOnDataAssemble<ParticleVariables, FinalizeVariablesAfterRead<DiscreteVariable>> finalize_variables_after_read_;
};

template <class ExecutionPolicy>
class ReloadParticleIOCK : public ReloadParticleIO
{
  public:
    template <typename... Args>
    ReloadParticleIOCK(Args &&...args) : ReloadParticleIO(std::forward<Args>(args)...){};
    virtual ~ReloadParticleIOCK() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        for (size_t i = 0; i < bodies_.size(); ++i)
        {
            BaseParticles &base_particles = bodies_[i]->getBaseParticles();
            prepare_variable_to_reload_(base_particles.EvolvingVariables(), ExecutionPolicy{});
        }
        ReloadParticleIO::writeToFile(iteration_step);
    };

  protected:
    OperationOnDataAssemble<ParticleVariables, PrepareVariablesToWrite<DiscreteVariable>> prepare_variable_to_reload_;
};
} // namespace SPH
#endif // IO_BASE_CK_H
