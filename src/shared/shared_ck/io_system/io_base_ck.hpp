#ifndef IO_BASE_CK_HPP
#define IO_BASE_CK_HPP

#include "io_base_ck.h"

#include "general_reduce_ck.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
void BodyStatesRecordingToVtpCK<ExecutionPolicy>::prepareToWrite()
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
//=================================================================================================//
template <class ExecutionPolicy>
void BodyStatesRecordingToVtpCK<ExecutionPolicy>::writeToFile()
{
    if (state_recording_)
    {
        prepareToWrite();
        BodyStatesRecordingToVtp::writeToFile();
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyStatesRecordingToVtpCK<ExecutionPolicy>::writeToFile(size_t iteration_step)
{
    if (state_recording_)
    {
        prepareToWrite();
        BodyStatesRecordingToVtp::writeToFile(iteration_step);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename DerivedVariableMethod, typename DynamicsIdentifier, typename... Args>
BodyStatesRecording &BodyStatesRecordingToVtpCK<ExecutionPolicy>::addDerivedVariableToWrite(
    DynamicsIdentifier &identifier, Args &&...args)
{
    SPHBody &sph_body = identifier.getSPHBody();
    if (isBodyIncluded(bodies_, &sph_body))
    {
        derived_variables_.push_back(
            derived_variables_keeper_.createPtr<StateDynamics<ParallelPolicy, DerivedVariableMethod>>(
                identifier, std::forward<Args>(args)...));
    }
    else
    {
        std::cout << "\n Error: the body:" << sph_body.getName()
                  << " is not in the recording body list" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *this;
}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename... Args>
RestartIOCK<ExecutionPolicy>::RestartIOCK(Args &&...args)
    : RestartIO(std::forward<Args>(args)...)
{
    if (this->summary_enabled_)
    {
        for (int k = 0; k < 3; ++k)
        {
            output_evolving_variables_bounds_[k].resize(real_bodies_.size());
            evolving_variables_names_[k].resize(real_bodies_.size());
        }

        for (size_t i = 0; i < real_bodies_.size(); ++i)
        {
            BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
            DiscreteVariables &evolving_variables = base_particles.EvolvingVariables();

            // scalar bounds
            constexpr int type_index_Real = DataTypeIndex<Real>::value;
            for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(evolving_variables))
            {
                evolving_variables_names_[0][i].push_back(variable->Name());
                output_evolving_variables_bounds_[0][i].push_back(
                    particle_dynamics_keeper_.template createPtr<
                        ReduceDynamicsCK<ExecutionPolicy, MaximumNorm<Real>>>(*real_bodies_[i], variable->Name()));
            }
            // vectors bounds
            constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
            for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(evolving_variables))
            {
                evolving_variables_names_[1][i].push_back(variable->Name());
                output_evolving_variables_bounds_[1][i].push_back(
                    particle_dynamics_keeper_.template createPtr<
                        ReduceDynamicsCK<ExecutionPolicy, MaximumNorm<Vecd>>>(*real_bodies_[i], variable->Name()));
            }

            // matrix bounds
            constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
            for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(evolving_variables))
            {
                evolving_variables_names_[2][i].push_back(variable->Name());
                output_evolving_variables_bounds_[2][i].push_back(
                    particle_dynamics_keeper_.template createPtr<
                        ReduceDynamicsCK<ExecutionPolicy, MaximumNorm<Matd>>>(*real_bodies_[i], variable->Name()));
            }
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void RestartIOCK<ExecutionPolicy>::writeToFile(size_t iteration_step)
{
    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
        prepare_variable_to_write_(base_particles.EvolvingVariables(), ExecutionPolicy{});
    }
    RestartIO::writeToFile(iteration_step);

    if (summary_enabled_)
    {
        reportEvolvingVariablesBounds(iteration_step);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void RestartIOCK<ExecutionPolicy>::reportEvolvingVariablesBounds(size_t restart_step)
{
    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        std::string body_name = real_bodies_[i]->getName();

        std::cout << "Evolving Variables Summary for " << body_name << ":\n";
        std::cout << "---------------------------------------------\n";
        for (UnsignedInt j = 0; j < output_evolving_variables_bounds_[0][i].size(); ++j)
        {
            std::pair<Real, UnsignedInt> bound = output_evolving_variables_bounds_[0][i][j]->exec();
            std::cout << std::fixed << std::setprecision(9)
                      << "Evolving scalar variable bound: " << evolving_variables_names_[0][i][j]
                      << " = " << bound.first << ", particle_index = " << bound.second << "\n";
        }
        std::cout << "---------------------------------------------\n";
        for (UnsignedInt j = 0; j < output_evolving_variables_bounds_[1][i].size(); ++j)
        {
            std::pair<Real, UnsignedInt> bound = output_evolving_variables_bounds_[1][i][j]->exec();
            std::cout << std::fixed << std::setprecision(9)
                      << "Evolving vector variable bound: " << evolving_variables_names_[1][i][j]
                      << " = " << bound.first << ", particle_index = " << bound.second << "\n";
        }
        std::cout << "---------------------------------------------\n";
        for (UnsignedInt j = 0; j < output_evolving_variables_bounds_[2][i].size(); ++j)
        {
            std::pair<Real, UnsignedInt> bound = output_evolving_variables_bounds_[2][i][j]->exec();
            std::cout << std::fixed << std::setprecision(9)
                      << "Evolving matrix variable bound: " << evolving_variables_names_[2][i][j]
                      << " = " << bound.first << ", particle_index = " << bound.second << "\n";
        }
        std::cout << "---------------------------------------------\n";
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void RestartIOCK<ExecutionPolicy>::readFromFile(size_t iteration_step)
{
    RestartIO::readFromFile(iteration_step);

    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
        finalize_variables_after_read_(base_particles.EvolvingVariables(), ExecutionPolicy{});
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void ReloadParticleIOCK<ExecutionPolicy>::writeToFile(size_t iteration_step)
{
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        BaseParticles &base_particles = bodies_[i]->getBaseParticles();
        prepare_variable_to_reload_(base_particles.EvolvingVariables(), ExecutionPolicy{});
    }
    ReloadParticleIO::writeToFile(iteration_step);
}
//=================================================================================================//
} // namespace SPH
#endif // IO_BASE_CK_HPP
