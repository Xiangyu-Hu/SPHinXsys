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
 * @file 	io_observation.h
 * @brief 	Classes for input and output with vtk (Paraview) files.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#ifndef IO_OBSERVATION_H
#define IO_OBSERVATION_H

#include "io_base.h"

#include "io_plt.h"

namespace SPH
{
/**
 * @class ObservedQuantityRecording
 * @brief write files for observed quantity
 */
template <typename VariableType>
class ObservedQuantityRecording : public BodyStatesRecording,
                                  public ObservingAQuantity<VariableType>
{
  protected:
    SPHBody &observer_;
    PltEngine plt_engine_;
    BaseParticles &base_particles_;
    std::string dynamics_identifier_name_;
    const std::string quantity_name_;
    std::string filefullpath_output_;

  public:
    VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    ObservedQuantityRecording(const std::string &quantity_name, BaseContactRelation &contact_relation)
        : BodyStatesRecording(contact_relation.getSPHBody()),
          ObservingAQuantity<VariableType>(contact_relation, quantity_name),
          observer_(contact_relation.getSPHBody()), plt_engine_(),
          base_particles_(observer_.getBaseParticles()),
          dynamics_identifier_name_(contact_relation.getSPHBody().getName()),
          quantity_name_(quantity_name)
    {
        /** Output for .dat file. */
        filefullpath_output_ = io_environment_.output_folder_ + "/" + dynamics_identifier_name_ + "_" + quantity_name + ".dat";
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "run_time"
                 << "   ";
        for (size_t i = 0; i != base_particles_.TotalRealParticles(); ++i)
        {
            std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
            plt_engine_.writeAQuantityHeader(out_file, this->interpolated_quantities_[i], quantity_name_i);
        }
        out_file << "\n";
        out_file.close();
    };
    virtual ~ObservedQuantityRecording(){};

    virtual void writeWithFileName(const std::string &sequence) override
    {
        this->exec();
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_.getValue() << "   ";
        for (size_t i = 0; i != base_particles_.TotalRealParticles(); ++i)
        {
            plt_engine_.writeAQuantity(out_file, this->interpolated_quantities_[i]);
        }
        out_file << "\n";
        out_file.close();
    };

    VariableType *getObservedQuantity()
    {
        return this->interpolated_quantities_;
    }
};

/**
 * @class ReducedQuantityRecording
 * @brief write reduced quantity of a body
 */
template <class LocalReduceMethodType>
class ReducedQuantityRecording : public BaseIO
{
  protected:
    PltEngine plt_engine_;
    ReduceDynamics<LocalReduceMethodType> reduce_method_;
    std::string dynamics_identifier_name_;
    const std::string quantity_name_;
    std::string filefullpath_output_;

  public:
    /*< deduce variable type from reduce method. */
    using VariableType = typename LocalReduceMethodType::ReturnType;
    VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    template <class DynamicsIdentifier, typename... Args>
    ReducedQuantityRecording(DynamicsIdentifier &identifier, Args &&...args)
        : BaseIO(identifier.getSPHBody().getSPHSystem()), plt_engine_(),
          reduce_method_(identifier, std::forward<Args>(args)...),
          dynamics_identifier_name_(reduce_method_.DynamicsIdentifierName()),
          quantity_name_(reduce_method_.QuantityName())
    {
        /** output for .dat file. */
        filefullpath_output_ = io_environment_.output_folder_ + "/" + dynamics_identifier_name_ + "_" + quantity_name_ + ".dat";
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "\"run_time\""
                 << "   ";
        plt_engine_.writeAQuantityHeader(out_file, reduce_method_.Reference(), quantity_name_);
        out_file << "\n";
        out_file.close();
    };
    virtual ~ReducedQuantityRecording(){};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_.getValue() << "   ";
        plt_engine_.writeAQuantity(out_file, reduce_method_.exec());
        out_file << "\n";
        out_file.close();
    };
};
} // namespace SPH
#endif // IO_OBSERVATION_H
