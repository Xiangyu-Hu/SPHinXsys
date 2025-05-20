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
 * @brief 	Classes for observation recording.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#ifndef IO_OBSERVATION_H
#define IO_OBSERVATION_H

#include "io_plt.hpp"
#include "general_interpolation.h" 

namespace SPH
{

class BaseQuantityRecording : public BaseIO
{
  public:
    BaseQuantityRecording(SPHSystem &sph_system,
                          const std::string &dynamics_identifier_name);
    void setFullPath(const std::string &quantity_name);

  protected:
    PltEngine plt_engine_;
    std::string quantity_name_;
    std::string dynamics_identifier_name_;
    std::string filefullpath_output_;
};

template <typename...>
class ObservedQuantityRecording;

template <typename DataType>
class ObservedQuantityRecording<DataType> : public BaseQuantityRecording
{
  protected:
    SPHBody &observer_;
    BaseParticles &base_particles_;
    ObservingAQuantity<DataType> observation_method_;
    DiscreteVariable<DataType> *dv_interpolated_quantities_;
    size_t number_of_observe_;

  public:
    DataType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    ObservedQuantityRecording(const std::string &quantity_name, BaseContactRelation &contact_relation)
        : BaseQuantityRecording(contact_relation.getSPHBody().getSPHSystem(),
                                contact_relation.getSPHBody().getName()),
          observer_(contact_relation.getSPHBody()),
          base_particles_(observer_.getBaseParticles()),
          observation_method_(contact_relation, quantity_name),
          dv_interpolated_quantities_(observation_method_.dvInterpolatedQuantities()),
          number_of_observe_(base_particles_.TotalRealParticles())
    {
        setFullPath(quantity_name);
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "run_time" << "   ";
        DataType *interpolated_quantities = getObservedQuantity();
        for (size_t i = 0; i != number_of_observe_; ++i)
        {
            std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
            plt_engine_.writeAQuantityHeader(out_file, interpolated_quantities[i], quantity_name_i);
        }
        out_file << "\n";
        out_file.close();
    };
    virtual ~ObservedQuantityRecording() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_->getValue() << "   ";
        observation_method_.exec();
        DataType *interpolated_quantities = getObservedQuantity();
        for (size_t i = 0; i != number_of_observe_; ++i)
        {
            plt_engine_.writeAQuantity(out_file, interpolated_quantities[i]);
        }
        out_file << "\n";
        out_file.close();
    };

    DataType *getObservedQuantity()
    {
        return this->dv_interpolated_quantities_->Data();
    };

    size_t NumberOfObservedQuantity()
    {
        return number_of_observe_;
    };
};

template <typename...>
class ReducedQuantityRecording;
/**
 * @class ReducedQuantityRecording
 * @brief write reduced quantity of a body
 */
template <class LocalReduceMethodType>
class ReducedQuantityRecording<LocalReduceMethodType> : public BaseQuantityRecording
{
  public:
    using VariableType = typename LocalReduceMethodType::FinishDynamics::OutputType;
    VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

  protected:
    ReduceDynamics<LocalReduceMethodType> reduce_method_;
    VariableType reduced_quantity_;

  public:
    template <class DynamicsIdentifier, typename... Args>
    ReducedQuantityRecording(DynamicsIdentifier &identifier, Args &&...args)
        : BaseQuantityRecording(identifier.getSPHBody().getSPHSystem(),
                                identifier.getName()),
          reduce_method_(identifier, std::forward<Args>(args)...),
          reduced_quantity_(ZeroData<VariableType>::value)
    {
        quantity_name_ = reduce_method_.QuantityName();
        setFullPath(quantity_name_);
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "\"run_time\"" << "   ";
        plt_engine_.writeAQuantityHeader(out_file, reduced_quantity_, quantity_name_);
        out_file << "\n";
        out_file.close();
    };
    virtual ~ReducedQuantityRecording() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_->getValue() << "   ";
        reduced_quantity_ = reduce_method_.exec();
        plt_engine_.writeAQuantity(out_file, reduced_quantity_);
        out_file << "\n";
        out_file.close();
    };

    VariableType *getObservedQuantity()
    {
        return &reduced_quantity_;
    };

    size_t NumberOfObservedQuantity()
    {
        return 1;
    };
};

template <typename DataType>
class SingularVariableRecording : public BaseQuantityRecording
{
  protected:
    SingularVariable<DataType> *variable_;

  public:
    typedef DataType VariableType;
    VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    SingularVariableRecording(SPHSystem &sph_system, SingularVariable<DataType> *variable)
        : BaseQuantityRecording(sph_system, "SingularVariable"),
          variable_(variable)
    {
        setFullPath(variable->Name());
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "\"run_time\"" << "   ";
        plt_engine_.writeAQuantityHeader(out_file, variable_->getValue(), quantity_name_);
        out_file << "\n";
        out_file.close();
    };

    virtual ~SingularVariableRecording() {};

    virtual void writeToFile(size_t iteration_step = 0) override
    {
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_->getValue() << "   ";
        plt_engine_.writeAQuantity(out_file, variable_->getValue());
        out_file << "\n";
        out_file.close();
    };
};
} // namespace SPH
#endif // IO_OBSERVATION_H