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
 * @file 	io_observation_ck.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef IO_OBSERVATION_CK_H
#define IO_OBSERVATION_CK_H

#include "io_base.h"

#include "execution_policy.h"
#include "interpolation_dynamics.hpp"
#include "io_plt.h"

namespace SPH
{
/**
 * @class ObservedQuantityRecordingCK
 * @brief write files for observed quantity
 */
template <class ExecutionPolicy, typename DataType>
class ObservedQuantityRecordingCK : public BodyStatesRecording,
                                    public ObservingAQuantityCK<ExecutionPolicy, DataType>
{
  protected:
    SPHBody &observer_;
    PltEngine plt_engine_;
    BaseParticles &base_particles_;
    std::string dynamics_identifier_name_;
    const std::string quantity_name_;
    std::string filefullpath_output_;

  public:
    DataType type_indicator_; /*< this is an indicator to identify the variable type. */

  public:
    ObservedQuantityRecordingCK(const std::string &quantity_name, ContactRelation &contact_relation)
        : BodyStatesRecording(contact_relation.getSPHBody()),
          ObservingAQuantityCK<ExecutionPolicy, DataType>(contact_relation, quantity_name),
          observer_(contact_relation.getSPHBody()), plt_engine_(),
          base_particles_(observer_.getBaseParticles()),
          dynamics_identifier_name_(contact_relation.getSPHBody().getName()),
          quantity_name_(quantity_name)
    {
        DataType *interpolated_quantities_ = this->dv_interpolated_quantities_->DataField();
        /** Output for .dat file. */
        filefullpath_output_ = io_environment_.output_folder_ + "/" + dynamics_identifier_name_ + "_" + quantity_name + ".dat";
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << "run_time"
                 << "   ";

        for (size_t i = 0; i != base_particles_.TotalRealParticles(); ++i)
        {
            std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
            plt_engine_.writeAQuantityHeader(out_file, interpolated_quantities_[i], quantity_name_i);
        }
        out_file << "\n";
        out_file.close();
    };
    virtual ~ObservedQuantityRecordingCK(){};

    virtual void writeWithFileName(const std::string &sequence) override
    {
        this->exec();
        //        dv_interpolated_quantities_->synchronizeWithDevice(ExecutionPolicy{});
        DataType *interpolated_quantities_ = this->dv_interpolated_quantities_->DataField();
        std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
        out_file << sv_physical_time_.getValue() << "   ";
        for (size_t i = 0; i != base_particles_.TotalRealParticles(); ++i)
        {
            plt_engine_.writeAQuantity(out_file, interpolated_quantities_[i]);
        }
        out_file << "\n";
        out_file.close();
    };

    DataType *getObservedQuantity()
    {
        return this->interpolated_quantities_;
    }
};

} // namespace SPH
#endif // IO_OBSERVATION_CK_H
