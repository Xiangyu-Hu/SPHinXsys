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
 * @file 	general_interpolation.h
 * @brief 	There are the classes for interpolation algorithm.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef GENERAL_INTERPOLATION_H
#define GENERAL_INTERPOLATION_H

#include "general_dynamics.h"

namespace SPH
{
typedef DataDelegateContact<BaseParticles, BaseParticles> InterpolationContactData;

template<class DataType>
class BaseInterpolationKernel {

  public:
    BaseInterpolationKernel(NeighborhoodDevice **contact_configuration_device,
                            size_t contact_size, DataType *interpolated_quantities)
        : contact_configuration_(contact_configuration_device), contact_size_(contact_size),
          interpolated_quantities_(interpolated_quantities) {}

    void interaction(size_t index_i, Real dt = 0.0) {
        DataType observed_quantity(0);
        DeviceReal ttl_weight(0);

        for (size_t k = 0; k < contact_size_; ++k)
        {
            const auto *Vol_k = contact_Vol_[k];
            const auto *data_k = contact_data_[k];
            const auto &contact_neighborhood = contact_configuration_[k][index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size(); ++n)
            {
                const auto index_j = contact_neighborhood.j_[n];
                const auto weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

                observed_quantity += weight_j * data_k[index_j];
                ttl_weight += weight_j;
            }
        }
        interpolated_quantities_[index_i] = observed_quantity / (ttl_weight + TinyReal);
    }

    void setContactVol(DeviceReal **contact_Vol) { contact_Vol_ = contact_Vol; }
    void setContactData(DataType **contact_data) { contact_data_ = contact_data; }

  private:
    NeighborhoodDevice** contact_configuration_;
    DeviceReal** contact_Vol_;
    DataType** contact_data_;
    size_t contact_size_;
    DataType* interpolated_quantities_;
};

/**
 * @class BaseInterpolation
 * @brief Base class for interpolation.
 */
template <typename DataType, bool DeviceAllocation = false>
class BaseInterpolation : public LocalDynamics, public InterpolationContactData
{
  public:
    StdLargeVec<DataType> *interpolated_quantities_;

    explicit BaseInterpolation(BaseContactRelation &contact_relation, const std::string &variable_name)
        : LocalDynamics(contact_relation.getSPHBody()), InterpolationContactData(contact_relation),
          interpolated_quantities_(nullptr)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
            StdLargeVec<DataType> *contact_data =
                this->contact_particles_[k]->template getVariableByName<DataType>(variable_name);
            contact_data_.push_back(contact_data);
        }
    };
    virtual ~BaseInterpolation(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        DataType observed_quantity = ZeroData<DataType>::value;
        Real ttl_weight(0);

        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
            StdLargeVec<DataType> &data_k = *(contact_data_[k]);
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

                observed_quantity += weight_j * data_k[index_j];
                ttl_weight += weight_j;
            }
        }
        (*interpolated_quantities_)[index_i] = observed_quantity / (ttl_weight + TinyReal);
    };

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<DataType> *> contact_data_;
};

template <typename DataType>
class BaseInterpolation<DataType, true> : public LocalDynamics, public InterpolationContactData,
                                          public DeviceExecutable<BaseInterpolation<DataType, true>,
                                              BaseInterpolationKernel<typename DataTypeEquivalence<DataType>::device_t>>
{
  public:
    StdLargeVec<DataType> *interpolated_quantities_;

    explicit BaseInterpolation(BaseContactRelation &contact_relation, const std::string &variable_name)
        : LocalDynamics(contact_relation.getSPHBody()), InterpolationContactData(contact_relation),
          DeviceExecutable<BaseInterpolation<DataType, true>,
                           BaseInterpolationKernel<typename DataTypeEquivalence<DataType>::device_t>>(
              this, this->contact_configuration_device_->data(), this->contact_configuration_device_->size(),
              particles_->registerDeviceVariable<typename DataTypeEquivalence<DataType>::device_t>(
                          variable_name, particles_->total_real_particles_)),
          interpolated_quantities_(nullptr)
    {
            using DataTypeDevice = typename DataTypeEquivalence<DataType>::device_t;

            contact_Vol_device_ = makeSharedDevice<StdSharedVec<DeviceReal*>>(this->contact_particles_.size(),
                                                                               execution::executionQueue.getQueue());
            contact_data_device_ = makeSharedDevice<StdSharedVec<DataTypeDevice*>>(this->contact_particles_.size(),
                                                                                    execution::executionQueue.getQueue());

            for (size_t k = 0; k != this->contact_particles_.size(); ++k)
            {
                contact_Vol_device_->at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceReal>("Volume");
                DataTypeDevice *contact_data =
                    this->contact_particles_[k]->template getDeviceVariableByName<DataTypeDevice>(variable_name);
                contact_data_device_->at(k) = contact_data;
            }

            // Set device contact volume and data that have just been initialized
            auto* device_kernel = this->getDeviceProxy().getKernel();
            device_kernel->setContactVol(contact_Vol_device_->data());
            device_kernel->setContactData(contact_data_device_->data());
    };
    virtual ~BaseInterpolation(){};

  protected:
    SharedPtr<StdSharedVec<DeviceReal*>> contact_Vol_device_;
    SharedPtr<StdSharedVec<typename DataTypeEquivalence<DataType>::device_t*>> contact_data_device_;
};

/**
 * @class InterpolatingAQuantity
 * @brief Interpolate a given member data in the particles of a general body
 */
template <typename DataType>
class InterpolatingAQuantity : public BaseInterpolation<DataType>
{
  public:
    explicit InterpolatingAQuantity(BaseContactRelation &contact_relation,
                                    const std::string &interpolated_variable, const std::string &target_variable)
        : BaseInterpolation<DataType>(contact_relation, target_variable)
    {
        this->interpolated_quantities_ =
            this->particles_->template getVariableByName<DataType>(interpolated_variable);
    };
    virtual ~InterpolatingAQuantity(){};
};

/**
 * @class ObservingAQuantity
 * @brief Observing a variable from contact bodies.
 */
template <typename DataType, typename ExecutionPolicy = ParallelPolicy>
class ObservingAQuantity :
    public InteractionDynamics<BaseInterpolation<DataType, std::conditional_t<std::is_same_v<ExecutionPolicy, ParallelSYCLDevicePolicy>,
                                                                              std::true_type, std::false_type>::value>, ExecutionPolicy>
{
  public:
    explicit ObservingAQuantity(BaseContactRelation &contact_relation, const std::string &variable_name)
        : InteractionDynamics<BaseInterpolation<DataType, std::conditional_t<std::is_same_v<ExecutionPolicy, ParallelSYCLDevicePolicy>,
            std::true_type, std::false_type>::value>, ExecutionPolicy>(contact_relation, variable_name)
    {
        this->interpolated_quantities_ = registerObservedQuantity(variable_name);
    };
    virtual ~ObservingAQuantity(){};

  protected:
    StdLargeVec<DataType> observed_quantities_;

    /** Register the  observed variable if the variable name is new.
     * If the variable is registered already, the registered variable will be returned. */
    StdLargeVec<DataType> *registerObservedQuantity(const std::string &variable_name)
    {
        BaseParticles *particles = this->particles_;
        DiscreteVariable<DataType> *variable = findVariableByName<DataType>(particles->AllDiscreteVariables(), variable_name);
        if (variable == nullptr)
        {
            particles->registerVariable(observed_quantities_, variable_name, ZeroData<DataType>::value);
            return &observed_quantities_;
        }
        return particles->getVariableByName<DataType>(variable_name);
    };
};

/**
 * @class CorrectInterpolationKernelWeights
 * @brief  correct kernel weights for interpolation between general bodies
 */
class CorrectInterpolationKernelWeights : public LocalDynamics,
                                          public InterpolationContactData
{
  public:
    explicit CorrectInterpolationKernelWeights(BaseContactRelation &contact_relation);
    virtual ~CorrectInterpolationKernelWeights(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd weight_correction = Vecd::Zero();
        Matd local_configuration = Eps * Matd::Identity();

        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                Vecd gradW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

                weight_correction += weight_j * r_ji;
                local_configuration += r_ji * gradW_ijV_j.transpose();
            }
        }

        // correction matrix for interacting configuration
        Matd B = local_configuration.inverse();
        Vecd normalized_weight_correction = B * weight_correction;
        // Add the kernel weight correction to W_ij_ of neighboring particles.
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                contact_neighborhood.W_ij_[n] -= normalized_weight_correction.dot(contact_neighborhood.e_ij_[n]) *
                                                 contact_neighborhood.dW_ijV_j_[n] / Vol_k[index_j];
            }
        }
    };

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
};
} // namespace SPH
#endif // GENERAL_INTERPOLATION_H
