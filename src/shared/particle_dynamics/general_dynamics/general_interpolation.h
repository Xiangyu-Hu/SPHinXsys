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

#include "base_general_dynamics.h"

namespace SPH
{
/**
 * @class BaseInterpolation
 * @brief Base class for interpolation.
 */
template <typename DataType>
class BaseInterpolation : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit BaseInterpolation(BaseContactRelation &contact_relation, const std::string &variable_name)
        : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
          dv_interpolated_quantities_(nullptr), interpolated_quantities_(nullptr)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
            DataType *contact_data =
                this->contact_particles_[k]->template getVariableDataByName<DataType>(variable_name);
            contact_data_.push_back(contact_data);
        }
    };
    virtual ~BaseInterpolation() {};
    DiscreteVariable<DataType> *dvInterpolatedQuantities() { return dv_interpolated_quantities_; };

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        DataType observed_quantity = ZeroData<DataType>::value;
        Real ttl_weight(0);

        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real *Vol_k = contact_Vol_[k];
            DataType *data_k = contact_data_[k];
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

                observed_quantity += weight_j * data_k[index_j];
                ttl_weight += weight_j;
            }
        }
        interpolated_quantities_[index_i] = observed_quantity / (ttl_weight + TinyReal);
    };

  protected:
    DiscreteVariable<DataType> *dv_interpolated_quantities_;
    DataType *interpolated_quantities_;
    StdVec<Real *> contact_Vol_;
    StdVec<DataType *> contact_data_;
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
        this->dv_interpolated_quantities_ =
            this->particles_->template getVariableByName<DataType>(interpolated_variable);
        this->interpolated_quantities_ = this->dv_interpolated_quantities_->Data();
    };
    virtual ~InterpolatingAQuantity() {};
};

/**
 * @class ObservingAQuantity
 * @brief Observing a variable from contact bodies.
 */
template <typename DataType>
class ObservingAQuantity : public InteractionDynamics<BaseInterpolation<DataType>>
{
  public:
    explicit ObservingAQuantity(BaseContactRelation &contact_relation, const std::string &variable_name)
        : InteractionDynamics<BaseInterpolation<DataType>>(contact_relation, variable_name)
    {
        this->dv_interpolated_quantities_ = this->particles_->template registerStateVariableOnly<DataType>(variable_name);
        this->interpolated_quantities_ = this->dv_interpolated_quantities_->Data();
    };
    virtual ~ObservingAQuantity() {};
};

/**
 * @class CorrectInterpolationKernelWeights
 * @brief  correct kernel weights for interpolation between general bodies
 * TODO: this formulation is not correct, need to be fixed.
 */
class CorrectInterpolationKernelWeights : public LocalDynamics,
                                          public DataDelegateContact
{
  public:
    explicit CorrectInterpolationKernelWeights(BaseContactRelation &contact_relation);
    virtual ~CorrectInterpolationKernelWeights() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd weight_correction = Vecd::Zero();
        Matd local_configuration = Eps * Matd::Identity();

        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];

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
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                contact_neighborhood.W_ij_[n] -= normalized_weight_correction.dot(contact_neighborhood.e_ij_[n]) *
                                                 contact_neighborhood.dW_ij_[n];
            }
        }
    };

  protected:
    StdVec<Real *> contact_Vol_;
};
} // namespace SPH
#endif // GENERAL_INTERPOLATION_H