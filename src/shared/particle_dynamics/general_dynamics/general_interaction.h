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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	general_interaction.h
 * @brief 	This is the interaction dynamics applicable for all type bodies
 * @author	Yaru Ren and Xiangyu Hu
 */

#ifndef GENERAL_INTERACTION_H
#define GENERAL_INTERACTION_H

#include "base_general_dynamics.h"

namespace SPH
{
class KernelCorrectionMatrixInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    explicit KernelCorrectionMatrixInner(BaseInnerRelation &inner_relation, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrixInner(){};

  protected:
    Real alpha_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class KernelCorrectionMatrixComplex : public KernelCorrectionMatrixInner, public GeneralDataDelegateContactOnly
{
  public:
    explicit KernelCorrectionMatrixComplex(ComplexRelation &complex_relation, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrixComplex(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;

    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class KernelGradientCorrectionInner
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient
 */
class KernelGradientCorrectionInner : public LocalDynamics, public GeneralDataDelegateInner
{
    ParticlesPairAverageInner<Matd> average_correction_matrix_;

  public:
    explicit KernelGradientCorrectionInner(KernelCorrectionMatrixInner &kernel_correction_inner);
    virtual ~KernelGradientCorrectionInner(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    template <class PairAverageType>
    void correctKernelGradient(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i)
    {
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            size_t index_j = neighborhood.j_[n];
            Vecd displacement = neighborhood.r_ij_[n] * neighborhood.e_ij_[n];

            Vecd corrected_direction = average_correction_matrix(index_i, index_j) * neighborhood.e_ij_[n];
            Real direction_norm = corrected_direction.norm();
            neighborhood.dW_ijV_j_[n] *= direction_norm;
            neighborhood.e_ij_[n] = corrected_direction / (direction_norm + Eps);
            neighborhood.r_ij_[n] = displacement.dot(neighborhood.e_ij_[n]);
        }
    };
};

/**
 * @class KernelGradientCorrectionComplex
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient in complex topology
 */
class KernelGradientCorrectionComplex : public KernelGradientCorrectionInner, public GeneralDataDelegateContactOnly
{
    StdVec<ParticlesPairAverageContact<Matd>> contact_average_correction_matrix_;

  public:
    KernelGradientCorrectionComplex(KernelCorrectionMatrixComplex &kernel_correction_complex);
    virtual ~KernelGradientCorrectionComplex(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ParticleSmoothing
 * @brief computing smoothed variable field by averaging with neighbors
 */
template <typename VariableType>
class ParticleSmoothing : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    explicit ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name)
        : LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
          W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)),
          smoothed_(*particles_->template getVariableByName<VariableType>(variable_name))
    {
        particles_->registerVariable(temp_, variable_name + "_temp");
    }

    virtual ~ParticleSmoothing(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real weight = W0_;
        VariableType summation = W0_ * smoothed_[index_i];
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            summation += inner_neighborhood.W_ij_[n] * smoothed_[index_j];
            weight += inner_neighborhood.W_ij_[n];
        }
        temp_[index_i] = summation / (weight + TinyReal);
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        smoothed_[index_i] = temp_[index_i];
    };

  protected:
    const Real W0_;
    StdLargeVec<VariableType> &smoothed_, temp_;
};
} // namespace SPH
#endif // GENERAL_INTERACTION_H
