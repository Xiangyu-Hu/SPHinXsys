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
 * @file transport_velocity_correction.h
 * @brief The particle positions are corrected for more uniformed distribution
 * when there is negative pressure in the flow.
 * @details Note that the default coefficient is for using the dual time criteria method:
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * C Zhang, M Rezavand, X Hu - Journal of Computational Physics,
 * Volume 404, 1 March 2020, 109135.
 * If single (acoustic) time step is used, the coefficient should be decrease
 * to about 1/4 of the default value.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef TRANSPORT_VELOCITY_CORRECTION_H
#define TRANSPORT_VELOCITY_CORRECTION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
    template <typename ErrorDataType, typename ParameterADataType, typename ParameterCDataType>
    struct ErrorAndParameters
    {
        ErrorDataType error_;
        ParameterADataType a_;
        ParameterCDataType c_;

        ErrorAndParameters<ErrorDataType, ParameterADataType, ParameterCDataType>() :
            error_(ZeroData<ErrorDataType>::value),
            a_(ZeroData<ParameterADataType>::value),
            c_(ZeroData<ParameterCDataType>::value) {};
    };

/**
 * @class TransportVelocityCorrectionInner
 * @brief The particle positions are corrected for more uniformed distribution
 * when there is negative pressure in the flow.
 * @details Note that the default coefficient is for using the dual time criteria method:
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * C Zhang, M Rezavand, X Hu - Journal of Computational Physics,
 * Volume 404, 1 March 2020, 109135.
 * If single (acoustic) time step is used, the coefficient should be decrease
 * to about 1/4 of the default value.
 */

template <class ParticleScopeType>
class TransportVelocityCorrectionInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInner(BaseInnerRelation &inner_relation, Real coefficient = 0.2)
        : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
          pos_(particles_->pos_), indicator_(*particles_->getVariableByName<int>("Indicator")),
          smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
          coefficient_(coefficient), checkWithinScope(particles_){};
    virtual ~TransportVelocityCorrectionInner(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (checkWithinScope(index_i))
        {
            Vecd acceleration_trans = Vecd::Zero();
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

                // acceleration for transport velocity
                acceleration_trans -= 2.0 * nablaW_ijV_j;
            }

            pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
        }
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
    ParticleScopeType checkWithinScope;
};

/**
 * @class TransportVelocityCorrectionComplex
 * @brief  transport velocity correction considering the contribution from contact bodies
 */
template <class ParticleScopeType>
class TransportVelocityCorrectionComplex
    : public BaseInteractionComplex<TransportVelocityCorrectionInner<ParticleScopeType>, FluidContactOnly>
{
  public:
    template <typename... Args>
    TransportVelocityCorrectionComplex(Args &&...args)
        : BaseInteractionComplex<TransportVelocityCorrectionInner<ParticleScopeType>, FluidContactOnly>(
              std::forward<Args>(args)...){};
    virtual ~TransportVelocityCorrectionComplex(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        TransportVelocityCorrectionInner<ParticleScopeType>::interaction(index_i, dt);

        if (this->checkWithinScope(index_i))
        {
            Vecd acceleration_trans = Vecd::Zero();
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

                    // acceleration for transport velocity
                    acceleration_trans -= 2.0 * nablaW_ijV_j;
                }
            }

            this->pos_[index_i] += this->coefficient_ * this->smoothing_length_sqr_ * acceleration_trans;
        }
    };
};

/**
 * @class TransportVelocityCorrectionInnerAdaptive
 * @brief transport velocity correction
 */
class TransportVelocityCorrectionInnerAdaptive : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionInnerAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
};

/**
 * @class TransportVelocityCorrectionComplex<AllParticles>Adaptive
 * @brief  transport velocity correction considering the contribution from contact bodies
 */
class TransportVelocityCorrectionComplexAdaptive
    : public BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactOnly>
{
  public:
    template <typename... Args>
    TransportVelocityCorrectionComplexAdaptive(Args &&...args)
        : BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactOnly>(
              std::forward<Args>(args)...){};
    virtual ~TransportVelocityCorrectionComplexAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class TransportVelocityConsistencyCorrectionInner
 */
template <class ParticleScopeType>
class TransportVelocityConsistencyInner : public LocalDynamics, public FluidDataInner
{
public:
    explicit TransportVelocityConsistencyInner(BaseInnerRelation& inner_relation, Real coefficient = 0.2)
        : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
        pos_(particles_->pos_), B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")),
        indicator_(*particles_->getVariableByName<int>("Indicator")),
        smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
        coefficient_(coefficient), checkWithinScope(particles_) {};
    virtual ~TransportVelocityConsistencyInner() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (checkWithinScope(index_i))
        {
            Vecd acceleration_trans = Vecd::Zero();
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];;
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

                // acceleration for transport velocity
                acceleration_trans -= (B_[index_i] + B_[index_j]) * nablaW_ijV_j;
            }

            pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
        }
    };

protected:
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<int>& indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
    ParticleScopeType checkWithinScope;
};

template <class ParticleScopeType>
class TransportVelocityConsistencyComplex
    : public BaseInteractionComplex<TransportVelocityConsistencyInner<ParticleScopeType>, FluidContactWall>
{
public:
    template <typename... Args>
    TransportVelocityConsistencyComplex(Args &&...args)
        : BaseInteractionComplex<TransportVelocityConsistencyInner<ParticleScopeType>, FluidContactWall>(
            std::forward<Args>(args)...) 
    {
        for (size_t k = 0; k != FluidContactWall::contact_particles_.size(); ++k)
        {
            wall_B_.push_back(&(FluidContactWall::contact_particles_[k]->B_));
        }
    };
    virtual ~TransportVelocityConsistencyComplex() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        TransportVelocityConsistencyInner<ParticleScopeType>::interaction(index_i, dt);

        if (this->checkWithinScope(index_i))
        {
            Vecd acceleration_trans = Vecd::Zero();
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                StdLargeVec<Matd>& wall_B_k = *(wall_B_[k]);
                Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

                    // acceleration for transport velocity
                    acceleration_trans -= (this->B_[index_i] + this->B_[index_i]) * nablaW_ijV_j;
                }
            }

            this->pos_[index_i] += this->coefficient_ * this->smoothing_length_sqr_ * acceleration_trans;
        }
    };

protected:
    StdVec<StdLargeVec<Matd>*> wall_B_;
};

/**
 * @class TransportVelocityConsistencyInnerImplicit
 */
template <class ParticleScopeType>
class TransportVelocityConsistencyInnerImplicit : public LocalDynamics, public FluidDataInner
{
public:
    explicit TransportVelocityConsistencyInnerImplicit(BaseInnerRelation& inner_relation, Real coefficient = 0.2)
        : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
        kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
        Vol_(particles_->Vol_), pos_(particles_->pos_), 
        B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")),
        indicator_(*particles_->getVariableByName<int>("Indicator")),
        smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 1)),
        coefficient_(coefficient), checkWithinScope(particles_) {};
    virtual ~TransportVelocityConsistencyInnerImplicit() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (checkWithinScope(index_i))
        {
            ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
            updateStates(index_i, dt, error_and_parameters);
        }
    };

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = smoothing_length_sqr_)
    {
        ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
        Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Matd parameter_b = (B_[index_i] + B_[index_j]) *
                inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) *
                Vol_[index_j] * dt * dt;

            error_and_parameters.error_ += (B_[index_i] + B_[index_j]) *
                inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] * dt * dt;
            error_and_parameters.a_ -= parameter_b;
            error_and_parameters.c_ += parameter_b * parameter_b;
        }

        Matd evolution = Matd::Identity();
        error_and_parameters.a_ -= evolution;
        return error_and_parameters;
    };
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
    {
        Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
        Vecd parameter_k = coefficient_ * parameter_l.inverse() * error_and_parameters.error_;

        pos_[index_i] += error_and_parameters.a_ * parameter_k;

        Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Matd parameter_b = (B_[index_i] + B_[index_j]) *
                inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
                kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) *
                Vol_[index_j] * dt * dt;
            pos_[index_j] -= parameter_b * parameter_k;
        }
    };
    
    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<int>& indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
    ParticleScopeType checkWithinScope;
};

/**
 * @class TransportVelocityConsistencyComplexImplicit
 */
template <class ParticleScopeType>
class TransportVelocityConsistencyComplexImplicit
    : public BaseInteractionComplex<TransportVelocityConsistencyInnerImplicit<ParticleScopeType>, FluidContactWall>
{
public:
    template <typename... Args>
    TransportVelocityConsistencyComplexImplicit(Args &&...args)
        : BaseInteractionComplex<TransportVelocityConsistencyInnerImplicit<ParticleScopeType>, FluidContactWall>(
            std::forward<Args>(args)...)
    {
        for (size_t k = 0; k != FluidContactWall::contact_particles_.size(); ++k)
        {
            wall_B_.push_back(&(FluidContactWall::contact_particles_[k]->B_));
            contact_Vol_.push_back(&(FluidContactWall::contact_particles_[k]->Vol_));
        }
    };
    virtual ~TransportVelocityConsistencyComplexImplicit() {};

    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = this->smoothing_length_sqr_)
    {
        ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = 
            TransportVelocityConsistencyInnerImplicit<ParticleScopeType>::computeErrorAndParameters(index_i, dt);

        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
            StdLargeVec<Matd>& B_k = *(wall_B_[k]);
            Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Matd parameter_b = (this->B_[index_i] + this->B_[index_i]) *
                    contact_neighborhood.e_ij_[n] * contact_neighborhood.e_ij_[n].transpose() *
                    this->kernel_->d2W(contact_neighborhood.r_ij_[n], contact_neighborhood.e_ij_[n]) *
                    Vol_k[index_j] * dt * dt;

                error_and_parameters.error_ += (this->B_[index_i] + this->B_[index_i]) *
                                                contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n] * dt * dt;
                error_and_parameters.a_ -= parameter_b;
            }
        }

        return error_and_parameters;
    };

protected:
    StdVec<StdLargeVec<Matd>*> wall_B_;
    StdVec<StdLargeVec<Real>*> contact_Vol_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_H
