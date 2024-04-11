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
 * @file 	general_refinement.h
 * @brief 	TODO: these classes should be improved on easy understanding and proper memory usage.
 * @author	Yijie Sun and Xiangyu Hu
 */
#ifndef GENERAL_REFINEMENT_H
#define GENERAL_REFINEMENT_H

#include "base_general_dynamics.h"
#include <array>

namespace SPH
{
class LevelSetComplexShape;
class ParticleSplitAndMerge;

/**
 * @class ComputeDensityErrorInner
 * @brief compute error of particle splitting and merging
 */
class ComputeDensityErrorInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    ComputeDensityErrorInner(BaseInnerRelation &inner_relation);
    virtual ~ComputeDensityErrorInner(){};

    Vecd getPositionFromDensityError(const StdVec<size_t> &original_indices, const StdVec<Vecd> &initial_new_positions,
                                     const StdVec<size_t> &new_indices, Real min_distance, Real max_distance);
    virtual void initializeDensityError();

    StdLargeVec<Real> density_error_;
    StdLargeVec<bool> tag_split_;

  protected:
    ParticleSplitAndMerge &particle_adaptation_;
    Real rho0_, inv_sigma0_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &h_ratio_, &Vol_, &rho_, &mass_;
    Vecd E_cof_ = Vecd::Zero();
    Real sigma_E_ = 0.0;
    Real E_cof_sigma_ = 0.0;
    StdVec<Vecd> grad_new_indices_;
    StdVec<Vecd> dW_new_indices_;
    StdVec<Real> sign_new_indices_;

    virtual Vecd computeKernelGradient(size_t index_rho);
    virtual Real computeNewGeneratedParticleDensity(size_t index_rho, const Vecd &position);
    virtual Vecd getPosition(const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions, const StdVec<size_t> &new_indices);
    virtual void densityErrorOfNewGeneratedParticles(const StdVec<size_t> &new_indices, const StdVec<Vecd> &new_positions);
    virtual void densityErrorOfNeighborParticles(const StdVec<size_t> &new_indices, const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions);
    virtual Real computeKernelWeightBetweenParticles(Real h_ratio, Vecd displacement, Real Vol_ratio = 1.0);
    virtual Vecd computeKernelWeightGradientBetweenParticles(Real h_ratio_min, Vecd displacement, Real Vol);
    virtual void computeDensityErrorOnNeighborParticles(Neighborhood &neighborhood, size_t index_rho,
                                                        const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions);
    virtual Vecd positionLimitation(Vecd displacement, Real min_distance, Real max_distance);
};

/**
 * @class ComputeDensityErrorWithWall
 * @brief compute error of particle splitting and merging
 */
class ComputeDensityErrorWithWall : public ComputeDensityErrorInner, public GeneralDataDelegateContactOnly
{
  public:
    ComputeDensityErrorWithWall(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : ComputeDensityErrorInner(inner_relation), GeneralDataDelegateContactOnly(contact_relation)
    {
        for (size_t k = 0; k != contact_bodies_.size(); ++k)
        {
            contact_Vol_.push_back(&(contact_bodies_[k]->getBaseParticles().VolumetricMeasures()));
        }
    };
    virtual ~ComputeDensityErrorWithWall(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;

    virtual Vecd computeKernelGradient(size_t index_rho) override;
    virtual Real computeNewGeneratedParticleDensity(size_t index_rho, const Vecd &position) override;
    virtual void densityErrorOfNeighborParticles(const StdVec<size_t> &new_indices, const StdVec<size_t> &original_indices,
                                                 const StdVec<Vecd> &new_positions) override;
};

/**
 * @class ParticleRefinementWithPrescribedArea
 * @brief particle split in prescribed area.
 */
class ParticleRefinementWithPrescribedArea : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    ParticleRefinementWithPrescribedArea(SPHBody &sph_body, Shape &refinement_region);
    virtual ~ParticleRefinementWithPrescribedArea(){};

  protected:
    BoundingBox refinement_region_bounds_;
    ParticleSplitAndMerge &particle_adaptation_;
    Real inv_rho0_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &rho_;
    StdLargeVec<Real> &mass_;
    StdLargeVec<Real> &h_ratio_; /**< the ratio between reference smoothing length to variable smoothing length */

    virtual bool checkLocation(const BoundingBox &refinement_region_bounds, Vecd position, Real volume);
};

/**
 * @class ParticleRefinementWithPrescribedArea
 * @brief particle split in prescribed area.
 */
class ParticleSplitWithPrescribedArea : public ParticleRefinementWithPrescribedArea
{
  public:
    ParticleSplitWithPrescribedArea(SPHBody &sph_body, Shape &refinement_region);
    virtual ~ParticleSplitWithPrescribedArea(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        if (splitCriteria(index_i))
        {
            StdVec<size_t> new_indices;
            splittingModel(index_i, new_indices);
        }
    };

    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> total_split_error_;
    StdVec<Vecd> split_position_;
    StdVec<std::pair<size_t, size_t>> split_index_;
    size_t particle_number_change = 0;

    virtual void setupDynamics(Real dt) override;
    virtual bool splitCriteria(size_t index_i);
    virtual void splittingModel(size_t index_i, StdVec<size_t> &new_indices);
    virtual Vecd getSplittingPosition(const StdVec<size_t> &new_indices);
    virtual void updateNewlySplittingParticle(size_t index_i, size_t index_j, Vecd pos_split);
};

/**
 * @class SplitWithMinimumDensityErrorInner
 * @brief split particles with minimum density error.
 */
class SplitWithMinimumDensityErrorInner : public ParticleSplitWithPrescribedArea
{
  public:
    SplitWithMinimumDensityErrorInner(BaseInnerRelation &inner_relation, Shape &refinement_region)
        : ParticleSplitWithPrescribedArea(inner_relation.getSPHBody(), refinement_region),
          compute_density_error(inner_relation)
    {
        particles_->registerVariable(total_split_error_, "SplitDensityError", Real(0));
    };
    virtual ~SplitWithMinimumDensityErrorInner(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    ComputeDensityErrorInner compute_density_error;

    virtual Vecd getSplittingPosition(const StdVec<size_t> &new_indices) override;
    virtual void setupDynamics(Real dt) override;
};

/**
 * @class SplitWithMinimumDensityErrorWithWall
 * @brief split particles with minimum density error.
 */
class SplitWithMinimumDensityErrorWithWall : public SplitWithMinimumDensityErrorInner
{
  public:
    SplitWithMinimumDensityErrorWithWall(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation,
                                         Shape &refinement_region)
        : SplitWithMinimumDensityErrorInner(inner_relation, refinement_region),
          compute_density_error(inner_relation, contact_relation){};
    virtual ~SplitWithMinimumDensityErrorWithWall(){};

  protected:
    ComputeDensityErrorWithWall compute_density_error;
};

/**
 * @class ParticleMergeWithPrescribedArea
 * @brief merging particle for a body in prescribed area.
 */
class ParticleMergeWithPrescribedArea : public ParticleRefinementWithPrescribedArea,
                                        public DataDelegateInner<BaseParticles, DataDelegateEmptyBase>
{
  public:
    ParticleMergeWithPrescribedArea(BaseInnerRelation &inner_relation, Shape &refinement_region);
    virtual ~ParticleMergeWithPrescribedArea(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        if (!tag_merged_[index_i])
        {
            StdVec<size_t> merge_indices; // three particles for merging to two
            if (mergeCriteria(index_i, merge_indices))
            {
                merge_indices.push_back(index_i);
                tag_merged_[index_i] = true;
                mergingModel(merge_indices);
            }
        }
    };

  protected:
    ParticleData &all_particle_data_;
    StdLargeVec<Vecd> &vel_n_;
    StdLargeVec<bool> tag_merged_;
    StdLargeVec<Real> total_merge_error_;

    virtual void setupDynamics(Real dt) override;
    virtual void mergingModel(const StdVec<size_t> &merge_indices);
    virtual bool mergeCriteria(size_t index_i, StdVec<size_t> &merge_indices);
    bool findMergeParticles(size_t index_i, StdVec<size_t> &merge_indices, Real search_size, Real search_distance);
    virtual void updateMergedParticleInformation(size_t merged_index, const StdVec<size_t> &merge_indices);

    struct MergeParticleDataValue
    {
        template <typename DataType>
        void operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper,
                        size_t merged_index, const StdVec<size_t> &merge_indices, StdVec<Real> merge_mass)
        {
            Real total_mass = 0.0;
            for (size_t k = 0; k != merge_indices.size(); ++k)
                total_mass += merge_mass[k];

            for (size_t i = 0; i != data_keeper.size(); ++i)
            {
                DataType particle_data_temp = ZeroData<DataType>::value;
                for (size_t k = 0; k != merge_indices.size(); ++k)
                    particle_data_temp += merge_mass[k] * (*data_keeper[i])[merge_indices[k]];

                (*data_keeper[i])[merged_index] = particle_data_temp / (total_mass + TinyReal);
            }
        };
    };

    OperationOnDataAssemble<ParticleData, MergeParticleDataValue> merge_particle_value_;
};

/**
 * @class ParticleMergeWithPrescribedArea
 * @brief merging particles with minimum density error.
 */
class MergeWithMinimumDensityErrorInner : public ParticleMergeWithPrescribedArea
{
  public:
    MergeWithMinimumDensityErrorInner(BaseInnerRelation &inner_relation, Shape &refinement_region)
        : ParticleMergeWithPrescribedArea(inner_relation, refinement_region),
          compute_density_error(inner_relation){};
    virtual ~MergeWithMinimumDensityErrorInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        ParticleMergeWithPrescribedArea::interaction(index_i, dt);
    };

  protected:
    ComputeDensityErrorInner compute_density_error;
    StdLargeVec<Real> merge_error_;
    Real rotation = 0.0;
    size_t merge_change_number = 0;

    virtual void setupDynamics(Real dt) override;
    virtual void mergingModel(const StdVec<size_t> &merge_indices) override;
    virtual bool mergeCriteria(size_t index_i, StdVec<size_t> &merge_indices) override;
    virtual Vecd getMergingPosition(const StdVec<size_t> &new_indices, const StdVec<size_t> &merge_indices);
    virtual Real angularMomentumConservation(size_t index_center, const StdVec<size_t> &merge_indices);
    virtual void kineticEnergyConservation(const StdVec<size_t> &merge_indices);
    virtual void updateNewlyMergingParticle(size_t index_center, const StdVec<size_t> &new_indices, Vecd pos_split);
};

/**
 * @class ParticleMergeWithPrescribedArea
 * @brief merging particles with minimum density error.
 */
class MergeWithMinimumDensityErrorWithWall : public MergeWithMinimumDensityErrorInner
{
  public:
    MergeWithMinimumDensityErrorWithWall(BaseInnerRelation &inner_relation,
                                         BaseContactRelation &contact_relation, Shape &refinement_region)
        : MergeWithMinimumDensityErrorInner(inner_relation, refinement_region),
          compute_density_error(inner_relation, contact_relation){};
    virtual ~MergeWithMinimumDensityErrorWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        MergeWithMinimumDensityErrorInner::interaction(index_i, dt);
    };

  protected:
    ComputeDensityErrorWithWall compute_density_error;
};
} // namespace SPH
#endif // GENERAL_REFINEMENT_H
