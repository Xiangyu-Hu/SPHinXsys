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
 * @file 	adaptation.h
 * @brief 	Adaptation scheme for particle in multi-resolution scenario.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ADAPTATION_H
#define ADAPTATION_H

#include "base_data_type_package.h"
#include "base_kernel.h"
#include "base_mesh.hpp"
#include "execution_policy.h"
#include "sphinxsys_containers.h"

namespace SPH
{

class SPHSystem;
class Shape;
class BaseParticles;
class BodyRegionByCell;
class LevelSet;
class BaseCellLinkedList;

/**
 * @class SPHAdaptation
 * @brief Base class for all adaptations.
 * The base class defines essential global parameters. It is also used for single-resolution method.
 * In the constructor parameter, system_refinement_ratio defines the relation between present resolution to the system reference resolution.
 * The derived classes are defined for more complex adaptations.
 */
class SPHAdaptation
{
  protected:
    Real h_spacing_ratio_;         /**< ratio of reference kernel smoothing length to particle spacing */
    Real system_refinement_ratio_; /**< ratio of system resolution to body resolution, set to 1.0 by default */
    int local_refinement_level_;   /**< refinement level respect to reference particle spacing */
    Real spacing_ref_;             /**< reference particle spacing used to determine local particle spacing */
    Real h_ref_;                   /**< reference smoothing length */
    UniquePtr<Kernel> kernel_ptr_; /**< unique pointer of kernel function owned this class */
    Real sigma0_ref_;              /**< Reference number density dependent on h_spacing_ratio_ and kernel function */
    Real spacing_min_;             /**< minimum particle spacing determined by local refinement level */
    Real Vol_min_;                 /**< minimum particle volume measure determined by local refinement level */
    Real h_ratio_max_;             /**< the ratio between the reference smoothing length to the minimum smoothing length */

  public:
    explicit SPHAdaptation(Real resolution_ref, Real h_spacing_ratio = 1.3, Real system_refinement_ratio = 1.0);
    explicit SPHAdaptation(SPHSystem &sph_system, Real h_spacing_ratio = 1.3, Real system_refinement_ratio = 1.0);
    virtual ~SPHAdaptation() {};

    int LocalRefinementLevel() { return local_refinement_level_; };
    Real SmoothingLengthSpacingRatio() { return h_spacing_ratio_; };
    Real ReferenceSpacing() { return spacing_ref_; };
    Real MinimumSpacing() { return spacing_min_; };
    Real ReferenceSmoothingLength() { return h_ref_; };
    Real MinimumSmoothingLength() { return h_ref_ / h_ratio_max_; };
    Kernel *getKernel() { return kernel_ptr_.get(); };
    Real LatticeNumberDensity() { return sigma0_ref_; };
    Real NumberDensityScaleFactor(Real smoothing_length_ratio);
    virtual Real SmoothingLengthRatio(size_t particle_index_i) { return 1.0; };
    void resetAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio = 1.0);
    virtual void initializeAdaptationVariables(BaseParticles &base_particles) {};
    Real SmoothingLengthByLevel(int level) { return h_ref_ / pow(2.0, level); };
    DiscreteVariable<Real> *AdaptiveSmoothingLength(BaseParticles &base_particles);

    virtual UniquePtr<BaseCellLinkedList> createCellLinkedList(const BoundingBoxd &domain_bounds, BaseParticles &base_particles);
    UniquePtr<BaseCellLinkedList> createRefinedCellLinkedList(int level, const BoundingBoxd &domain_bounds, BaseParticles &base_particles);
    virtual UniquePtr<LevelSet> createLevelSet(Shape &shape, Real refinement_ratio);

    template <class KernelType, typename... Args>
    void resetKernel(Args &&...args)
    {
        kernel_ptr_.reset(new KernelType(h_ref_, std::forward<Args>(args)...));
        sigma0_ref_ = computeLatticeNumberDensity(Vecd());
    };

    class UnitSmoothingLengthRatio
    {
      public:
        explicit UnitSmoothingLengthRatio(SPHAdaptation &adaptation) {};
        Real operator()(size_t /*particle_index_i*/) { return 1.0; };
    };

    typedef UnitSmoothingLengthRatio SmoothingRatioType;

  protected:
    Real computeLatticeNumberDensity(Vec2d zero);
    Real computeLatticeNumberDensity(Vec3d zero);
    virtual Real MostRefinedSpacing(Real coarse_particle_spacing, int local_refinement_level);
    Real MostRefinedSpacingRegular(Real coarse_particle_spacing, int local_refinement_level);
};

/**
 * @class AdaptiveSmoothingLength
 * @brief Base class for particle with local refinement.
 * @details Different refinement strategies will be used in derived classes.
 */
class AdaptiveSmoothingLength : public SPHAdaptation
{
  public:
    Real *h_ratio_; /**< the ratio between reference smoothing length to variable smoothing length */
    int *level_;    /**< the mesh level of the particle */

    AdaptiveSmoothingLength(Real resolution_ref, Real h_spacing_ratio_, Real system_refinement_ratio, int local_refinement_level);
    AdaptiveSmoothingLength(SPHSystem &sph_system, Real h_spacing_ratio_, Real system_refinement_ratio, int local_refinement_level);
    virtual ~AdaptiveSmoothingLength() {};

    virtual Real SmoothingLengthRatio(size_t particle_index_i) override
    {
        return h_ratio_[particle_index_i];
    };

    virtual void initializeAdaptationVariables(BaseParticles &base_particles) override;
    virtual UniquePtr<BaseCellLinkedList> createCellLinkedList(const BoundingBoxd &domain_bounds, BaseParticles &base_particles) override;
    virtual UniquePtr<LevelSet> createLevelSet(Shape &shape, Real refinement_ratio) override;

    class ContinuousSmoothingLengthRatio
    {
        Real *h_ratio_;

      public:
        explicit ContinuousSmoothingLengthRatio(AdaptiveSmoothingLength &adaptation)
            : h_ratio_(adaptation.h_ratio_) {};
        Real operator()(size_t index_i) { return h_ratio_[index_i]; };
    };

    typedef ContinuousSmoothingLengthRatio SmoothingRatioType;

  protected:
    Real finest_spacing_bound_;   /**< the adaptation bound for finest particles */
    Real coarsest_spacing_bound_; /**< the adaptation bound for coarsest particles */
};

/**
 * @class AdaptiveByShape
 * @brief Adaptive resolutions within a SPH body according to the distance to the body surface.
 */
class AdaptiveByShape : public AdaptiveSmoothingLength
{
  public:
    template <typename... Args>
    AdaptiveByShape(Args &&...args)
        : AdaptiveSmoothingLength(std::forward<Args>(args)...){};

    virtual ~AdaptiveByShape() {};
    virtual Real getLocalSpacing(Shape &shape, const Vecd &position) = 0;

  protected:
    Real smoothedSpacing(const Real &measure, const Real &transition_thickness);
};

/**
 * @class AdaptiveNearSurface
 * @brief Adaptive resolutions within a SPH body according to the distance to the body surface.
 */
class AdaptiveNearSurface : public AdaptiveByShape
{
  public:
    template <typename... Args>
    AdaptiveNearSurface(Args &&...args)
        : AdaptiveByShape(std::forward<Args>(args)...){};
    virtual ~AdaptiveNearSurface() {};

    virtual Real getLocalSpacing(Shape &shape, const Vecd &position) override;
};

/**
 * @class AdaptiveWithinShape
 * @brief Adaptive resolutions within a SPH body according to the distance to the body surface.
 */
class AdaptiveWithinShape : public AdaptiveByShape
{
  public:
    template <typename... Args>
    AdaptiveWithinShape(Args &&...args)
        : AdaptiveByShape(std::forward<Args>(args)...){};
    virtual ~AdaptiveWithinShape() {};

    virtual Real getLocalSpacing(Shape &shape, const Vecd &position) override;
};
} // namespace SPH
#endif // ADAPTATION_H