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
 * @file anisotropic_adaptation.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef ANISOTROPIC_ADAPTATION_H
#define ANISOTROPIC_ADAPTATION_H

#include "base_adaptation.h"

#include "base_particles.hpp"
#include "cell_linked_list.h"
#include "level_set_shape.h"

namespace SPH
{

template <class IsotropicAdaptationType>
class Anisotropy : public IsotropicAdaptationType
{
  public:
    template <typename... Args>
    Anisotropy(const Vecd &scaling, const Vecd &orientation, Args &&...args)
        : IsotropicAdaptationType(std::forward<Args>(args)...),
          scaling_ref_(scaling), orientation_ref_(orientation),
          deformation_matrix_ref_(Matd::Identity()),
          spacing_ref_min_(this->spacing_ref_ * scaling_ref_.minCoeff()),
          h_ref_min_(this->h_ref_ * scaling_ref_.minCoeff()),
          h_ref_max_(this->h_ref_ * scaling_ref_.maxCoeff()),
          dv_scaling_(nullptr), dv_orientation_(nullptr),
          dv_deformation_matrix_(nullptr), dv_deformation_det_(nullptr)
    {
        deformation_matrix_ref_ =
            scaling_ref_.cwiseInverse().asDiagonal() * RotationMatrix(Vecd::UnitX(), orientation_ref_);
        this->spacing_min_ = spacing_ref_min_;
    }
    virtual ~Anisotropy() {};
    Real MinRefSmoothingLength() const { return h_ref_min_; }
    Real MinCutOffRadius() const { return kernel_ptr_->KernelSize() * h_ref_min_; }
    Real MaxCutOffRadius() const { return kernel_ptr_->KernelSize() * h_ref_max_; }

    virtual void initializeAdaptationVariables(BaseParticles &particles) override
    {
        dv_scaling_ = particles.registerStateVariable<Vecd>("AnisotropicScaling", scaling_ref_);
        dv_orientation_ = particles.registerStateVariable<Vecd>("AnisotropicOrientation", orientation_ref_);
        dv_deformation_matrix_ = particles.registerStateVariable<Matd>("AnisotropicMatrix", deformation_matrix_ref_);
        dv_deformation_det_ = particles.registerStateVariable<Real>("AnisotropicDeterminate", deformation_matrix_ref_.determinant());
        particles.addVariableToWrite<Vecd>(dv_scaling_);
        particles.addVariableToWrite<Vecd>(dv_orientation_);
    };

    typedef Anisotropy<IsotropicAdaptationType> CellLinkedListIdentifier;
    virtual UniquePtr<BaseCellLinkedList> createCellLinkedList(
        const BoundingBoxd &domain_bounds, BaseParticles &particles) override
    {
        return makeUnique<CellLinkedList<CellLinkedListIdentifier>>(
            domain_bounds, MinCutOffRadius(), particles, *this);
    };

    virtual UniquePtr<LevelSet> createLevelSet(Shape &shape, Real refinement) const override
    {
        // estimate the required minimum grid spacing for the level set
        int total_levels = (int)log10(shape.getBounds().MinimumDimension() / spacing_ref_min_) + 2;
        Real coarsest_spacing = spacing_ref_min_ * pow(2.0, total_levels - 1);
        LevelSet coarser_level_sets(shape.getBounds(), coarsest_spacing / refinement,
                                    total_levels - 1, shape, *this, refinement);
        // return the finest level set only
        return makeUnique<LevelSet>(shape.getBounds(), &coarser_level_sets, shape, *this, refinement);
    };

  protected:
    Vecd scaling_ref_;
    Vecd orientation_ref_;
    Matd deformation_matrix_ref_;
    Real spacing_ref_min_;
    Real h_ref_min_;
    Real h_ref_max_;

    DiscreteVariable<Vecd> *dv_scaling_, *dv_orientation_;
    DiscreteVariable<Matd> *dv_deformation_matrix_;
    DiscreteVariable<Real> *dv_deformation_det_;
};
} // namespace SPH
#endif // ANISOTROPIC_ADAPTATION_H