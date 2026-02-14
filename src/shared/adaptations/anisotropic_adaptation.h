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

namespace SPH
{

class AnisotropicAdaptation : public SPHAdaptation
{
  public:
    AnisotropicAdaptation(const Vecd &scaling, const Vecd &orientation, Real global_resolution,
                          Real h_spacing_ratio = 1.3, Real refinement_to_global = 1.0);
    virtual ~AnisotropicAdaptation() {};
    Real MinRefSmoothingLength() const { return h_ref_min_; }
    Real MinCutOffRadius() const { return kernel_ptr_->KernelSize() * h_ref_min_; }
    Real MaxCutOffRadius() const { return kernel_ptr_->KernelSize() * h_ref_max_; }
    virtual void initializeAdaptationVariables(BaseParticles &particles) override;

    typedef AnisotropicAdaptation CellLinkedListIdentifier;
    virtual UniquePtr<BaseCellLinkedList> createCellLinkedList(
        const BoundingBoxd &domain_bounds, BaseParticles &particles) override;
    virtual UniquePtr<LevelSet> createLevelSet(Shape &shape, Real refinement) const override;

  protected:
    Vecd scaling_ref_;
    Vecd orientation_ref_;
    Matd deformation_matrix_ref_;
    Real spacing_ref_min_;
    Real h_ref_min_;
    Real h_ref_max_;

    DiscreteVariable<Vecd> *dv_scaling_, *dv_orientation_;
};
} // namespace SPH
#endif // ANISOTROPIC_ADAPTATION_H