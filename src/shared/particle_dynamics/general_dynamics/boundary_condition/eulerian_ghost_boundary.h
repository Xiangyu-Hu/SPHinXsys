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
 * @file eulerian_ghost_boundary.h
 * @brief This is the particle dynamics for domain bounding
 * @author Zhentong Wang, Xiangyu Hu
 */
#ifndef Eulerian_GHOST_BOUNDARY_H
#define Eulerian_GHOST_BOUNDARY_H

#include "base_general_dynamics.h"
#include "particle_reserve.h"
#include "unstructured_mesh.h"
namespace SPH
{
struct RealAndGhostParticleData
{
    size_t real_index_, ghost_index_;
    Vecd e_ij_ghost_;
};

/**
 * @class GhostCreationInESPH
 * @brief The extra storage for the boundary particles, similar with the ghost particle,
 * is established to obey the zero-order consistency.
 */
class GhostCreationInESPH : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit GhostCreationInESPH(BaseInnerRelation &inner_relation, Ghost<ReserveSizeFactor> &ghost_boundary);
    virtual ~GhostCreationInESPH(){};
    std::vector<RealAndGhostParticleData> real_and_ghost_particle_data_;
    void ghostGenerationAndAddToConfiguration();

  protected:
    Ghost<ReserveSizeFactor> &ghost_boundary_;
    std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
    NeighborBuilderInnerInFVM get_inner_neighbor_;
    int *indicator_;
    Real *Vol_;
    Vecd *pos_;

  public:
    std::pair<size_t, size_t> &ghost_bound_;
};

//----------------------------------------------------------------------
//	GhostBoundaryConditionSetupInESPH
//----------------------------------------------------------------------
class GhostBoundaryConditionSetupInESPH : public LocalDynamics, public DataDelegateInner
{
  public:
    GhostBoundaryConditionSetupInESPH(BaseInnerRelation &inner_relation, GhostCreationInESPH &ghost_creation);
    virtual ~GhostBoundaryConditionSetupInESPH(){};
    virtual void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyGivenValueInletFlow(size_t ghost_index){};
    virtual void applyOutletBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyTopBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyFarFieldBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyPressureOutletBC(size_t ghost_index, size_t index_i){};
    virtual void applySymmetryBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyVelocityInletFlow(size_t ghost_index, size_t index_i){};
    virtual void setupBoundaryTypes(){};
    void resetBoundaryConditions();

  protected:
    Real *rho_, *Vol_, *mass_;
    Vecd *vel_, *pos_, *mom_;
    std::pair<size_t, size_t> &ghost_bound_;
    std::vector<RealAndGhostParticleData> real_and_ghost_particle_data_;
    int *boundary_type_;
    Real W0_;
};
/**
 * @class GhostKernelGradientUpdate
 * @brief Here, this class is used to update the boundary particle
          after implementing the kernel correction to strictly follow the zero-order consistency
 */
class GhostKernelGradientUpdate : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit GhostKernelGradientUpdate(BaseInnerRelation &inner_relation);
    virtual ~GhostKernelGradientUpdate(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_;
    Vecd *kernel_gradient_original_summation_;
    int *indicator_;
};
} // namespace SPH
#endif // Eulerian_GHOST_BOUNDARY_H
