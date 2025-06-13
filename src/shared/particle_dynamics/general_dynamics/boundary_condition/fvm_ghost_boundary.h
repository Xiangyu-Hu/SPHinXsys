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
 * @file fvm_ghost_boundary.h
 * @brief This is the particle dynamics for domain bounding
 * @author Xiangyu Hu
 */

#ifndef FVM_GHOST_BOUNDARY_H
#define FVM_GHOST_BOUNDARY_H

#include "base_general_dynamics.h"
#include "particle_reserve.h"
#include "unstructured_mesh.h"

namespace SPH
{
/**
 * @class BaseGhostCreation
 * @brief Base class for the ghost particle
 */
class GhostCreationFromMesh : public LocalDynamics
{
  public:
    GhostCreationFromMesh(RealBody &real_body, ANSYSMesh &ansys_mesh,
                          Ghost<ReserveSizeFactor> &ghost_boundary);
    virtual ~GhostCreationFromMesh() {};

  protected:
    Ghost<ReserveSizeFactor> &ghost_boundary_;
    std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
    StdVec<Vecd> &node_coordinates_;
    StdVec<StdVec<StdVec<size_t>>> &mesh_topology_;
    Vecd *pos_;
    Real *Vol_;
    void addGhostParticleAndSetInConfiguration();

  public:
    std::pair<size_t, size_t> &ghost_bound_;
    StdVec<StdVec<size_t>> each_boundary_type_with_all_ghosts_index_;
    StdVec<StdVec<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    StdVec<StdVec<size_t>> each_boundary_type_contact_real_index_;
};

//----------------------------------------------------------------------
//	BoundaryConditionSetupInFVM
//----------------------------------------------------------------------
class BoundaryConditionSetupInFVM : public LocalDynamics, public DataDelegateInner
{
  public:
    BoundaryConditionSetupInFVM(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation);
    virtual ~BoundaryConditionSetupInFVM() {};
    virtual void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) {};
    virtual void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) {};
    virtual void applyGivenValueInletFlow(size_t ghost_index) {};
    virtual void applyOutletBoundary(size_t ghost_index, size_t index_i) {};
    virtual void applyTopBoundary(size_t ghost_index, size_t index_i) {};
    virtual void applyFarFieldBoundary(size_t ghost_index) {};
    virtual void applyPressureOutletBC(size_t ghost_index, size_t index_i) {};
    virtual void applySymmetryBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) {};
    virtual void applyVelocityInletFlow(size_t ghost_index, size_t index_i) {};
    // Common functionality for resetting boundary conditions
    void resetBoundaryConditions();

  protected:
    Real *rho_, *Vol_, *mass_, *p_;
    Vecd *vel_, *pos_, *mom_;
    std::pair<size_t, size_t> &ghost_bound_;
    StdVec<StdVec<size_t>> &each_boundary_type_with_all_ghosts_index_;
    StdVec<StdVec<Vecd>> &each_boundary_type_with_all_ghosts_eij_;
    StdVec<StdVec<size_t>> &each_boundary_type_contact_real_index_;
};
} // namespace SPH
#endif // FVM_GHOST_BOUNDARY_H
