#include "fvm_ghost_boundary.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
GhostCreationFromMesh::GhostCreationFromMesh(RealBody &real_body, ANSYSMesh &ansys_mesh,
                                             Ghost<ReserveSizeFactor> &ghost_boundary)
    : LocalDynamics(real_body),
      ghost_boundary_(ghost_boundary),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      ghost_bound_(ghost_boundary.GhostBound())
{
    ghost_boundary.checkParticlesReserved();
    each_boundary_type_with_all_ghosts_index_.resize(50);
    each_boundary_type_with_all_ghosts_eij_.resize(50);
    each_boundary_type_contact_real_index_.resize(50);
    addGhostParticleAndSetInConfiguration();
}
//=================================================================================================//
BoundaryConditionSetupInFVM::
    BoundaryConditionSetupInFVM(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      mom_(particles_->getVariableDataByName<Vecd>("Momentum")),
      ghost_bound_(ghost_creation.ghost_bound_),
      each_boundary_type_with_all_ghosts_index_(ghost_creation.each_boundary_type_with_all_ghosts_index_),
      each_boundary_type_with_all_ghosts_eij_(ghost_creation.each_boundary_type_with_all_ghosts_eij_),
      each_boundary_type_contact_real_index_(ghost_creation.each_boundary_type_contact_real_index_) {}
//=================================================================================================//
} // namespace SPH
