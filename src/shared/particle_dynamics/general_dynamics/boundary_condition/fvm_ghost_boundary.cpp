#include "fvm_ghost_boundary.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
GhostCreationFromMesh::GhostCreationFromMesh(RealBody &real_body, ANSYSMesh &ansys_mesh,
                                             Ghost<ReserveSizeFactor> &ghost_boundary)
    : GeneralDataDelegateSimple(real_body),
      ghost_boundary_(ghost_boundary),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      pos_(*particles_->getVariableByName<Vecd>("Position")),
      Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
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
    : GeneralDataDelegateInner(inner_relation),
      rho_(*particles_->getVariableByName<Real>("Density")),
      Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableByName<Real>("Mass")),
      p_(*particles_->getVariableByName<Real>("Pressure")),
      vel_(*particles_->getVariableByName<Vecd>("Velocity")),
      pos_(*particles_->getVariableByName<Vecd>("Position")),
      mom_(*particles_->getVariableByName<Vecd>("Momentum")),
      ghost_bound_(ghost_creation.ghost_bound_),
      each_boundary_type_with_all_ghosts_index_(ghost_creation.each_boundary_type_with_all_ghosts_index_),
      each_boundary_type_with_all_ghosts_eij_(ghost_creation.each_boundary_type_with_all_ghosts_eij_),
      each_boundary_type_contact_real_index_(ghost_creation.each_boundary_type_contact_real_index_) {}
//=================================================================================================//
} // namespace SPH
