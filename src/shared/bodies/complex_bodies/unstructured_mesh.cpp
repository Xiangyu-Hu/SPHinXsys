
#include "unstructured_mesh.h"

#include "base_particle_dynamics.h"

namespace SPH
{
//=================================================================================================//
BaseInnerRelationInFVM::BaseInnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh)
    : BaseInnerRelation(real_body), real_body_(&real_body),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      Vol_(*base_particles_.getVariableDataByName<Real>("VolumetricMeasure"))
{
    subscribeToBody();
    inner_configuration_.resize(base_particles_.RealParticlesBound(), Neighborhood());
};
//=============================================================================================//
} // namespace SPH
