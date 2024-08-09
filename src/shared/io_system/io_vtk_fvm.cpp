#include "io_vtk_fvm.h"

namespace SPH
{
//=================================================================================================//
BodyStatesRecordingInMeshToVtp::BodyStatesRecordingInMeshToVtp(SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecordingToVtp(body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_) {}
//=================================================================================================//
BodyStatesRecordingInMeshToVtu::BodyStatesRecordingInMeshToVtu(SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecordingToVtp(body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_), bounds_(body){};
//=================================================================================================//
} // namespace SPH
