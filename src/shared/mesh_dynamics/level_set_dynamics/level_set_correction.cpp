#include "level_set_correction.hpp"

namespace SPH
{
//=============================================================================================//
MarkCutInterfaces::MarkCutInterfaces(MeshWithGridDataPackagesType &data_mesh, Real perturbation_ratio)
    : BaseMeshLocalDynamics(data_mesh),
      perturbation_ratio_(perturbation_ratio),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
MarkNearInterface::MarkNearInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
ReinitializeLevelSet::ReinitializeLevelSet(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
RedistanceInterface::RedistanceInterface(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
DiffuseLevelSetSign::DiffuseLevelSetSign(
    MeshWithGridDataPackagesType &data_mesh, SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getMeshVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      sv_count_modified_(sv_count_modified) {}
//=================================================================================================//
} // namespace SPH
