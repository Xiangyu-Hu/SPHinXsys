#include "level_set_correction.hpp"

namespace SPH
{
//=============================================================================================//
MarkCutInterfaces::MarkCutInterfaces(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level, Real perturbation_ratio)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      perturbation_ratio_(perturbation_ratio),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
MarkNearInterface::MarkNearInterface(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
ReinitializeLevelSet::ReinitializeLevelSet(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
RedistanceInterface::RedistanceInterface(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getPackageVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
DiffuseLevelSetSign::DiffuseLevelSetSign(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
    SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      sv_count_modified_(sv_count_modified) {}
//=================================================================================================//
} // namespace SPH
