#include "level_set_correction.hpp"

namespace SPH
{
//=============================================================================================//
MarkCutInterfaces::MarkCutInterfaces(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level, Real perturbation_ratio)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      perturbation_ratio_(perturbation_ratio),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      pmv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
MarkNearInterface::MarkNearInterface(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      pmv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
ReinitializeLevelSet::ReinitializeLevelSet(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      pmv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
RedistanceInterface::RedistanceInterface(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      pmv_phi_gradient_(*data_mesh.getPackageVariable<Vecd>("LevelSetGradient")),
      pmv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
DiffuseLevelSetSign::DiffuseLevelSetSign(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
    SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      pmv_near_interface_id_(*data_mesh.getPackageVariable<int>("NearInterfaceID")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      sv_count_modified_(sv_count_modified) {}
//=============================================================================================//
LevelSetSignFromFine::LevelSetSignFromFine(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      pmv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      dv_pkg_1d_cell_index_(data_mesh.getPackage1DCellIndex()),
      fine_resolution_level_(resolution_level + 1)
{
    if (fine_resolution_level_ > data_mesh.ResolutionLevels() - 1)
    {
        std::cout << "\nError: in" << type_name<LevelSetSignFromFine>()
                  << ", the fine resolution level: " << fine_resolution_level_
                  << " is not exist !" << std::endl;
        exit(1);
    }
}
//=================================================================================================//
} // namespace SPH
