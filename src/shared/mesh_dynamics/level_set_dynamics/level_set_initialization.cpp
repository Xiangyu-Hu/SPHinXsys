#include "level_set_initialization.hpp"

namespace SPH
{
//=================================================================================================//
InitialCellTagging::InitialCellTagging(MeshWithGridDataPackagesType &data_mesh, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh, 0), shape_(shape),
      occupied_data_pkgs_(data_mesh.getOccupiedDataPackages()),
      mcv_cell_pkg_index_(data_mesh.getCellPackageIndex()),
      mcv_cell_contain_id_(*data_mesh.registerMeshCellVariable<int>(
          "CellContainID", // default value is 2, indicating not near interface
          [&](UnsignedInt index)
          { return 2; })) // the value is not touched in this class
{
    data_mesh.addMeshCellVariableToWrite<int>("CellContainID");
}
//=================================================================================================//
void InitialCellTagging::UpdateKernel::update(const Arrayi &cell_index)
{
    UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
    cell_pkg_index_[index_1d] = 1;  // outside far field by default
    cell_contain_id_[index_1d] = 2; // default value is 2, indicating not near interface
    Vecd cell_position = index_handler_.CellPositionFromIndex(cell_index);
    if (ABS(shape_->findSignedDistance(cell_position)) < grid_spacing_)
    {
        cell_pkg_index_[index_1d] = 2;                               // indicate initially tagged temporarily
        occupied_data_pkgs_->push_back(std::make_pair(index_1d, 1)); // core package
    }
}
//=================================================================================================//
InitialCellTaggingFromCoarse::InitialCellTaggingFromCoarse(
    MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level,
    MeshWithGridDataPackagesType &coarse_mesh, UnsignedInt coarse_resolution_level, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      coarse_mesh_(coarse_mesh), coarse_resolution_level_(coarse_resolution_level), shape_(shape),
      occupied_data_pkgs_(data_mesh.getOccupiedDataPackages()),
      mcv_cell_pkg_index_(data_mesh.getCellPackageIndex()),
      mcv_cell_contain_id_(*data_mesh.registerMeshCellVariable<int>(
          "CellContainID", // default value is 2, indicating not near interface
          [&](UnsignedInt index)
          { return 2; })),
      mcv_cell_pkg_index_coarse_(coarse_mesh.getCellPackageIndex()),
      dv_pkg_type_coarse_(coarse_mesh.getPackageType())
{
    data_mesh.addMeshCellVariableToWrite<int>("CellContainID");
}
//=================================================================================================//
void InitialCellTaggingFromCoarse::UpdateKernel::update(const Arrayi &cell_index)
{
    Vecd cell_position = index_handler_.CellPositionFromIndex(cell_index);
    Real phi = probe_coarse_phi_(cell_position);
    UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
    cell_pkg_index_[index_1d] = phi < 0.0 ? 0 : 1;
    cell_contain_id_[index_1d] = 2;
    if (ABS(phi) > far_field_distance_)
    {
        cell_contain_id_[index_1d] = phi < 0.0 ? -1 : 1;
    }

    if (coarse_index_handler_.isWithinCorePackage(cell_pkg_index_coarse_, pkg_type_coarse_, cell_position))
    {
        if (ABS(shape_->findSignedDistance(cell_position)) < grid_spacing_)
        {
            cell_pkg_index_[index_1d] = 2;                               // indicate initially tagged temporarily
            occupied_data_pkgs_->push_back(std::make_pair(index_1d, 1)); // core package
        }
    }
}
//=================================================================================================//
InnerCellTagging::InnerCellTagging(MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      occupied_data_pkgs_(data_mesh.getOccupiedDataPackages()),
      mcv_cell_pkg_index_(data_mesh.getCellPackageIndex()),
      num_pkgs_offsets_(data_mesh.getNumPackageOffsets()) {}
//=================================================================================================//
void InnerCellTagging::finishDynamics(Real dt)
{
    num_pkgs_offsets_[this->resolution_level_ + 1] = occupied_data_pkgs_.size();
}
//=================================================================================================//
void InnerCellTagging::UpdateKernel::update(const Arrayi &cell_index)
{
    if (isNearInitiallyTagged(cell_index) && !isInitiallyTagged(cell_index))
    {
        UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
        occupied_data_pkgs_->push_back(std::make_pair(index_1d, 0)); // inner package
    }
}
//=============================================================================================//
bool InnerCellTagging::UpdateKernel::isInitiallyTagged(const Arrayi &cell_index)
{
    UnsignedInt index_1d = index_handler_.LinearCellIndex(cell_index);
    return cell_pkg_index_[index_1d] == 2;
}
//=============================================================================================//
bool InnerCellTagging::UpdateKernel::isNearInitiallyTagged(const Arrayi &cell_index)
{
    return mesh_any_of(
        Arrayi::Zero().max(cell_index - Arrayi::Ones()),
        index_handler_.AllCells().min(cell_index + 2 * Arrayi::Ones()),
        [&](const Arrayi &index)
        {
            return isInitiallyTagged(index);
        });
}
//=================================================================================================//
InitializeCellNeighborhood::InitializeCellNeighborhood(
    MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      dv_pkg_1d_cell_index_(data_mesh.getPackage1DCellIndex()),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      mcv_cell_pkg_index_(data_mesh.getCellPackageIndex())
{
    data_mesh.setSingularPackages(
        &dv_cell_neighborhood_, resolution_level,
        [&](UnsignedInt l, UnsignedInt k)
        {
            return CellNeighborhood::Constant(data_mesh.NumSingularPackages() * l + k);
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::UpdateKernel::update(const UnsignedInt &package_index)
{
    CellNeighborhood &current = cell_neighborhood_[package_index];
    Arrayi cell_index = index_handler_.DimensionalCellIndex(pkg_1d_cell_index_[package_index]);
    mesh_for_each(
        -Arrayi::Ones(), Arrayi::Ones() * 2,
        [&](const Arrayi &index)
        {
            current(index + Arrayi::Ones()) =
                index_handler_.PackageIndexFromCellIndex(cell_pkg_index_, cell_index + index);
        });
}
//=============================================================================================//
InitializeBasicPackageData::InitializeBasicPackageData(
    MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level, Shape &shape)
    : BaseMeshLocalDynamics(data_mesh, resolution_level), shape_(shape),
      dv_pkg_1d_cell_index_(data_mesh.getPackage1DCellIndex()),
      mv_phi_(*data_mesh.registerMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      mv_near_interface_id_(*data_mesh.registerMeshVariable<int>("NearInterfaceID"))
{
    Real far_field_distance = index_handler_.GridSpacing() * (Real)index_handler_.BufferWidth();
    data_mesh.addMeshVariableToWrite<Real>("LevelSet");
    initializeSingularPackages(0, -far_field_distance);
    initializeSingularPackages(1, far_field_distance);
}
//=============================================================================================//
void InitializeBasicPackageData::initializeSingularPackages(
    const UnsignedInt package_index, Real far_field_level_set)
{
    auto &phi = mv_phi_.Data()[package_index];
    auto &near_interface_id = mv_near_interface_id_.Data()[package_index];
    auto &phi_gradient = mv_phi_gradient_.Data()[package_index];

    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &data_index)
                  {
                      phi(data_index) = far_field_level_set;
                      near_interface_id(data_index) = far_field_level_set < 0.0 ? -2 : 2;
                      phi_gradient(data_index) = Vecd::Ones();
                  });
}
//=============================================================================================//
void InitializeBasicPackageData::UpdateKernel::update(const UnsignedInt &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = index_handler_.DimensionalCellIndex(pkg_1d_cell_index_[package_index]);
    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &data_index)
                  {
                      Vecd position = index_handler_.DataPositionFromIndex(cell_index, data_index);
                      phi(data_index) = shape_->findSignedDistance(position);
                      near_interface_id(data_index) = phi(data_index) < 0.0 ? -2 : 2;
                  });
}
//=============================================================================================//
NearInterfaceCellTagging::NearInterfaceCellTagging(
    MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      dv_pkg_1d_cell_index_(data_mesh.getPackage1DCellIndex()),
      mcv_cell_contain_id_(*data_mesh.getMeshCellVariable<int>("CellContainID")),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")) {}
//=============================================================================================//
CellContainDiffusion::CellContainDiffusion(
    MeshWithGridDataPackagesType &data_mesh, UnsignedInt resolution_level,
    SingularVariable<UnsignedInt> &sv_count_modified)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mcv_cell_contain_id_(*data_mesh.getMeshCellVariable<int>("CellContainID")),
      mcv_cell_package_index_(data_mesh.getCellPackageIndex()),
      sv_count_modified_(sv_count_modified) {}
//=============================================================================================//
FinishDataPackages::FinishDataPackages(
    MeshWithGridDataPackagesType &mesh_data, UnsignedInt resolution_level, Shape &shape)
    : BaseDynamics<void>(), mesh_data_(mesh_data), resolution_level_(resolution_level), shape_(shape),
      initialize_cell_neighborhood{mesh_data_, resolution_level_},
      initialize_basic_data_for_a_package{mesh_data_, resolution_level_, shape_},
      near_interface_cell_tagging{mesh_data_, resolution_level_},
      cell_contain_diffusion{mesh_data_, resolution_level_, sv_count_modified_} {}
//=============================================================================================//
void FinishDataPackages::exec(Real dt)
{
    initialize_basic_data_for_a_package.exec();

    near_interface_cell_tagging.exec();
    while (sv_count_modified_.getValue() > 0)
    {
        sv_count_modified_.setValue(0);
        cell_contain_diffusion.exec();
    }

    initialize_cell_neighborhood.exec();
}
//=================================================================================================//
} // namespace SPH
