#ifndef MESH_WITH_DATA_PACKAGES_HXX
#define MESH_WITH_DATA_PACKAGES_HXX

#include "mesh_with_data_packages.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
PackageMesh<PKG_SIZE>::PackageMesh(
    BoundingBoxd tentative_bounds, Real grid_spacing,
    UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset)
    : Mesh(tentative_bounds, grid_spacing, buffer_width, linear_cell_index_offset),
      data_spacing_(grid_spacing / Real(PKG_SIZE)) {}
//=============================================================================================//
template <int PKG_SIZE>
Vecd PackageMesh<PKG_SIZE>::DataLowerBoundInCell(const Arrayi &cell_index) const
{
    return CellLowerCornerPosition(cell_index) + 0.5 * data_spacing_ * Vecd::Ones();
}
//=============================================================================================//
template <int PKG_SIZE>
Arrayi PackageMesh<PKG_SIZE>::DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position) const
{
    return floor((position - DataLowerBoundInCell(cell_index)).array() / data_spacing_)
        .template cast<int>()
        .max(Arrayi::Zero())
        .min((PKG_SIZE - 1) * Arrayi::Ones());
}
//=============================================================================================//
template <int PKG_SIZE>
Vecd PackageMesh<PKG_SIZE>::DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index) const
{
    return DataLowerBoundInCell(cell_index) + data_index.cast<Real>().matrix() * data_spacing_;
}
//=============================================================================================//
template <int PKG_SIZE>
UnsignedInt PackageMesh<PKG_SIZE>::PackageIndexFromCellIndex(
    UnsignedInt *cell_pkg_index, const Arrayi &cell_index) const
{
    UnsignedInt index_1d = LinearCellIndex(cell_index);
    return cell_pkg_index[index_1d];
}
//=============================================================================================//
template <int PKG_SIZE>
bool PackageMesh<PKG_SIZE>::isWithinCorePackage(
    UnsignedInt *cell_pkg_index, int *pkg_type, const Vecd &position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    UnsignedInt pkg_index = PackageIndexFromCellIndex(cell_pkg_index, cell_index);
    return pkg_type[pkg_index] == 1;
}
//=============================================================================================//
template <int PKG_SIZE>
Mesh PackageMesh<PKG_SIZE>::getGlobalMesh() const
{
    return Mesh(MeshLowerBound() + 0.5 * Vecd::Constant(data_spacing_),
                data_spacing_, AllCells() * PKG_SIZE);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType PackageMesh<PKG_SIZE>::DataValueFromGlobalIndex(
    PackageData<DataType, PKG_SIZE> *pkg_data, const Arrayi &global_grid_index,
    UnsignedInt *cell_package_index) const
{
    Arrayi cell_index_on_mesh = global_grid_index / PKG_SIZE;
    Arrayi local_index = global_grid_index - cell_index_on_mesh * PKG_SIZE;
    UnsignedInt package_index = PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh);
    return pkg_data[package_index](local_index);
}
//=============================================================================================//
template <int PKG_SIZE>
MeshWithGridDataPackages<PKG_SIZE>::MeshWithGridDataPackages(
    const std::string &name, UnsignedInt resolution_levels, BoundingBoxd tentative_bounds,
    Real reference_grid_spacing, UnsignedInt buffer_width, UnsignedInt num_singular_pkgs)
    : MultiResolutionMeshField<PackageMesh<PKG_SIZE>>(
          name, resolution_levels, tentative_bounds, reference_grid_spacing, buffer_width),
      num_boundary_pkgs_(num_singular_pkgs),
      dv_pkg_1d_cell_index_(nullptr), dv_pkg_type_(nullptr), cell_neighborhood_(nullptr),
      mcv_cell_pkg_index_(this->template registerMeshCellVariable<UnsignedInt>("CellPackageIndex"))
{
    for (UnsignedInt i = 0; i != this->resolution_levels_ * num_boundary_pkgs_; i++)
    {
        occupied_data_pkgs_.push_back(std::make_pair(0, -1)); // for data alignment
    }
    num_pkgs_offsets_.resize(this->resolution_levels_ + 1);
    num_pkgs_offsets_[0] = occupied_data_pkgs_.size();
    pkgs_bound_ = occupied_data_pkgs_.size();
}
//=============================================================================================//
template <int PKG_SIZE>
MeshWithGridDataPackages<PKG_SIZE>::MeshWithGridDataPackages(
    BoundingBoxd tentative_bounds, Real data_spacing, UnsignedInt buffer_size, UnsignedInt num_singular_pkgs)
    : MultiResolutionMeshField<PackageMesh<PKG_SIZE>>(
          "SparseStorageMesh", 1, tentative_bounds, data_spacing * Real(4), buffer_size),
      num_boundary_pkgs_(num_singular_pkgs),
      dv_pkg_1d_cell_index_(nullptr), dv_pkg_type_(nullptr), cell_neighborhood_(nullptr),
      mcv_cell_pkg_index_(this->template registerMeshCellVariable<UnsignedInt>("CellPackageIndex"))
{
    for (UnsignedInt i = 0; i != this->resolution_levels_ * num_boundary_pkgs_; i++)
    {
        occupied_data_pkgs_.push_back(std::make_pair(0, -1)); // for data alignment
    }
    num_pkgs_offsets_.resize(this->resolution_levels_ + 1);
    num_pkgs_offsets_[0] = occupied_data_pkgs_.size();
    pkgs_bound_ = occupied_data_pkgs_.size();
};
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<CellNeighborhood> &MeshWithGridDataPackages<PKG_SIZE>::getCellNeighborhood()
{
    return checkOrganized("getCellNeighborhood", *cell_neighborhood_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<UnsignedInt> &MeshWithGridDataPackages<PKG_SIZE>::getPackage1DCellIndex()
{
    return checkOrganized("getPackage1DCellIndex", *dv_pkg_1d_cell_index_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<int> &MeshWithGridDataPackages<PKG_SIZE>::getPackageType()
{
    return checkOrganized("getPackageType", *dv_pkg_type_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename T>
T &MeshWithGridDataPackages<PKG_SIZE>::checkOrganized(std::string func_name, T &value)
{
    if (!is_organized_)
    {
        std::cout << "\n Error: the mesh is not organized! (called from " << func_name << ")" << std::endl;
        exit(1);
    }
    return value;
}
//=============================================================================================//
template <int PKG_SIZE>
template <class ExecutionPolicy>
void MeshWithGridDataPackages<PKG_SIZE>::syncMeshVariablesToWrite(ExecutionPolicy &ex_policy)
{
    sync_mesh_variable_data_(mesh_variables_to_write_, ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <class ExecutionPolicy>
void MeshWithGridDataPackages<PKG_SIZE>::syncMeshVariablesToProbe(ExecutionPolicy &ex_policy)
{
    sync_mesh_variable_data_(mesh_variables_to_probe_, ex_policy);
    dv_pkg_1d_cell_index_->prepareForOutput(ex_policy);
    dv_pkg_type_->prepareForOutput(ex_policy);
    mcv_cell_pkg_index_->prepareForOutput(ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DiscreteVariable<PackageData<DataType, PKG_SIZE>> *MeshWithGridDataPackages<PKG_SIZE>::
    registerMeshVariable(const std::string &variable_name)
{
    if (!is_organized_)
    {
        std::cout << "\n Error: the mesh variable '" << variable_name
                  << "' is registered before the data packages are organized!" << std::endl;
        exit(1);
    }
    return registerVariable<MeshVariable, DataType>(
        all_mesh_variables_, mesh_variable_ptrs_, variable_name, pkgs_bound_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *MeshWithGridDataPackages<PKG_SIZE>::registerMetaVariable(
    const std::string &variable_name, Args &&...args)
{
    if (!is_organized_)
    {
        std::cout << "\n Error: the meta variable '" << variable_name
                  << "' is registered before the data packages are organized!" << std::endl;
        exit(1);
    }
    return registerVariable<DiscreteVariable, DataType>(
        all_meta_variables_, meta_variable_ptrs_, variable_name, pkgs_bound_, std::forward<Args>(args)...);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DiscreteVariable<PackageData<DataType, PKG_SIZE>> *MeshWithGridDataPackages<PKG_SIZE>::
    getMeshVariable(const std::string &variable_name)
{
    MeshVariable<DataType> *variable =
        findVariableByName<DataType, MeshVariable>(all_mesh_variables_, variable_name);
    if (variable == nullptr)
    {
        std::cout << "\n Error: the mesh variable '" << variable_name << "' is not exist!" << std::endl;
        exit(1);
    }
    return variable;
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DiscreteVariable<DataType> *MeshWithGridDataPackages<PKG_SIZE>::
    getMetaVariable(const std::string &variable_name)
{
    MetaVariable<DataType> *variable =
        findVariableByName<DataType, MetaVariable>(all_meta_variables_, variable_name);
    if (variable == nullptr)
    {
        std::cout << "\n Error: the meta variable '" << variable_name << "' is not exist!" << std::endl;
        exit(1);
    }
    return variable;
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addMeshVariableToWrite(const std::string &variable_name)
{
    addVariableToList<MeshVariable, DataType>(
        mesh_variables_to_write_, getMeshVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addMeshVariableToProbe(const std::string &variable_name)
{
    addVariableToList<MeshVariable, DataType>(
        mesh_variables_to_probe_, getMeshVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addEvolvingMetaVariable(const std::string &variable_name)
{
    addVariableToList<MetaVariable, DataType>(
        evolving_meta_variables_, getMetaVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::writeMeshFieldToPlt(
    const std::string &partial_file_name, size_t sequence)
{
    MultiResolutionMeshField<PackageMesh<PKG_SIZE>>::writeMeshFieldToPlt(partial_file_name, sequence);
    for (UnsignedInt l = 0; l != this->meshes_.size(); ++l)
    {
        std::string full_file_name = partial_file_name + "_MeshVariables_" + std::to_string(l) +
                                     std::to_string(sequence) + ".dat"; // level and sequence
        std::ofstream out_file(full_file_name.c_str(), std::ios::app);
        writeMeshVariablesToPltByMesh(l, out_file);
        out_file.close();
    }
}
//=============================================================================================//
template <int PKG_SIZE>
template <class DiscreteVariableType, class BoundaryDataFunction>
void MeshWithGridDataPackages<PKG_SIZE>::setBoundaryData(
    DiscreteVariableType *variable, UnsignedInt resolution_level,
    const BoundaryDataFunction &boundary_data_function)
{
    using ContainedDataType = typename DiscreteVariableType::ContainedDataType;
    for (UnsignedInt k = 0; k != num_boundary_pkgs_; k++)
        variable->setValue(
            resolution_level * num_boundary_pkgs_ + k,
            ContainedDataType(boundary_data_function(k)));
}
//=============================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::organizeOccupiedPackages()
{
    UnsignedInt coarsest_pkgs = occupied_data_pkgs_.size() - num_pkgs_offsets_[0];
    /* prodict the allocation size*/
    UnsignedInt scale = math::pow(2, Dimensions - 1);
    for (UnsignedInt i = 0; i != this->resolution_levels_; i++)
    {
        pkgs_bound_ += coarsest_pkgs * math::pow(scale, i);
    }
    pkgs_bound_ += pkgs_bound_ / 4; // add extra space for packages

    is_organized_ = true;

    dv_pkg_1d_cell_index_ = registerMetaVariable<UnsignedInt>(
        "Package1DCellIndex",
        [&](UnsignedInt i)
        { return i < occupied_data_pkgs_.size() ? occupied_data_pkgs_[i].first : 0; });
    dv_pkg_type_ = registerMetaVariable<int>(
        "PackageType",
        [&](UnsignedInt i)
        { return i < occupied_data_pkgs_.size() ? occupied_data_pkgs_[i].second : 0; });
    addEvolvingMetaVariable<UnsignedInt>("Package1DCellIndex");
    addEvolvingMetaVariable<int>("PackageType");
    cell_neighborhood_ = unique_variable_ptrs_.createPtr<
        MetaVariable<CellNeighborhood>>("CellNeighborhood", pkgs_bound_);
}
//=============================================================================================//
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_HXX