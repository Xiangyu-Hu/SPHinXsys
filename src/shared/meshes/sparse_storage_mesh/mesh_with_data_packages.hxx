#ifndef MESH_WITH_DATA_PACKAGES_HXX
#define MESH_WITH_DATA_PACKAGES_HXX

#include "mesh_with_data_packages.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
MeshWithGridDataPackages<PKG_SIZE>::MeshWithGridDataPackages(
    BoundingBoxd tentative_bounds, Real data_spacing, UnsignedInt buffer_size, UnsignedInt num_singular_pkgs)
    : index_handler_(tentative_bounds, data_spacing * PKG_SIZE, buffer_size, 0, data_spacing),
      num_singular_pkgs_(num_singular_pkgs), sv_num_grid_pkgs_("NumGridPackages", num_singular_pkgs),
      dv_pkg_1d_cell_index_(nullptr), dv_pkg_type_(nullptr), cell_neighborhood_(nullptr),
      bmv_cell_pkg_index_(registerBKGMeshVariable<UnsignedInt>("CellPackageIndex")),
      global_mesh_(index_handler_.MeshLowerBound() + 0.5 * data_spacing * Vecd::Ones(),
                   data_spacing, index_handler_.AllCells() * PKG_SIZE)
{
    for (UnsignedInt i = 0; i != num_singular_pkgs_; i++)
    {
        occupied_data_pkgs_.push_back(std::make_pair(0, 0)); // for data alignment
    }
};
//=============================================================================================//
template <int PKG_SIZE>
SingularVariable<UnsignedInt> &MeshWithGridDataPackages<PKG_SIZE>::svNumGridPackages()
{
    return checkOrganized("svNumGridPackages", sv_num_grid_pkgs_);
}
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
MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    IndexHandler(BoundingBoxd tentative_bounds, Real grid_spacing,
                 UnsignedInt buffer_width, UnsignedInt linear_cell_index_offset, Real data_spacing)
    : Mesh(tentative_bounds, grid_spacing, buffer_width, linear_cell_index_offset),
      data_spacing_(data_spacing) {}
//=============================================================================================//
template <int PKG_SIZE>
Vecd MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    DataLowerBoundInCell(const Arrayi &cell_index) const
{
    return CellLowerCornerPosition(cell_index) + 0.5 * data_spacing_ * Vecd::Ones();
}
//=============================================================================================//
template <int PKG_SIZE>
Arrayi MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    DataIndexFromPosition(const Arrayi &cell_index, const Vecd &position) const
{
    return floor((position - DataLowerBoundInCell(cell_index)).array() / data_spacing_)
        .template cast<int>()
        .max(Arrayi::Zero())
        .min((PKG_SIZE - 1) * Arrayi::Ones());
}
//=============================================================================================//
template <int PKG_SIZE>
Vecd MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    DataPositionFromIndex(const Arrayi &cell_index, const Arrayi &data_index) const
{
    return DataLowerBoundInCell(cell_index) + data_index.cast<Real>().matrix() * data_spacing_;
}
//=============================================================================================//
template <int PKG_SIZE>
UnsignedInt MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    PackageIndexFromCellIndex(UnsignedInt *cell_package_index, const Arrayi &cell_index) const
{
    UnsignedInt index_1d = LinearCellIndex(cell_index);
    return cell_package_index[index_1d];
}
//=============================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    isWithinCorePackage(UnsignedInt *cell_package_index, int *pkg_type, const Vecd &position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    UnsignedInt pkg_index = PackageIndexFromCellIndex(cell_package_index, cell_index);
    return pkg_type[pkg_index] == 1;
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
void MeshWithGridDataPackages<PKG_SIZE>::syncBKGMeshVariablesToWrite(ExecutionPolicy &ex_policy)
{
    sync_bkg_mesh_variable_data_(bkg_mesh_variables_to_write_, ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <class ExecutionPolicy>
void MeshWithGridDataPackages<PKG_SIZE>::syncMeshVariablesToProbe(ExecutionPolicy &ex_policy)
{
    sync_mesh_variable_data_(mesh_variables_to_probe_, ex_policy);
    dv_pkg_1d_cell_index_->prepareForOutput(ex_policy);
    dv_pkg_type_->prepareForOutput(ex_policy);
    bmv_cell_pkg_index_->prepareForOutput(ex_policy);
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
DiscreteVariable<DataType> *MeshWithGridDataPackages<PKG_SIZE>::registerBKGMeshVariable(
    const std::string &variable_name, Args &&...args)
{
    return registerVariable<BKGMeshVariable, DataType>(
        all_bkg_mesh_variables_, bkg_mesh_variable_ptrs_, variable_name,
        index_handler_.NumberOfCells(), std::forward<Args>(args)...);
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
    getBKGMeshVariable(const std::string &variable_name)
{
    BKGMeshVariable<DataType> *variable =
        findVariableByName<DataType, BKGMeshVariable>(all_bkg_mesh_variables_, variable_name);
    if (variable == nullptr)
    {
        std::cout << "\n Error: the BKG mesh variable '" << variable_name << "' is not exist!" << std::endl;
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
DataType MeshWithGridDataPackages<PKG_SIZE>::
    DataValueFromGlobalIndex(PackageData<DataType, PKG_SIZE> *pkg_data,
                             const Arrayi &global_grid_index, UnsignedInt *cell_package_index)
{
    Arrayi cell_index_on_mesh_ = global_grid_index / PKG_SIZE;
    Arrayi local_index = global_grid_index - cell_index_on_mesh_ * PKG_SIZE;
    UnsignedInt package_index = index_handler_.PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    return pkg_data[package_index](local_index);
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
void MeshWithGridDataPackages<PKG_SIZE>::addBKGMeshVariableToWrite(const std::string &variable_name)
{
    addVariableToList<DiscreteVariable, DataType>(
        bkg_mesh_variables_to_write_, getBKGMeshVariable<DataType>(variable_name));
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
void MeshWithGridDataPackages<PKG_SIZE>::organizeOccupiedPackages()
{
    sv_num_grid_pkgs_.setValue(occupied_data_pkgs_.size());
    pkgs_bound_ = sv_num_grid_pkgs_.getValue();
    is_organized_ = true;

    dv_pkg_1d_cell_index_ = registerMetaVariable<UnsignedInt>(
        "Package1DCellIndex",
        [&](UnsignedInt i)
        { return occupied_data_pkgs_[i].first; });
    dv_pkg_type_ = registerMetaVariable<int>(
        "PackageType",
        [&](UnsignedInt i)
        { return occupied_data_pkgs_[i].second; });
    addEvolvingMetaVariable<UnsignedInt>("Package1DCellIndex");
    addEvolvingMetaVariable<int>("PackageType");
    cell_neighborhood_ = unique_variable_ptrs_.createPtr<
        MetaVariable<CellNeighborhood>>("CellNeighborhood", pkgs_bound_);
}
//=============================================================================================//
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_HXX