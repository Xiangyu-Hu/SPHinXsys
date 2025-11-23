#ifndef MESH_WITH_DATA_PACKAGES_HXX
#define MESH_WITH_DATA_PACKAGES_HXX

#include "mesh_with_data_packages.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
MeshWithGridDataPackages<PKG_SIZE>::MeshWithGridDataPackages(
    BoundingBox tentative_bounds, Real data_spacing, UnsignedInt buffer_size, UnsignedInt num_singular_pkgs)
    : Mesh(tentative_bounds, PKG_SIZE * data_spacing, buffer_size),
      global_mesh_(mesh_lower_bound_ + 0.5 * data_spacing * Vecd::Ones(), data_spacing, all_cells_ * PKG_SIZE),
      num_singular_pkgs_(num_singular_pkgs), sv_num_grid_pkgs_("NumGridPackages", num_singular_pkgs),
      dv_pkg_cell_index_("PackageCellIndex", num_singular_pkgs_),
      dv_pkg_type_("PackageType", num_singular_pkgs_),
      cell_neighborhood_("CellNeighborhood", num_singular_pkgs_),
      bmv_cell_pkg_index_(*registerBKGMeshVariable<UnsignedInt>("CellPackageIndex")),
      data_spacing_(data_spacing), index_handler_(*this, data_spacing){};
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
    return checkOrganized("getCellNeighborhood", cell_neighborhood_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<UnsignedInt> &MeshWithGridDataPackages<PKG_SIZE>::getCellPackageIndex()
{
    return checkOrganized("getCellPackageIndex", bmv_cell_pkg_index_);
}
//=============================================================================================//
template <int PKG_SIZE>
ConcurrentVec<std::pair<UnsignedInt, int>> &MeshWithGridDataPackages<PKG_SIZE>::getOccupiedDataPackages()
{
    return checkOrganized("getOccupiedDataPackages", occupied_data_pkgs_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<Arrayi> &MeshWithGridDataPackages<PKG_SIZE>::getPackageCellIndex()
{
    return checkOrganized("getPackageCellIndex", dv_pkg_cell_index_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<int> &MeshWithGridDataPackages<PKG_SIZE>::getPackageType()
{
    return checkOrganized("getPackageType", dv_pkg_type_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addMeshVariableToWrite(const std::string &variable_name)
{
    addVariableToList<MeshVariable, DataType>(mesh_variables_to_write_, all_mesh_variables_, variable_name);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addMeshVariableToProbe(const std::string &variable_name)
{
    addVariableToList<MeshVariable, DataType>(mesh_variables_to_probe_, all_mesh_variables_, variable_name);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addBKGMeshVariableToWrite(const std::string &variable_name)
{
    addVariableToList<DiscreteVariable, DataType>(
        bkg_mesh_variables_to_write_, all_bkg_mesh_variables_, variable_name);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::addEvolvingMetaVariable(const std::string &variable_name)
{
    addVariableToList<DiscreteVariable, DataType>(
        evolving_meta_variables_, all_meta_variables_, variable_name);
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
template <template <typename> typename ContainerType, typename DataType>
void MeshWithGridDataPackages<PKG_SIZE>::
    addVariableToList(DataContainerAddressAssemble<ContainerType> &variable_set,
                      DataContainerAddressAssemble<ContainerType> &all_variable_set,
                      const std::string &variable_name)
{
    ContainerType<DataType> *variable =
        findVariableByName<DataType, ContainerType>(all_variable_set, variable_name);

    if (variable == nullptr)
    {
        std::cout << "\n Error: the mesh variable '" << variable_name << "' is  not exist!" << std::endl;
        exit(1);
    }

    ContainerType<DataType> *listed_variable =
        findVariableByName<DataType, ContainerType>(variable_set, variable_name);
    if (listed_variable == nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        std::get<type_index>(variable_set).push_back(variable);
    }
}
//=============================================================================================//
template <int PKG_SIZE>
template <template <typename> class MeshVariableType>
template <typename DataType, class ExecutionPolicy>
void MeshWithGridDataPackages<PKG_SIZE>::SyncMeshVariableData<MeshVariableType>::
operator()(DataContainerAddressKeeper<MeshVariableType<DataType>> &mesh_variables,
           ExecutionPolicy &ex_policy)
{
    for (UnsignedInt l = 0; l != mesh_variables.size(); l++)
    {
        MeshVariableType<DataType> *variable = mesh_variables[l];
        variable->prepareForOutput(ex_policy);
    }
}
//=============================================================================================//
template <int PKG_SIZE>
MeshWithGridDataPackages<PKG_SIZE>::IndexHandler::
    IndexHandler(const Mesh &mesh, Real data_spacing)
    : Mesh(mesh), data_spacing_(data_spacing) {}
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
    dv_pkg_cell_index_.prepareForOutput(ex_policy);
    dv_pkg_type_.prepareForOutput(ex_policy);
    bmv_cell_pkg_index_.prepareForOutput(ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <template <typename> typename ContainerType, typename DataType, typename... Args>
ContainerType<DataType> *MeshWithGridDataPackages<PKG_SIZE>::
    registerVariable(DataContainerAddressAssemble<ContainerType> &all_variable_set,
                     DataContainerUniquePtrAssemble<ContainerType> &all_variable_ptrs_,
                     const std::string &variable_name, Args &&...args)
{
    ContainerType<DataType> *variable =
        findVariableByName<DataType, ContainerType>(all_variable_set, variable_name);
    if (variable == nullptr)
    {
        return addVariableToAssemble<DataType, ContainerType>(
            all_variable_set, all_variable_ptrs_, variable_name, std::forward<Args>(args)...);
    }
    return variable;
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
typename MeshWithGridDataPackages<PKG_SIZE>::template MeshVariable<DataType> *
MeshWithGridDataPackages<PKG_SIZE>::registerMeshVariable(const std::string &variable_name)
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
typename MeshWithGridDataPackages<PKG_SIZE>::template BKGMeshVariable<DataType> *
MeshWithGridDataPackages<PKG_SIZE>::registerBKGMeshVariable(const std::string &variable_name, Args &&...args)
{
    return registerVariable<BKGMeshVariable, DataType>(
        all_bkg_mesh_variables_, bkg_mesh_variable_ptrs_, variable_name, AllCells().prod(),
        std::forward<Args>(args)...);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DiscreteVariable<DataType> *MeshWithGridDataPackages<PKG_SIZE>::registerMetaVariable(
    const std::string &variable_name)
{
    if (!is_organized_)
    {
        std::cout << "\n Error: the meta variable '" << variable_name
                  << "' is registered before the data packages are organized!" << std::endl;
        exit(1);
    }
    return registerVariable<DiscreteVariable, DataType>(
        all_meta_variables_, meta_variable_ptrs_, variable_name, pkgs_bound_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
typename MeshWithGridDataPackages<PKG_SIZE>::template MeshVariable<DataType> *
MeshWithGridDataPackages<PKG_SIZE>::getMeshVariable(const std::string &variable_name)
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
typename MeshWithGridDataPackages<PKG_SIZE>::template BKGMeshVariable<DataType> *
MeshWithGridDataPackages<PKG_SIZE>::getBKGMeshVariable(const std::string &variable_name)
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
DataType MeshWithGridDataPackages<PKG_SIZE>::
    DataValueFromGlobalIndex(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data,
                             const Arrayi &global_grid_index, UnsignedInt *cell_package_index)
{
    Arrayi cell_index_on_mesh_ = global_grid_index / PKG_SIZE;
    Arrayi local_index = global_grid_index - cell_index_on_mesh_ * PKG_SIZE;
    UnsignedInt package_index = index_handler_.PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    return pkg_data[package_index](local_index);
}
//=============================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::registerOccupied(UnsignedInt sort_index, int type)
{
    occupied_data_pkgs_.push_back(std::make_pair(sort_index, type));
}
//=============================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::organizeOccupiedPackages()
{
    sv_num_grid_pkgs_.setValue(occupied_data_pkgs_.size() + num_singular_pkgs_);
    pkgs_bound_ = sv_num_grid_pkgs_.getValue();
    cell_neighborhood_.reallocateData(par_host, pkgs_bound_);
    dv_pkg_cell_index_.reallocateData(par_host, pkgs_bound_);
    dv_pkg_type_.reallocateData(par_host, pkgs_bound_);
    is_organized_ = true;
}
//=============================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::isInnerDataPackage(const Arrayi &cell_index)
{
    UnsignedInt index_1d = LinearCellIndex(cell_index);
    /**
     * NOTE currently this func is only used in non-device mode;
     *      use the `DelegatedData` version when needed.
     */
    return bmv_cell_pkg_index_.Data()[index_1d] > 1;
}
//=============================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::isWithinCorePackage(
    UnsignedInt *cell_package_index, int *pkg_type, Vecd position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    UnsignedInt package_index = index_handler_.PackageIndexFromCellIndex(cell_package_index, cell_index);
    return pkg_type[package_index] == 1;
}
//=============================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::
    assignDataPackageIndex(const Arrayi &cell_index, const UnsignedInt package_index)
{
    UnsignedInt index_1d = LinearCellIndex(cell_index);
    /**
     * NOTE currently the `bmv_cell_pkg_index_` is only assigned in the host;
     *      use the `DelegatedData` version when needed.
     */
    bmv_cell_pkg_index_.Data()[index_1d] = package_index;
}
//=============================================================================================//
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_HXX