#ifndef SPARSE_MESH_FIELD_HXX
#define SPARSE_MESH_FIELD_HXX

#include "spares_mesh_field.h"

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
bool PackageMesh<PKG_SIZE>::isWithinPackageType(
    int target_type, UnsignedInt *cell_pkg_index, int *pkg_type, const Vecd &position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    UnsignedInt pkg_index = PackageIndexFromCellIndex(cell_pkg_index, cell_index);
    return pkg_type[pkg_index] == target_type;
}
//=============================================================================================//
template <int PKG_SIZE>
Mesh PackageMesh<PKG_SIZE>::GlobalMesh() const
{
    return Mesh(MeshLowerBound() + 0.5 * Vecd::Constant(data_spacing_),
                data_spacing_, AllCells() * PKG_SIZE);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType PackageMesh<PKG_SIZE>::ValueByGlobalMesh(
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
SparseMeshField<PKG_SIZE>::SparseMeshField(
    const std::string &name, UnsignedInt resolution_levels, BoundingBoxd tentative_bounds,
    Real reference_grid_spacing, UnsignedInt buffer_width, UnsignedInt num_boundary_pkgs)
    : MultiResolutionMeshField<PackageMesh<PKG_SIZE>>(
          name, resolution_levels, tentative_bounds, reference_grid_spacing, buffer_width),
      num_boundary_pkgs_(num_boundary_pkgs),
      dv_pkg_1d_cell_index_(nullptr), dv_pkg_type_(nullptr), cell_neighborhood_(nullptr),
      mcv_cell_pkg_index_(this->template registerCellVariable<UnsignedInt>("CellPackageIndex"))
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
DiscreteVariable<CellNeighborhood> &SparseMeshField<PKG_SIZE>::getCellNeighborhood()
{
    return checkOrganized("getCellNeighborhood", *cell_neighborhood_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<UnsignedInt> &SparseMeshField<PKG_SIZE>::getPackage1DCellIndex()
{
    return checkOrganized("getPackage1DCellIndex", *dv_pkg_1d_cell_index_);
}
//=============================================================================================//
template <int PKG_SIZE>
DiscreteVariable<int> &SparseMeshField<PKG_SIZE>::getPackageType()
{
    return checkOrganized("getPackageType", *dv_pkg_type_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename T>
T &SparseMeshField<PKG_SIZE>::checkOrganized(std::string func_name, T &value)
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
void SparseMeshField<PKG_SIZE>::syncPackageVariablesToWrite(ExecutionPolicy &ex_policy)
{
    sync_mesh_variable_data_(pkg_variables_to_write_, ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <class ExecutionPolicy>
void SparseMeshField<PKG_SIZE>::syncPackageVariablesToProbe(ExecutionPolicy &ex_policy)
{
    sync_mesh_variable_data_(pkh_variables_to_probe_, ex_policy);
    dv_pkg_1d_cell_index_->prepareForOutput(ex_policy);
    dv_pkg_type_->prepareForOutput(ex_policy);
    mcv_cell_pkg_index_->prepareForOutput(ex_policy);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DiscreteVariable<PackageData<DataType, PKG_SIZE>> *SparseMeshField<PKG_SIZE>::
    registerPackageVariable(const std::string &variable_name)
{
    if (!is_organized_)
    {
        std::cout << "\n Error: the mesh variable '" << variable_name
                  << "' is registered before the data packages are organized!" << std::endl;
        exit(1);
    }
    return registerVariable<PackageVariable, DataType>(
        all_pkg_variables_, pkg_variable_ptrs_, variable_name, pkgs_bound_);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *SparseMeshField<PKG_SIZE>::registerMetaVariable(
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
DiscreteVariable<PackageData<DataType, PKG_SIZE>> *SparseMeshField<PKG_SIZE>::
    getPackageVariable(const std::string &variable_name)
{
    PackageVariable<DataType> *variable =
        findVariableByName<DataType, PackageVariable>(all_pkg_variables_, variable_name);
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
DiscreteVariable<DataType> *SparseMeshField<PKG_SIZE>::
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
void SparseMeshField<PKG_SIZE>::addPackageVariableToWrite(const std::string &variable_name)
{
    addVariableToList<PackageVariable, DataType>(
        pkg_variables_to_write_, getPackageVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void SparseMeshField<PKG_SIZE>::addPackageVariableToProbe(const std::string &variable_name)
{
    addVariableToList<PackageVariable, DataType>(
        pkh_variables_to_probe_, getPackageVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
void SparseMeshField<PKG_SIZE>::addEvolvingMetaVariable(const std::string &variable_name)
{
    addVariableToList<MetaVariable, DataType>(
        evolving_meta_variables_, getMetaVariable<DataType>(variable_name));
}
//=============================================================================================//
template <int PKG_SIZE>
void SparseMeshField<PKG_SIZE>::writeMeshFieldToPlt(
    const std::string &partial_file_name, size_t sequence)
{
    MultiResolutionMeshField<PackageMesh<PKG_SIZE>>::writeMeshFieldToPlt(partial_file_name, sequence);
    for (UnsignedInt l = 0; l != this->resolution_levels_; ++l)
    {
        std::string full_file_name = partial_file_name + "_PackageVariables_" + std::to_string(l) +
                                     std::to_string(sequence) + ".dat"; // level and sequence
        std::ofstream out_file(full_file_name.c_str(), std::ios::app);
        writePackageVariablesToPltByMesh(l, out_file);
        out_file.close();
    }
}
//=============================================================================================//
template <int PKG_SIZE>
template <class DiscreteVariableType, class BoundaryDataFunction>
void SparseMeshField<PKG_SIZE>::setBoundaryData(
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
void SparseMeshField<PKG_SIZE>::organizeOccupiedPackages(UnsignedInt level)
{
    num_pkgs_offsets_[level + 1] = occupied_data_pkgs_.size();
    if (level == 0)
    {
        pkgs_bound_ = num_pkgs_offsets_[1];
        for (UnsignedInt i = 1; i != this->resolution_levels_; i++)
        {
            pkgs_bound_ += pkgs_bound_ * pow(2, Dimensions - 1) + pkgs_bound_ / 4;
        }
        pkgs_bound_ += pkgs_bound_ / 4; // add extra for safety
        is_organized_ = true;

        dv_pkg_1d_cell_index_ = registerMetaVariable<UnsignedInt>("Package1DCellIndex");
        dv_pkg_type_ = registerMetaVariable<int>("PackageType");
        addEvolvingMetaVariable<UnsignedInt>("Package1DCellIndex");
        addEvolvingMetaVariable<int>("PackageType");
        cell_neighborhood_ = unique_variable_ptrs_.createPtr<
            MetaVariable<CellNeighborhood>>("CellNeighborhood", pkgs_bound_);
    }
    dv_pkg_1d_cell_index_->fill(num_pkgs_offsets_[level], num_pkgs_offsets_[level + 1],
                                [&](UnsignedInt i)
                                { return occupied_data_pkgs_[i].first; });
    dv_pkg_type_->fill(num_pkgs_offsets_[level], num_pkgs_offsets_[level + 1],
                       [&](UnsignedInt i)
                       { return occupied_data_pkgs_[i].second; });
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
template <class ExecutionPolicy>
SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::ProbeMesh(
    const ExecutionPolicy &ex_policy, SparseMeshField<PKG_SIZE> &encloser, const std::string &name)
    : pkg_data_(encloser.getPackageVariable<DataType>(name)->DelegatedData(ex_policy)),
      index_handler_(encloser.caMeshes().DelegatedData(ex_policy)),
      resolution_levels_(encloser.ResolutionLevels()),
      cell_pkg_index_(encloser.getCellPackageIndex().DelegatedData(ex_policy)),
      pkg_type_(encloser.getPackageType().DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.getCellNeighborhood().DelegatedData(ex_policy)) {}
//=================================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::operator()(
    IndexHandler &index_handler, const Vecd &position)
{
    Arrayi cell_index = index_handler.CellIndexFromPosition(position);
    UnsignedInt package_index =
        index_handler.PackageIndexFromCellIndex(cell_pkg_index_, cell_index);
    return pkg_type_[package_index] != -1
               ? probePackageData(index_handler, package_index, cell_index, position)
               : pkg_data_[package_index](Arrayi::Zero());
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::probeInResolutionLevel(
    UnsignedInt level, const Vecd &position)
{
    return operator()(index_handler_[level], position);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::probeBetweenResolutionLevels(
    UnsignedInt coarser_level, Real coarser_weight, const Vecd &position)
{
    return operator()(index_handler_[coarser_level], position) * coarser_weight +
           operator()(index_handler_[coarser_level + 1], position) * (1.0 - coarser_weight);
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
UnsignedInt SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::locateResolutionLevelByPackageType(
    int target_type, const Vecd &position)
{
    UnsignedInt proble_level = 0;
    for (size_t level = resolution_levels_ - 1; level != 0; --level)
    {
        if (index_handler_[level].isWithinPackageType(target_type, cell_pkg_index_, pkg_type_, position))
        {
            return level; // jump out of the loop!
        }
    }
    return proble_level;
}
//=============================================================================================//
} // namespace SPH
#endif // SPARSE_MESH_FIELD_HXX