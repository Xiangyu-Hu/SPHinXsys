#include "cell_linked_list.h"
#include "adaptation.h"
#include "base_kernel.h"
#include "base_particles.h"
#include "mesh_iterators.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
BaseCellLinkedList::BaseCellLinkedList(BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : BaseMeshField("CellLinkedList"), kernel_(*sph_adaptation.getKernel()),
      number_of_split_cell_lists_(static_cast<UnsignedInt>(pow(3, Dimensions))),
      dv_particle_index_(nullptr), dv_cell_offset_(nullptr),
      cell_index_lists_(nullptr), cell_data_lists_(nullptr) {}
//=================================================================================================//
BaseCellLinkedList::~BaseCellLinkedList()
{
    delete[] cell_index_lists_;
    delete[] cell_data_lists_;
}
//=================================================================================================//
void BaseCellLinkedList::initialize(BaseParticles &base_particles)
{
    cell_offset_list_size_ = total_number_of_cells_ + 1;
    index_list_size_ = SMAX(base_particles.ParticlesBound(), cell_offset_list_size_);
    dv_particle_index_ = unique_variable_ptrs_.createPtr<DiscreteVariable<UnsignedInt>>("ParticleIndex", index_list_size_);
    dv_cell_offset_ = unique_variable_ptrs_.createPtr<DiscreteVariable<UnsignedInt>>("CellOffset", cell_offset_list_size_);
    cell_index_lists_ = new ConcurrentIndexVector[total_number_of_cells_];
    cell_data_lists_ = new ListDataVector[total_number_of_cells_];
}
//=================================================================================================//
void BaseCellLinkedList::clearCellLists()
{
    parallel_for(
        IndexRange(0, total_number_of_cells_),
        [&](const IndexRange &r)
        {
            for (UnsignedInt i = r.begin(); i != r.end(); ++i)
            {
                cell_index_lists_[i].clear();
            }
        },
        ap);
}
//=================================================================================================//
void BaseCellLinkedList::UpdateCellListData(BaseParticles &base_particles)
{
    Vecd *pos = base_particles.ParticlePositions();
    parallel_for(
        IndexRange(0, total_number_of_cells_),
        [&](const IndexRange &r)
        {
            for (UnsignedInt i = r.begin(); i != r.end(); ++i)
            {
                ListDataVector &cell_data_list = cell_data_lists_[i];
                cell_data_list.clear();
                ConcurrentIndexVector &cell_list = cell_index_lists_[i];
                for (UnsignedInt s = 0; s != cell_list.size(); ++s)
                {
                    UnsignedInt index = cell_list[s];
                    cell_data_list.emplace_back(std::make_pair(index, pos[index]));
                }
            }
        },
        ap);
}
//=================================================================================================//
void BaseCellLinkedList::tagBodyPartByCellByMesh(Mesh &mesh, UnsignedInt mesh_offset,
                                                 ConcurrentCellLists &cell_lists,
                                                 ConcurrentIndexVector &cell_indexes,
                                                 std::function<bool(Vecd, Real)> &check_included)
{
    mesh_parallel_for(
        MeshRange(Arrayi::Zero(), mesh.AllCells()),
        [&](const Arrayi &cell_index)
        {
            bool is_included = false;
            mesh_for_each(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                mesh.AllCells().min(cell_index + 2 * Arrayi::Ones()),
                [&](const Arrayi &neighbor_cell_index)
                {
                    if (check_included(mesh.CellPositionFromIndex(neighbor_cell_index), mesh.GridSpacing()))
                    {
                        is_included = true;
                    }
                });
            if (is_included == true)
            {
                UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                cell_lists.push_back(&cell_index_lists_[linear_index]);
                cell_indexes.push_back(linear_index);
            }
        });
}
//=================================================================================================//
void BaseCellLinkedList::UpdateCellLists(BaseParticles &base_particles)
{
    clearCellLists();
    Vecd *pos_n = base_particles.ParticlePositions();
    UnsignedInt total_real_particles = base_particles.TotalRealParticles();
    parallel_for(
        IndexRange(0, total_real_particles),
        [&](const IndexRange &r)
        {
            for (UnsignedInt i = r.begin(); i != r.end(); ++i)
            {
                insertParticleIndex(i, pos_n[i]);
            }
        },
        ap);

    UpdateCellListData(base_particles);
}
//=================================================================================================//
void BaseCellLinkedList::findNearestListDataEntryByMesh(Mesh &mesh, UnsignedInt mesh_offset,
                                                        Real &min_distance_sqr, ListData &nearest_entry,
                                                        const Vecd &position)
{
    Arrayi cell = mesh.CellIndexFromPosition(position);
    mesh_for_each(
        Arrayi::Zero().max(cell - Arrayi::Ones()),
        mesh.AllCells().min(cell + 2 * Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
            ListDataVector &target_particles = cell_data_lists_[linear_index];
            for (const ListData &list_data : target_particles)
            {
                Real distance_sqr = (position - std::get<1>(list_data)).squaredNorm();
                if (distance_sqr < min_distance_sqr)
                {
                    min_distance_sqr = distance_sqr;
                    nearest_entry = list_data;
                }
            }
        });
}
//=================================================================================================//
CellLinkedList::CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
                               BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : BaseCellLinkedList(base_particles, sph_adaptation), mesh_(nullptr)
{
    mesh_ = mesh_ptrs_keeper_.createPtr<Mesh>(tentative_bounds, grid_spacing, 2);
    meshes_.push_back(mesh_);
    mesh_offsets_.push_back(0);
    total_number_of_cells_ = mesh_->NumberOfCells();
    initialize(base_particles);
}
//=================================================================================================//
void CellLinkedList ::insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position)
{
    UnsignedInt linear_index = mesh_->LinearCellIndexFromPosition(particle_position);
    cell_index_lists_[linear_index].emplace_back(particle_index);
}
//=================================================================================================//
void CellLinkedList ::InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position)
{
    UnsignedInt linear_index = mesh_->LinearCellIndexFromPosition(particle_position);
    cell_data_lists_[linear_index].emplace_back(std::make_pair(particle_index, particle_position));
}
//=================================================================================================//
void CellLinkedList::tagBoundingCells(StdVec<CellLists> &cell_data_lists,
                                      const BoundingBox &bounding_bounds, int axis)
{
    tagBoundingCellsByMesh(*mesh_, 0, cell_data_lists, bounding_bounds, axis);
}
//=================================================================================================//
ListData CellLinkedList::findNearestListDataEntry(const Vecd &position)
{
    Real min_distance_sqr = MaxReal;
    ListData nearest_entry = std::make_pair(MaxSize_t, MaxReal * Vecd::Ones());
    findNearestListDataEntryByMesh(*mesh_, 0, min_distance_sqr, nearest_entry, position);
    return nearest_entry;
}
//=================================================================================================//
UnsignedInt CellLinkedList::computingSequence(Vecd &position, UnsignedInt index_i)
{
    return mesh_->transferMeshIndexToMortonOrder(mesh_->CellIndexFromPosition(position));
}
//=================================================================================================//
void CellLinkedList::tagBodyPartByCell(ConcurrentCellLists &cell_lists,
                                       ConcurrentIndexVector &cell_indexes,
                                       std::function<bool(Vecd, Real)> &check_included)
{
    tagBodyPartByCellByMesh(*mesh_, 0, cell_lists, cell_indexes, check_included);
}
//=================================================================================================//
void CellLinkedList::writeMeshFieldToPlt(std::ofstream &output_file)
{
    writeMeshFieldToPltByMesh(*mesh_, 0, output_file);
}
//=================================================================================================//
MultilevelCellLinkedList::MultilevelCellLinkedList(
    BoundingBox tentative_bounds, Real reference_grid_spacing, UnsignedInt total_levels,
    BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : BaseCellLinkedList(base_particles, sph_adaptation),
      h_ratio_(DynamicCast<ParticleWithLocalRefinement>(this, &sph_adaptation)->h_ratio_),
      level_(DynamicCast<ParticleWithLocalRefinement>(this, &sph_adaptation)->level_)
{
    meshes_.push_back(mesh_ptrs_keeper_.createPtr<Mesh>(tentative_bounds, reference_grid_spacing, 2));
    mesh_offsets_.push_back(0);
    total_number_of_cells_ = meshes_[0]->NumberOfCells();
    for (UnsignedInt level = 1; level != total_levels; ++level)
    {
        /** all mesh levels aligned at the lower bound of tentative_bounds */
        Real refined_spacing = meshes_[level - 1]->GridSpacing() / 2.0;
        meshes_.push_back(mesh_ptrs_keeper_.createPtr<Mesh>(tentative_bounds, refined_spacing, 2));
        mesh_offsets_.push_back(total_number_of_cells_);
        total_number_of_cells_ += meshes_[level]->NumberOfCells();
    }
    initialize(base_particles);
}
//=================================================================================================//
UnsignedInt MultilevelCellLinkedList::getMeshLevel(Real particle_cutoff_radius)
{
    for (UnsignedInt level = meshes_.size(); level != 0; --level)
    {
        if (particle_cutoff_radius - meshes_[level - 1]->GridSpacing() < SqrtEps)
            return level - 1;
    }

    std::cout << "\n Error: CellLinkedList level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void MultilevelCellLinkedList::insertParticleIndex(UnsignedInt particle_index, const Vecd &particle_position)
{
    UnsignedInt level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    level_[particle_index] = level;
    UnsignedInt linear_index = mesh_offsets_[level] + meshes_[level]->LinearCellIndexFromPosition(particle_position);
    cell_index_lists_[linear_index].emplace_back(particle_index);
}
//=================================================================================================//
void MultilevelCellLinkedList::InsertListDataEntry(UnsignedInt particle_index, const Vecd &particle_position)
{
    UnsignedInt level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    UnsignedInt linear_index = mesh_offsets_[level] + meshes_[level]->LinearCellIndexFromPosition(particle_position);
    cell_data_lists_[linear_index]
        .emplace_back(std::make_pair(particle_index, particle_position));
}
//=================================================================================================//
UnsignedInt MultilevelCellLinkedList::computingSequence(Vecd &position, UnsignedInt index_i)
{
    UnsignedInt level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[index_i]));
    return meshes_[level]->transferMeshIndexToMortonOrder(
        meshes_[level]->CellIndexFromPosition(position));
}
//=================================================================================================//
void MultilevelCellLinkedList::tagBodyPartByCell(ConcurrentCellLists &cell_lists,
                                                 ConcurrentIndexVector &cell_indexes,
                                                 std::function<bool(Vecd, Real)> &check_included)
{
    for (UnsignedInt l = 0; l != meshes_.size(); ++l)
    {
        tagBodyPartByCellByMesh(*meshes_[l], mesh_offsets_[l], cell_lists, cell_indexes, check_included);
    }
}
//=================================================================================================//
void MultilevelCellLinkedList::tagBoundingCells(StdVec<CellLists> &cell_data_lists,
                                                const BoundingBox &bounding_bounds, int axis)
{
    for (UnsignedInt l = 0; l != meshes_.size(); ++l)
    {
        tagBoundingCellsByMesh(*meshes_[l], mesh_offsets_[l], cell_data_lists, bounding_bounds, axis);
    }
}
//=================================================================================================//
void MultilevelCellLinkedList::writeMeshFieldToPlt(std::ofstream &output_file)
{
    for (UnsignedInt l = 0; l != meshes_.size(); ++l)
    {
        writeMeshFieldToPltByMesh(*meshes_[l], mesh_offsets_[l], output_file);
    }
}
//=================================================================================================//
} // namespace SPH
