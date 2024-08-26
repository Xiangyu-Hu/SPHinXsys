#include "cell_linked_list.h"
#include "adaptation.h"
#include "base_kernel.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "mesh_iterators.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
BaseCellLinkedList::
    BaseCellLinkedList(BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : BaseMeshField("CellLinkedList"), kernel_(*sph_adaptation.getKernel()),
      base_particles_(base_particles), dv_pos_(base_particles.getVariableByName<Vecd>("Position")) {}
//=================================================================================================//
SplitCellLists *BaseCellLinkedList::getSplitCellLists()
{
    std::cout << "\n Error: SplitCellList not defined!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return nullptr;
}
//=================================================================================================//
void BaseCellLinkedList::setUseSplitCellLists()
{
    std::cout << "\n Error: SplitCellList not defined!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
};
//=================================================================================================//
void BaseCellLinkedList::clearSplitCellLists(SplitCellLists &split_cell_lists)
{
    for (size_t i = 0; i < split_cell_lists.size(); i++)
        split_cell_lists[i].clear();
}
//=================================================================================================//
CellLinkedList::CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
                               BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : BaseCellLinkedList(base_particles, sph_adaptation), Mesh(tentative_bounds, grid_spacing, 2),
      use_split_cell_lists_(false), cell_index_lists_(nullptr), cell_data_lists_(nullptr),
      cell_offset_list_size_(NumberOfCells() + 1),
      index_list_size_(SMAX(base_particles.ParticlesBound(), cell_offset_list_size_)),
      dv_particle_index_(base_particles.registerDiscreteVariableOnly<UnsignedInt>("ParticleIndex", index_list_size_)),
      dv_cell_offset_(base_particles.registerDiscreteVariableOnly<UnsignedInt>("CellOffset", cell_offset_list_size_))
{
    allocateMeshDataMatrix();
    single_cell_linked_list_level_.push_back(this);
    size_t number_of_split_cell_lists = pow(3, Dimensions);
    split_cell_lists_.resize(number_of_split_cell_lists);
}
//=================================================================================================//
void CellLinkedList ::allocateMeshDataMatrix()
{
    size_t number_of_all_cells = transferMeshIndexTo1D(all_cells_, all_cells_);
    cell_index_lists_ = new ConcurrentIndexVector[number_of_all_cells];
    cell_data_lists_ = new ListDataVector[number_of_all_cells];
}
//=================================================================================================//
void CellLinkedList ::deleteMeshDataMatrix()
{
    delete[] cell_index_lists_;
    delete[] cell_data_lists_;
}
//=================================================================================================//
void CellLinkedList::clearCellLists()
{
    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](const Arrayi &cell_index)
                      {
                          getCellDataList(cell_index_lists_, cell_index).clear();
                      });
}
//=================================================================================================//
void CellLinkedList::UpdateCellListData(BaseParticles &base_particles)
{
    Vecd *pos = base_particles.ParticlePositions();
    mesh_parallel_for(
        MeshRange(Arrayi::Zero(), all_cells_),
        [&](const Arrayi &cell_index)
        {
            ListDataVector &cell_data_list = getCellDataList(cell_data_lists_, cell_index);
            cell_data_list.clear();
            ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
            for (size_t s = 0; s != cell_list.size(); ++s)
            {
                size_t index = cell_list[s];
                cell_data_list.emplace_back(std::make_pair(index, pos[index]));
            }
        });
}
//=================================================================================================//
void CellLinkedList::updateSplitCellLists(SplitCellLists &split_cell_lists)
{
    clearSplitCellLists(split_cell_lists);
    mesh_parallel_for(
        MeshRange(Arrayi::Zero(), all_cells_),
        [&](const Arrayi &cell_index)
        {
            ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
            size_t real_particles_in_cell = cell_list.size();
            if (real_particles_in_cell != 0)
            {
                split_cell_lists[transferMeshIndexTo1D(3 * Arrayi::Ones(), mod(cell_index, 3))]
                    .push_back(&cell_list);
            }
        });
}
//=================================================================================================//
void CellLinkedList::UpdateCellLists(BaseParticles &base_particles)
{
    clearCellLists();
    Vecd *pos_n = base_particles.ParticlePositions();
    size_t total_real_particles = base_particles.TotalRealParticles();
    parallel_for(
        IndexRange(0, total_real_particles),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                insertParticleIndex(i, pos_n[i]);
            }
        },
        ap);

    UpdateCellListData(base_particles);

    if (use_split_cell_lists_)
    {
        updateSplitCellLists(split_cell_lists_);
    }
}
//=================================================================================================//
void CellLinkedList ::insertParticleIndex(size_t particle_index, const Vecd &particle_position)
{
    Arrayi cell_index = CellIndexFromPosition(particle_position);
    getCellDataList(cell_index_lists_, cell_index).emplace_back(particle_index);
}
//=================================================================================================//
void CellLinkedList ::InsertListDataEntry(size_t particle_index, const Vecd &particle_position)
{
    Arrayi cell_index = CellIndexFromPosition(particle_position);
    getCellDataList(cell_data_lists_, cell_index)
        .emplace_back(std::make_pair(particle_index, particle_position));
}
//=================================================================================================//
ListData CellLinkedList::findNearestListDataEntry(const Vecd &position)
{
    Real min_distance_sqr = MaxReal;
    ListData nearest_entry = std::make_pair(MaxSize_t, MaxReal * Vecd::Ones());

    Arrayi cell = CellIndexFromPosition(position);
    mesh_for_each(
        Arrayi::Zero().max(cell - Arrayi::Ones()),
        all_cells_.min(cell + 2 * Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            ListDataVector &target_particles = getCellDataList(cell_data_lists_, cell_index);
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
    return nearest_entry;
}
//=================================================================================================//
void CellLinkedList::
    tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included)
{
    mesh_parallel_for(
        MeshRange(Arrayi::Zero(), all_cells_),
        [&](const Arrayi &cell_index)
        {
            bool is_included = false;
            mesh_for_each(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                all_cells_.min(cell_index + 2 * Arrayi::Ones()),
                [&](const Arrayi &neighbor_cell_index)
                {
                    if (check_included(CellPositionFromIndex(neighbor_cell_index), grid_spacing_))
                    {
                        is_included = true;
                    }
                });
            if (is_included == true)
                cell_lists.push_back(&getCellDataList(cell_index_lists_, cell_index));
        });
}
//=================================================================================================//
UnsignedInt CellLinkedList::computingSequence(Vecd &position, size_t index_i)
{
    return transferMeshIndexToMortonOrder(CellIndexFromPosition(position));
}
//=================================================================================================//
MultilevelCellLinkedList::MultilevelCellLinkedList(BoundingBox tentative_bounds,
                                                   Real reference_grid_spacing, size_t total_levels,
                                                   BaseParticles &base_particles, SPHAdaptation &sph_adaptation)
    : MultilevelMesh<BaseCellLinkedList, CellLinkedList>(
          tentative_bounds, reference_grid_spacing, total_levels, base_particles, sph_adaptation),
      h_ratio_(DynamicCast<ParticleWithLocalRefinement>(this, &sph_adaptation)->h_ratio_)
{
}
//=================================================================================================//
size_t MultilevelCellLinkedList::getMeshLevel(Real particle_cutoff_radius)
{
    for (size_t level = total_levels_; level != 0; --level)
        if (particle_cutoff_radius - mesh_levels_[level - 1]->GridSpacing() < Eps)
            return level - 1; // jump out the loop!

    std::cout << "\n Error: CellLinkedList level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void MultilevelCellLinkedList::insertParticleIndex(size_t particle_index, const Vecd &particle_position)
{
    size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    mesh_levels_[level]->insertParticleIndex(particle_index, particle_position);
}
//=================================================================================================//
void MultilevelCellLinkedList::InsertListDataEntry(size_t particle_index, const Vecd &particle_position)
{
    size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    mesh_levels_[level]->InsertListDataEntry(particle_index, particle_position);
}
//=================================================================================================//
void MultilevelCellLinkedList::UpdateCellLists(BaseParticles &base_particles)
{
    for (size_t level = 0; level != total_levels_; ++level)
        mesh_levels_[level]->clearCellLists();

    Vecd *pos_n = base_particles.ParticlePositions();
    size_t total_real_particles = base_particles.TotalRealParticles();
    // rebuild the corresponding particle list.
    parallel_for(
        IndexRange(0, total_real_particles),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                insertParticleIndex(i, pos_n[i]);
            }
        },
        ap);

    for (size_t level = 0; level != total_levels_; ++level)
    {
        mesh_levels_[level]->UpdateCellListData(base_particles);
    }
}
//=================================================================================================//
UnsignedInt MultilevelCellLinkedList::computingSequence(Vecd &position, size_t index_i)
{
    size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[index_i]));
    return mesh_levels_[level]->transferMeshIndexToMortonOrder(
        mesh_levels_[level]->CellIndexFromPosition(position));
}
//=================================================================================================//
void MultilevelCellLinkedList::
    tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included)
{
    for (size_t l = 0; l != total_levels_; ++l)
    {
        mesh_levels_[l]->tagBodyPartByCell(cell_lists, check_included);
    }
}
//=================================================================================================//
} // namespace SPH
