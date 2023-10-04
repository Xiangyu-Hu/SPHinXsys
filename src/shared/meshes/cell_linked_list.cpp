#include "cell_linked_list.h"
#include "adaptation.h"
#include "base_body.h"
#include "base_kernel.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
BaseCellLinkedList::
    BaseCellLinkedList(RealBody &real_body, SPHAdaptation &sph_adaptation)
    : BaseMeshField("CellLinkedList"),
      real_body_(real_body), kernel_(*sph_adaptation.getKernel()) {}
//=================================================================================================//
void BaseCellLinkedList::clearSplitCellLists(SplitCellLists &split_cell_lists)
{
    for (size_t i = 0; i < split_cell_lists.size(); i++)
        split_cell_lists[i].clear();
}

CellLinkedListKernel::CellLinkedListKernel(BaseParticles& particles, const DeviceVecd &meshLowerBound, DeviceReal gridSpacing,
                                           const DeviceArrayi &allGridPoints, const DeviceArrayi &allCells)
    : total_real_particles_(particles.total_real_particles_), list_data_pos_(particles.getDeviceVariableByName<DeviceVecd>("Position")),
      list_data_Vol_(particles.getDeviceVariableByName<DeviceReal>("Volume")), mesh_lower_bound_(meshLowerBound), grid_spacing_(gridSpacing),
      all_grid_points_(allGridPoints),  all_cells_(allCells), index_list_(allocateSharedData<size_t>(total_real_particles_)),
      index_head_list_(allocateSharedData<size_t>(allCells[0] * allCells[1])){}

void CellLinkedListKernel::clearCellLists() {
    // Only clear head list, since index list does not depend on its previous values
    copyDataToDevice(static_cast<size_t>(0), index_head_list_, all_cells_[0] * all_cells_[1]);
}

void CellLinkedListKernel::UpdateCellLists(SPH::BaseParticles &base_particles)
{
    clearCellLists();
    auto* pos_n = base_particles.getDeviceVariableByName<DeviceVecd>("Position");
    size_t total_real_particles = base_particles.total_real_particles_;
    executionQueue.getQueue().submit(
                  [&, mesh_lower_bound=mesh_lower_bound_, grid_spacing=grid_spacing_, all_grid_points=all_grid_points_,
                   all_cells=all_cells_, index_list=index_list_, index_head_list=index_head_list_](sycl::handler& cgh) {
                        cgh.parallel_for(total_real_particles, [=](sycl::id<1> id) {
                                 const size_t index_i = id.get(0);
                                 const auto cell_index = CellIndexFromPosition(pos_n[index_i], mesh_lower_bound,
                                                 grid_spacing, all_grid_points);
                                 const auto linear_cell_index = transferCellIndexTo1D(cell_index, all_cells);
                                 sycl::atomic_ref<size_t, sycl::memory_order_relaxed, sycl::memory_scope_device,
                                                  sycl::access::address_space::global_space>
                                     atomic_head_list(index_head_list[linear_cell_index]);
                                 /*
                                  * Insert index at the head of the list, the index previously at the top is then
                                  * used as the next one pointed by the new index.
                                  * Indices values are increased by 1 to let 0 be the value that indicates list termination.
                                  * If the cell list is empty (i.e. head == 0) then head will point to the new index and the
                                  * new index will point to 0 (i.e. the new index will be the first and last element of the cell list).
                                  *     index_head_list[linear_cell_index] = index_i+1  ---> index_list[index_i] = 0
                                  * Since the cell list order is not relevant, memory_order_relaxed will only ensure that each cell of
                                  * index_head_list gets a new index one at a time.
                                  */
                                 index_list[index_i] = atomic_head_list.exchange(index_i+1);
                        });
                  }).wait();
}

//=================================================================================================//
CellLinkedList::CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
                               RealBody &real_body, SPHAdaptation &sph_adaptation)
    : BaseCellLinkedList(real_body, sph_adaptation), Mesh(tentative_bounds, grid_spacing, 2),
      execution::DeviceExecutable<CellLinkedList, CellLinkedListKernel>(this, real_body_.getBaseParticles(),
                                                                        hostToDeviceVecd(mesh_lower_bound_), DeviceReal(grid_spacing_),
                                                                        hostToDeviceArrayi(all_grid_points_), hostToDeviceArrayi(all_cells_))
{
    allocateMeshDataMatrix();
    single_cell_linked_list_level_.push_back(this);
}
//=================================================================================================//
void CellLinkedList::UpdateCellLists(BaseParticles &base_particles)
{
    clearCellLists();
    StdLargeVec<Vecd> &pos_n = base_particles.pos_;
    size_t total_real_particles = base_particles.total_real_particles_;
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

    if (real_body_.getUseSplitCellLists())
    {
        updateSplitCellLists(real_body_.getSplitCellLists());
    }
}
//=================================================================================================//
StdLargeVec<size_t> &CellLinkedList::computingSequence(BaseParticles &base_particles)
{
//    return computingSequence(base_particles, execution::par);

    StdLargeVec<Vecd> &pos = base_particles.pos_;
    StdLargeVec<size_t> &sequence = base_particles.sequence_;
    size_t total_real_particles = base_particles.total_real_particles_;
    particle_for(execution::ParallelPolicy(), total_real_particles, [&](size_t i)
                 { sequence[i] = transferMeshIndexToMortonOrder(CellIndexFromPosition(pos[i])); });
    return sequence;
}
//=================================================================================================//
MultilevelCellLinkedList::MultilevelCellLinkedList(
    BoundingBox tentative_bounds, Real reference_grid_spacing,
    size_t total_levels, RealBody &real_body, SPHAdaptation &sph_adaptation)
    : MultilevelMesh<BaseCellLinkedList, CellLinkedList, RefinedMesh<CellLinkedList>>(
          tentative_bounds, reference_grid_spacing, total_levels, real_body, sph_adaptation),
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
void MultilevelCellLinkedList::
    insertParticleIndex(size_t particle_index, const Vecd &particle_position)
{
    size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    mesh_levels_[level]->insertParticleIndex(particle_index, particle_position);
}
//=================================================================================================//
void MultilevelCellLinkedList::
    InsertListDataEntry(size_t particle_index, const Vecd &particle_position, Real volumetric)
{
    size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
    mesh_levels_[level]->InsertListDataEntry(particle_index, particle_position, volumetric);
}
//=================================================================================================//
void MultilevelCellLinkedList::UpdateCellLists(BaseParticles &base_particles)
{
    for (size_t level = 0; level != total_levels_; ++level)
        mesh_levels_[level]->clearCellLists();

    StdLargeVec<Vecd> &pos_n = base_particles.pos_;
    size_t total_real_particles = base_particles.total_real_particles_;
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

    if (real_body_.getUseSplitCellLists())
    {
        updateSplitCellLists(real_body_.getSplitCellLists());
    }
}
//=================================================================================================//
StdLargeVec<size_t> &MultilevelCellLinkedList::computingSequence(BaseParticles &base_particles)
{
    StdLargeVec<Vecd> &pos = base_particles.pos_;
    StdLargeVec<size_t> &sequence = base_particles.sequence_;
    size_t total_real_particles = base_particles.total_real_particles_;
    particle_for(execution::ParallelPolicy(), total_real_particles,
                 [&](size_t i)
                 {
						 size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[i]));
						 sequence[i] = mesh_levels_[level]->transferMeshIndexToMortonOrder(
						 mesh_levels_[level]->CellIndexFromPosition(pos[i])); });

    return sequence;
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
