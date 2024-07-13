#include "cell_linked_list.h"
#include "adaptation.h"
#include "base_kernel.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
BaseCellLinkedList::
    BaseCellLinkedList(SPHAdaptation &sph_adaptation)
    : BaseMeshField("CellLinkedList"),
      kernel_(*sph_adaptation.getKernel()) {}
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

CellLinkedListKernel::CellLinkedListKernel(const DeviceVecd &meshLowerBound, DeviceReal gridSpacing,
                                           const DeviceArrayi &allGridPoints, const DeviceArrayi &allCells)
    : total_real_particles_(0), list_data_pos_(nullptr),  // will be initialized at first UpdateCellLists execution
      mesh_lower_bound_(allocateDeviceData<DeviceVecd>(1)), grid_spacing_(allocateDeviceData<DeviceReal>(1)),
      all_grid_points_(allocateDeviceData<DeviceArrayi>(1)), all_cells_(allocateDeviceData<DeviceArrayi>(1)),
      linear_cell_size_(VecdFoldProduct(allCells)), offset_cell_size_(allocateDeviceData<DeviceInt>(linear_cell_size_+1)),
      curr_cell_size_(allocateDeviceData<DeviceInt>(linear_cell_size_))
{
    copyDataToDevice(&meshLowerBound, mesh_lower_bound_, 1)
        .add(copyDataToDevice(&gridSpacing, grid_spacing_, 1))
        .add(copyDataToDevice(&allGridPoints, all_grid_points_, 1))
        .add(copyDataToDevice(&allCells, all_cells_, 1)).wait();
}

CellLinkedListKernel::~CellLinkedListKernel()
{
    freeDeviceData(mesh_lower_bound_);
    freeDeviceData(grid_spacing_);
    freeDeviceData(all_grid_points_);
    freeDeviceData(all_cells_);
}

execution::ExecutionEvent CellLinkedListKernel::clearCellLists()
{
    return std::move(copyDataToDevice(static_cast<DeviceInt>(0), offset_cell_size_, linear_cell_size_+1)
                .add(copyDataToDevice(static_cast<DeviceInt>(0), curr_cell_size_, linear_cell_size_)));
}

execution::ExecutionEvent CellLinkedListKernel::UpdateCellLists(const SPH::BaseParticles &base_particles)
{
    if(!list_data_pos_)  // initialize list data with base_particles at first execution
    {
        total_real_particles_ = base_particles.total_real_particles_;
        list_data_pos_ = base_particles.getDeviceVariableByName<DeviceVecd>("Position");
        particle_id_list_ = allocateDeviceData<DeviceInt>(total_real_particles_);
    }
    auto clear_event = clearCellLists();

    DeviceInt total_real_particles = base_particles.total_real_particles_;
    auto determine_cells_size_event =  executionQueue.getQueue().submit(
        [&, mesh_lower_bound = mesh_lower_bound_, grid_spacing = grid_spacing_, all_grid_points = all_grid_points_,
         all_cells = all_cells_, offset_cell_size = offset_cell_size_, pos_n = list_data_pos_](sycl::handler &cgh)
                      {
                             cgh.depends_on(clear_event.getEventList());
                             cgh.parallel_for(executionQueue.getUniformNdRange(total_real_particles),
                                              [=](sycl::nd_item<1> item) {
                                                  if(item.get_global_id() < total_real_particles)
                                                  {
                                                        const DeviceInt index_i = item.get_global_id();
                                                        const auto cell_index = CellIndexFromPosition(pos_n[index_i], *mesh_lower_bound,
                                                                                                     *grid_spacing, *all_grid_points);
                                                        const auto linear_cell_index = transferCellIndexTo1D(cell_index, *all_cells);
                                                        sycl::atomic_ref<DeviceInt, sycl::memory_order_relaxed, sycl::memory_scope_device,
                                                                        sycl::access::address_space::global_space>
                                                           atomic_cell_size(offset_cell_size[linear_cell_index]);
                                                        ++atomic_cell_size;
                                                  }
                                         }); });

    auto exclusive_scan_cells_size_event = executionQueue.getQueue().submit(
        [&, offset_cell_size = offset_cell_size_, num_cells = linear_cell_size_ + 1](sycl::handler &cgh) {
            cgh.depends_on(determine_cells_size_event);
            cgh.parallel_for(executionQueue.getUniformNdRange(executionQueue.getWorkGroupSize()),
                             [=](sycl::nd_item<1> item) {
                                          if(item.get_group_linear_id() == 0) {
                                              sycl::joint_exclusive_scan(item.get_group(), offset_cell_size,
                                                                         offset_cell_size + num_cells,
                                                                         offset_cell_size, sycl::plus<>{});
                                          }
                                      });
                 });

    return executionQueue.getQueue().submit(
        [&, mesh_lower_bound = mesh_lower_bound_, grid_spacing = grid_spacing_, all_grid_points = all_grid_points_,
        all_cells = all_cells_, particle_id_list = particle_id_list_, offset_cell_size = offset_cell_size_,
        curr_cell_size = curr_cell_size_, pos_n = list_data_pos_] (sycl::handler &cgh)
        {
            cgh.depends_on(exclusive_scan_cells_size_event);
            cgh.parallel_for(executionQueue.getUniformNdRange(total_real_particles),
                             [=](sycl::nd_item<1> item)
                             {
                                 if(item.get_global_id() < total_real_particles)
                                 {
                                     const DeviceInt index_i = item.get_global_id();
                                     const auto cell_index = CellIndexFromPosition(pos_n[index_i], *mesh_lower_bound,
                                                                                   *grid_spacing, *all_grid_points);
                                     const auto linear_cell_index = transferCellIndexTo1D(cell_index, *all_cells);
                                     sycl::atomic_ref<DeviceInt, sycl::memory_order_relaxed, sycl::memory_scope_device,
                                                      sycl::access::address_space::global_space>
                                         atomic_current_cell_size(curr_cell_size[linear_cell_index]);

                                     particle_id_list[offset_cell_size[linear_cell_index] + atomic_current_cell_size++] = index_i;
                                 }
                             });
        });
}
//=================================================================================================//
DeviceInt *CellLinkedListKernel::computingSequence(BaseParticles &baseParticles)
{
    auto *pos = baseParticles.getDeviceVariableByName<DeviceVecd>("Position");
    auto *sequence = baseParticles.sequence_device_;
    DeviceInt total_real_particles = baseParticles.total_real_particles_;
    executionQueue.getQueue().submit(
                                 [&, mesh_lower_bound = mesh_lower_bound_, grid_spacing = grid_spacing_,
                                  all_grid_points = all_grid_points_](sycl::handler &cgh)
                                 {
                                     cgh.parallel_for(executionQueue.getUniformNdRange(total_real_particles), [=](sycl::nd_item<1> item)
                                                      {
                                                          DeviceInt i = item.get_global_id();
                                                          if(i < total_real_particles)
                                                              sequence[i] = BaseMesh::transferMeshIndexToMortonOrder(
                                                                  CellIndexFromPosition(pos[i], *mesh_lower_bound,
                                                                                        *grid_spacing, *all_grid_points)); });
                                 })
        .wait_and_throw();
    return sequence;
}
//=================================================================================================//
CellLinkedList::CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
                               SPHAdaptation &sph_adaptation)
    : BaseCellLinkedList(sph_adaptation), Mesh(tentative_bounds, grid_spacing, 2),
      use_split_cell_lists_(false), device_kernel(hostToDeviceVecd(mesh_lower_bound_), DeviceReal(grid_spacing_),
                                                  hostToDeviceArrayi(all_grid_points_), hostToDeviceArrayi(all_cells_))
{
    allocateMeshDataMatrix();
    single_cell_linked_list_level_.push_back(this);
    size_t number_of_split_cell_lists = pow(3, Dimensions);
    split_cell_lists_.resize(number_of_split_cell_lists);
}
//=================================================================================================//
execution::ExecutionEvent CellLinkedList::UpdateCellLists(BaseParticles &base_particles)
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

    if (use_split_cell_lists_)
    {
        updateSplitCellLists(split_cell_lists_);
    }

    return {};
}
//=================================================================================================//
size_t *CellLinkedList::computingSequence(BaseParticles &base_particles)
{
    StdLargeVec<Vecd> &pos = base_particles.pos_;
    StdLargeVec<size_t> &sequence = base_particles.sequence_;
    size_t total_real_particles = base_particles.total_real_particles_;
    particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles), [&](size_t i)
                 { sequence[i] = transferMeshIndexToMortonOrder(CellIndexFromPosition(pos[i])); });
    return sequence.data();
}
//=================================================================================================//
MultilevelCellLinkedList::MultilevelCellLinkedList(
    BoundingBox tentative_bounds, Real reference_grid_spacing,
    size_t total_levels, SPHAdaptation &sph_adaptation)
    : MultilevelMesh<BaseCellLinkedList, CellLinkedList, RefinedMesh<CellLinkedList>>(
          tentative_bounds, reference_grid_spacing, total_levels, sph_adaptation),
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
execution::ExecutionEvent MultilevelCellLinkedList::UpdateCellLists(BaseParticles &base_particles)
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

    return {};
}
//=================================================================================================//
size_t *MultilevelCellLinkedList::computingSequence(BaseParticles &base_particles)
{
    StdLargeVec<Vecd> &pos = base_particles.pos_;
    StdLargeVec<size_t> &sequence = base_particles.sequence_;
    size_t total_real_particles = base_particles.total_real_particles_;
    particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
                 [&](size_t i)
                 {
						 size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[i]));
						 sequence[i] = mesh_levels_[level]->transferMeshIndexToMortonOrder(
						 mesh_levels_[level]->CellIndexFromPosition(pos[i])); });

    return sequence.data();
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
