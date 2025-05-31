#include "body_partition.h"

namespace SPH
{
//=================================================================================================//
BodyPartition::BodyPartition(SPHBody &sph_body, UnsignedInt partition_adapt_level)
    : sph_body_(sph_body), sph_adaptation_(sph_body.getSPHAdaptation()),
      base_particles_(sph_body.getBaseParticles()),
      cell_linked_list_created_(false), present_adapt_level_(partition_adapt_level),
      dv_adapt_level_(base_particles_.registerStateVariableOnly<int>("AdaptLevel")) {}
//=================================================================================================//
std::string BodyPartition::getName()
{
    return sph_body_.getName() + "Partition" + std::to_string(present_adapt_level_);
}
//=================================================================================================//
BodyPartitionTemporal::BodyPartitionTemporal(SPHBody &sph_body, UnsignedInt partition_adapt_level)
    : BodyPartition(sph_body, partition_adapt_level) {}
//=================================================================================================//
BaseCellLinkedList &BodyPartitionTemporal::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_.createCellLinkedList(sph_body_.getSPHSystemBounds(), base_particles_);
        cell_linked_list_created_ = true;
        cell_linked_list_ptr_.get()->setName(getName() + "CellLinkedList");
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
BodyPartitionSpatial::BodyPartitionSpatial(SPHBody &sph_body, UnsignedInt present_adapt_level_)
    : BodyPartition(sph_body, present_adapt_level_) {}
//=================================================================================================//
BaseCellLinkedList &BodyPartitionSpatial::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_.createRefinedCellLinkedList(
                present_adapt_level_, sph_body_.getSPHSystemBounds(), base_particles_);
        cell_linked_list_created_ = true;
        cell_linked_list_ptr_.get()->setName(getName() + "CellLinkedList");
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
} // namespace SPH
