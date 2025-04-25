#include "body_partition.h"

namespace SPH
{
//=================================================================================================//
BodyPartition::BodyPartition(SPHBody &sph_body, UnsignedInt adaptation_level)
    : BodyPartByID(sph_body), cell_linked_list_created_(false),
      adaptation_level_(adaptation_level),
      sph_adaptation_(sph_body.getSPHAdaptation()) {}
//=================================================================================================//
BodyPartitionTemporal::BodyPartitionTemporal(SPHBody &sph_body, UnsignedInt adaptation_level)
    : BodyPartition(sph_body, adaptation_level) {}
//=================================================================================================//
BaseCellLinkedList &BodyPartitionTemporal::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_.createCellLinkedList(sph_body_.getSPHSystemBounds(), base_particles_);
        cell_linked_list_created_ = true;
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
BodyPartitionSpatial::BodyPartitionSpatial(SPHBody &sph_body, UnsignedInt adaptation_level)
    : BodyPartition(sph_body, adaptation_level) {}
//=================================================================================================//
BaseCellLinkedList &BodyPartitionSpatial::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_.createRefinedCellLinkedList(
                adaptation_level_, sph_body_.getSPHSystemBounds(), base_particles_);
        cell_linked_list_created_ = true;
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
} // namespace SPH
