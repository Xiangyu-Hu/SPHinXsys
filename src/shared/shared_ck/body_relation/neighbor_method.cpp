#include "neighbor_method.h"

#include "body_partition.h"
namespace SPH
{
//=================================================================================================//
Real SmoothingLength<Base>::getSmoothingLength(
    const Fixed &fixed, BodyPartitionSpatial &body_partition_spatial)
{
    return body_partition_spatial.getReferenceSmoothingLength();
}
//=================================================================================================//
} // namespace SPH
