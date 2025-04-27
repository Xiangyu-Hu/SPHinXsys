#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
Neighbor<>::NeighborCriterion::NeighborCriterion(Neighbor<> &neighbor)
    : source_pos_(neighbor.source_pos_), target_pos_(neighbor.target_pos_),
      inv_h_(neighbor.inv_h_), kernel_size_square_(neighbor.kernel_size_square_) {}
//=================================================================================================//
} // namespace SPH
