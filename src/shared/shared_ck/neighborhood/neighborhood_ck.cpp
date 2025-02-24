#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
Neighbor<>::NeighborCriterion::NeighborCriterion(Neighbor<> &neighbor)
    : source_pos_(neighbor.source_pos_), target_pos_(neighbor.target_pos_),
      cut_radius_square_(neighbor.getKernel().CutOffRadiusSqr()) {}
//=================================================================================================//
} // namespace SPH
