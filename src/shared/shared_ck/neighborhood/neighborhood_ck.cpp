#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
Neighbor<>::SearchMethod::SearchMethod(Vecd *source_pos, Vecd *target_pos, Real cut_radius_square)
    : source_pos_(source_pos), target_pos_(target_pos),
      cut_radius_square_(cut_radius_square) {}
//=================================================================================================//
} // namespace SPH

