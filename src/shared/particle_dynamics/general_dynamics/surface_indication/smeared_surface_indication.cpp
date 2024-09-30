#include "smeared_surface_indication.h"

namespace SPH
{
//=================================================================================================//
SmearedSurfaceIndication::SmearedSurfaceIndication(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      indicator_(particles_->getVariableDataByName<int>("Indicator")),
      smeared_surface_(particles_->registerStateVariable<int>("SmearedSurface")) {}
//=================================================================================================//
void SmearedSurfaceIndication::interaction(size_t index_i, Real dt)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    smeared_surface_[index_i] = is_near_surface;
}
//=================================================================================================//
} // namespace SPH
