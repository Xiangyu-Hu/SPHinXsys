
#include "outflow_boundary.h"

namespace SPH
{
    OutflowConditionWithFace::OutflowConditionWithFace(FluidBody &fluid_body, BodyRegionByCellsWithFace &body_region)
        : PartSimpleDynamicsByCellsWithFace(fluid_body, body_region),
          DataDelegateSimple<FluidBody, FluidParticles, Fluid>(fluid_body),
          pos_n_(fluid_body.base_particles_->pos_n_)
    {
    }

    void OutflowConditionWithFace::Update(size_t index_i, Real dt)
    {
        auto d = body_region_.getSignedDistance(pos_n_[index_i]);
        if (d < 0.0)
        {
            base_particles_->switchToBufferParticle(index_i);
        }
        else if (d < body_region_.getRegionWidth())
        {
            UpdateOutletRegionParticles(index_i, dt);
        }
    }

    ModifiedDoNothingConditionWithFace::ModifiedDoNothingConditionWithFace(FluidBody &fluid_body, BodyRegionByCellsWithFace &body_region)
        : OutflowConditionWithFace(fluid_body, body_region),
          vel_n_(fluid_body.base_particles_->vel_n_),
          cell_linked_list_(dynamic_cast<CellLinkedList *>(fluid_body.cell_linked_list_)),
          relation_inner_(&fluid_body)
    {
    }


} // namespace SPH