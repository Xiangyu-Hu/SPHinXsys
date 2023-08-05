#include "wetting_coupled_spatial_temporal_method.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
NonWettingSurfaceIndication::
    NonWettingSurfaceIndication(BaseInnerRelation &inner_relation,
                                               BaseContactRelation &contact_relation, Real threshold, Real criterion)
    : FreeSurfaceIndicationComplex(inner_relation, contact_relation, threshold), wetting_criterion(criterion)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_phi_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Phi"));
    }
}
//=================================================================================================//
NonWettingSurfaceIndication::
    NonWettingSurfaceIndication(ComplexRelation &complex_relation, Real threshold, Real criterion)
    : NonWettingSurfaceIndication(complex_relation.getInnerRelation(),
                                                 complex_relation.getContactRelation(), threshold, criterion) {}

//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//