#include "fluid_surface_complex.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
FreeSurfaceIndicationComplex::
    FreeSurfaceIndicationComplex(BaseInnerRelation &inner_relation,
                                 BaseContactRelation &contact_relation, Real threshold)
    : FreeSurfaceIndicationInner(inner_relation, threshold), FluidContactData(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
    }
}
//=================================================================================================//
FreeSurfaceIndicationComplex::
    FreeSurfaceIndicationComplex(ComplexRelation &complex_relation, Real threshold)
    : FreeSurfaceIndicationComplex(complex_relation.getInnerRelation(),
                                   complex_relation.getContactRelation(), threshold) {}
//=================================================================================================//
ColorFunctionGradientComplex::ColorFunctionGradientComplex(BaseInnerRelation &inner_relation,
                                                           BaseContactRelation &contact_relation)
    : ColorFunctionGradientInner(inner_relation), FluidContactData(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
ColorFunctionGradientComplex::ColorFunctionGradientComplex(ComplexRelation &complex_relation)
    : ColorFunctionGradientComplex(complex_relation.getInnerRelation(),
                                   complex_relation.getContactRelation()) {}
//=================================================================================================//
SurfaceNormWithWall::SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle)
    : LocalDynamics(contact_relation.getSPHBody()), FSIContactData(contact_relation),
      contact_angle_(contact_angle),
      surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
      surface_norm_(*particles_->getVariableByName<Vecd>("SurfaceNormal")),
      pos_div_(*particles_->getVariableByName<Real>("PositionDivergence"))
{
    particle_spacing_ = contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    smoothing_length_ = contact_relation.getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_n_.push_back(&(contact_particles_[k]->n_));
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//