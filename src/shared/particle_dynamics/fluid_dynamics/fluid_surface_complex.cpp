#include "fluid_surface_complex.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
ColorFunctionGradientComplex::ColorFunctionGradientComplex(BaseInnerRelation &inner_relation,
                                                           BaseContactRelation &contact_relation)
    : ColorFunctionGradientInner(inner_relation), FluidContactOnly(contact_relation)
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
      indicator_(*particles_->getVariableByName<int>("Indicator")),
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