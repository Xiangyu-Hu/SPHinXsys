#include "solid_dynamics_variable.h"

namespace SPH
{
//=============================================================================================//
Displacement::Displacement(SPHBody &sph_body)
    : BaseDerivedVariable<Vecd>(sph_body, "Displacement"), SolidDataSimple(sph_body),
      LocalDynamics(sph_body),
      pos_(*base_particles_.getVariableByName<Vecd>("Position")),
      pos0_(*particles_->registerSharedVariableFrom<Vecd>("InitialPosition", "Position")) {}
//=============================================================================================//
void Displacement::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = pos_[index_i] - pos0_[index_i];
}
//=============================================================================================//
OffsetInitialPosition::OffsetInitialPosition(SPHBody &sph_body, Vecd &offset)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body),
      offset_(offset),
      pos_(*base_particles_.getVariableByName<Vecd>("Position")),
      pos0_(*particles_->registerSharedVariableFrom<Vecd>("InitialPosition", "Position")) {}
//=============================================================================================//
void OffsetInitialPosition::update(size_t index_i, Real dt)
{
    pos_[index_i] += offset_;
    pos0_[index_i] += offset_;
}
//=============================================================================================//
TranslationAndRotation::TranslationAndRotation(SPHBody &sph_body, Transform &transform)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body), transform_(transform),
      pos_(*base_particles_.getVariableByName<Vecd>("Position")),
      pos0_(*particles_->registerSharedVariableFrom<Vecd>("InitialPosition", "Position")) {}
//=============================================================================================//
void TranslationAndRotation::update(size_t index_i, Real dt)
{
    pos_[index_i] = transform_.shiftFrameStationToBase(pos_[index_i]);
    pos0_[index_i] = transform_.shiftFrameStationToBase(pos0_[index_i]);
}
//=============================================================================================//
GreenLagrangeStrain::GreenLagrangeStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Matd>(sph_body, "GreenLagrangeStrain"), ElasticSolidDataSimple(sph_body),
      LocalDynamics(sph_body), F_(*particles_->getVariableByName<Matd>("DeformationGradient")) {}
//=============================================================================================//
void GreenLagrangeStrain::update(size_t index_i, Real dt)
{
    Matd F = F_[index_i];
    derived_variable_[index_i] = 0.5 * (F.transpose() * F - Matd::Identity());
}
//=============================================================================================//
VonMisesStress::VonMisesStress(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStress"), ElasticSolidDataSimple(sph_body),
      LocalDynamics(sph_body), rho0_(sph_body_.base_material_->ReferenceDensity()),
      rho_(*particles_->getVariableByName<Real>("Density")),
      F_(*particles_->getVariableByName<Matd>("DeformationGradient")),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial())) {}
//=============================================================================================//
VonMisesStrain::VonMisesStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStrain"),
      ElasticSolidDataSimple(sph_body), LocalDynamics(sph_body),
      F_(*particles_->getVariableByName<Matd>("DeformationGradient")) {}
//=============================================================================================//
VonMisesStrainDynamic::VonMisesStrainDynamic(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStrainDynamic"),
      ElasticSolidDataSimple(sph_body), LocalDynamics(sph_body),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial())),
      poisson_ratio_(elastic_solid_.PoissonRatio()),
      F_(*particles_->getVariableByName<Matd>("DeformationGradient")) {}
//=============================================================================================//
MidSurfaceVonMisesStress::MidSurfaceVonMisesStress(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "MidSurfaceVonMisesStress"), ShellSolidDataSimple(sph_body),
      LocalDynamics(sph_body),
      mid_surface_cauchy_stress_(*particles_->getVariableByName<Matd>("MidSurfaceCauchyStress")) {}
//=================================================================================================//
} // namespace SPH
