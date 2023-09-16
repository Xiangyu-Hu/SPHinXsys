#include "solid_particles_variable.h"
#include "elastic_solid.h"

namespace SPH
{
//=============================================================================================//
Displacement::Displacement(SPHBody &sph_body)
    : BaseDerivedVariable<Vecd>(sph_body, "Displacement"), SolidDataSimple(sph_body),
      LocalDynamics(sph_body), pos_(particles_->pos_), pos0_(particles_->pos0_) {}
//=============================================================================================//
void Displacement::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = pos_[index_i] - pos0_[index_i];
}
//=============================================================================================//
OffsetInitialPosition::
    OffsetInitialPosition(SPHBody &sph_body, Vecd &offset)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body),
      offset_(offset), pos_(particles_->pos_), pos0_(particles_->pos0_) {}
//=============================================================================================//
void OffsetInitialPosition::update(size_t index_i, Real dt)
{
    pos_[index_i] += offset_;
    pos0_[index_i] += offset_;
}
//=============================================================================================//
TranslationAndRotation::TranslationAndRotation(SPHBody &sph_body, Transform &transform)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body), transform_(transform),
      pos_(particles_->pos_), pos0_(particles_->pos0_) {}
//=============================================================================================//
void TranslationAndRotation::update(size_t index_i, Real dt)
{
    pos_[index_i] = transform_.shiftFrameStationToBase(pos_[index_i]);
    pos0_[index_i] = transform_.shiftFrameStationToBase(pos0_[index_i]);
}
//=============================================================================================//
NormalDirectionFromBodyShape::NormalDirectionFromBodyShape(SPHBody &sph_body)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body), body_shape_(*sph_body.body_shape_),
      pos_(particles_->pos_), n_(particles_->n_), n0_(particles_->n0_) {}
//=============================================================================================//
void NormalDirectionFromBodyShape::update(size_t index_i, Real dt)
{
    Vecd normal_direction = body_shape_.findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
}
//=============================================================================================//
NormalDirectionFromShapeAndOp::
    NormalDirectionFromShapeAndOp(SPHBody &sph_body, const std::string &shape_name)
    : SolidDataSimple(sph_body), LocalDynamics(sph_body),
      shape_and_op_(DynamicCast<ComplexShape>(this, sph_body.body_shape_)->getShapeAndOpByName(shape_name)),
      shape_(shape_and_op_->first),
      switch_sign_(shape_and_op_->second == ShapeBooleanOps::add ? 1.0 : -1.0),
      pos_(particles_->pos_), n_(particles_->n_), n0_(particles_->n0_) {}
//=============================================================================================//
void NormalDirectionFromShapeAndOp::update(size_t index_i, Real dt)
{
    Vecd normal_direction = switch_sign_ * shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
}
//=============================================================================================//
GreenLagrangeStrain::GreenLagrangeStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Matd>(sph_body, "GreenLagrangeStrain"), ElasticSolidDataSimple(sph_body),
      LocalDynamics(sph_body), F_(particles_->F_) {}
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
      rho_(particles_->rho_), F_(particles_->F_), elastic_solid_(particles_->elastic_solid_) {}
//=============================================================================================//
VonMisesStrain::VonMisesStrain(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStrain"),
      ElasticSolidDataSimple(sph_body), LocalDynamics(sph_body) {}
//=============================================================================================//
void VonMisesStrain::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = particles_->getVonMisesStrain(index_i);
}
//=============================================================================================//
VonMisesStrainDynamic::VonMisesStrainDynamic(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "VonMisesStrainDynamic"),
      ElasticSolidDataSimple(sph_body), LocalDynamics(sph_body),
      poisson_ratio_(particles_->elastic_solid_.PoissonRatio()) {}
//=============================================================================================//
void VonMisesStrainDynamic::update(size_t index_i, Real dt)
{
    derived_variable_[index_i] = particles_->getVonMisesStrainDynamic(index_i, poisson_ratio_);
}
//=============================================================================================//
MidSurfaceVonMisesStressofShells::MidSurfaceVonMisesStressofShells(SPHBody &sph_body)
    : BaseDerivedVariable<Real>(sph_body, "MidSurfaceVonMisesStressofShells"), ShellSolidDataSimple(sph_body),
      LocalDynamics(sph_body), mid_surface_cauchy_stress_(particles_->mid_surface_cauchy_stress_) {}
//=================================================================================================//
} // namespace SPH
