/**
 * @file solid_particles_variable.cpp
 * @brief Definition of functions declared in solid_particles_variable.h
 * @author	Xiangyu Hu
 */
#include "solid_particles_variable.h"
#include "elastic_solid.h"

namespace SPH
{
    //=============================================================================================//
    Displacement::Displacement(SPHBody &sph_body)
        : BaseDerivedVariable<Vecd>(sph_body, "Displacement"), SolidDataSimple(sph_body),
          pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_) {}
    //=============================================================================================//
    void Displacement::update(size_t index_i, Real dt)
    {
        derived_variable_[index_i] = pos_n_[index_i] - pos_0_[index_i];
    }
    //=============================================================================================//
    OffsetInitialPosition::
        OffsetInitialPosition(SPHBody &sph_body, Vecd &offset)
        : SolidDataSimple(sph_body), offset_(offset),
          pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_) {}
    //=============================================================================================//
    void OffsetInitialPosition::update(size_t index_i, Real dt)
    {
        pos_n_[index_i] += offset_;
        pos_0_[index_i] += offset_;
    }
    //=============================================================================================//
    TranslationAndRotation::
        TranslationAndRotation(SPHBody &sph_body, Transformd &transform)
        : SolidDataSimple(sph_body), transform_(transform),
          pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_) {}
	//=============================================================================================//
	void TranslationAndRotation::update(size_t index_i, Real dt)
	{
		pos_n_[index_i] = transform_.shiftFrameStationToBase(pos_n_[index_i]);
		pos_0_[index_i] = transform_.shiftFrameStationToBase(pos_0_[index_i]);
	}
   //=============================================================================================//
    NormalDirectionFromBodyShape::
        NormalDirectionFromBodyShape(SPHBody &sph_body)
        : SolidDataSimple(sph_body), body_shape_(*sph_body.body_shape_),
          pos_n_(particles_->pos_n_), n_(particles_->n_), n_0_(particles_->n_0_) {}
    //=============================================================================================//
    void NormalDirectionFromBodyShape::update(size_t index_i, Real dt)
    {
        Vecd normal_direction = body_shape_.findNormalDirection(pos_n_[index_i]);
        n_[index_i] = normal_direction;
        n_0_[index_i] = normal_direction;
    }
    //=============================================================================================//
    NormalDirectionFromShapeAndOp::
        NormalDirectionFromShapeAndOp(SPHBody &sph_body, const std::string &shape_name)
        : SolidDataSimple(sph_body),
          shape_and_op_(DynamicCast<ComplexShape>(this, sph_body.body_shape_)->getShapeAndOpByName(shape_name)),
          shape_(shape_and_op_->first),
          switch_sign_(shape_and_op_->second == ShapeBooleanOps::add ? 1.0 : -1.0),
          pos_n_(particles_->pos_n_), n_(particles_->n_), n_0_(particles_->n_0_) {}
    //=============================================================================================//
    void NormalDirectionFromShapeAndOp::update(size_t index_i, Real dt)
    {
        Vecd normal_direction = switch_sign_ * shape_->findNormalDirection(pos_n_[index_i]);
        n_[index_i] = normal_direction;
        n_0_[index_i] = normal_direction;
    }
    //=============================================================================================//
    VonMisesStress::VonMisesStress(SPHBody &sph_body)
        : BaseDerivedVariable<Real>(sph_body, "VonMisesStress"), ElasticSolidDataSimple(sph_body),
          rho0_(particles_->rho0_), rho_n_(particles_->rho_n_),
          F_(particles_->F_), stress_PK1_(particles_->stress_PK1_) {}
    //=============================================================================================//
    VonMisesStrain::VonMisesStrain(SPHBody &sph_body)
        : BaseDerivedVariable<Real>(sph_body, "VonMisesStrain"),
          ElasticSolidDataSimple(sph_body), F_(particles_->F_) {}
    //=================================================================================================//
}
