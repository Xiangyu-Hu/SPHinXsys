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
        : BaseDerivedVariable<Vecd>(sph_body, "Displacemant"), SolidDataSimple(sph_body),
          pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_) {}
    //=============================================================================================//
    void Displacement::operator()(size_t index_i, Real dt)
    {
        derived_variable_[index_i] = pos_n_[index_i] - pos_0_[index_i];
    }
    //=============================================================================================//
    VonMisesStress::VonMisesStress(SPHBody &sph_body)
        : BaseDerivedVariable<Real>(sph_body, "VonMiesStress"), ElasticSolidDataSimple(sph_body),
          rho0_(particles_->rho0_), rho_n_(particles_->rho_n_),
          F_(particles_->F_), stress_PK1_(particles_->stress_PK1_) {}
    //=============================================================================================//
    VonMisesStrain::VonMisesStrain(SPHBody &sph_body)
        : BaseDerivedVariable<Real>(sph_body, "VonMiesStrain"),
          ElasticSolidDataSimple(sph_body), F_(particles_->F_) {}
    //=================================================================================================//
}
