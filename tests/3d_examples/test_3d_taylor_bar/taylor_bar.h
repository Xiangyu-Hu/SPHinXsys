/**
 * @file 	taylor_bar.h
 * @brief 	This is the case setup for plastic taylor bar.
 * @author 	xiaojing tang Chi Zhang and Xiangyu Hu
 */
#pragma once

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 0.0032; /**< X-direction domain. */
Real PW = 0.0324; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 10.0;
/** YOU can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real SL = particle_spacing_ref * 4.0; /**< Length of the holder is one layer particle. */
Real inner_circle_radius = PL;

Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; /**< Reference density. */
Real poisson = 0.35;  /**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;

class Wall : public ComplexShape
{
  public:
    explicit Wall(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
        Vecd translation_holder(0.0, 0.0, -0.5 * SL);
        add<TriangleMeshShapeBrick>(halfsize_holder, resolution, translation_holder);
    }
};

/** Define the body. */
class Column : public ComplexShape
{
  public:
    explicit Column(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(0.0, 0.0, 0.6 * PW);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), inner_circle_radius,
                                       0.5 * PW, resolution, translation_column);
    }
};

/**
 * application dependent initial condition
 */

class InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][2] = -227.0;
    }
};

// define an observer body
class ColumnObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit ColumnObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        positions_.push_back(Vecd(0.0, 0.0, PW));
        positions_.push_back(Vecd(PL, 0.0, 0.0));
    }
};
