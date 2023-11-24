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
//Real PL = 0.00391; /**< X-direction domain. */
//Real PW = 0.02346; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 6.0;
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
Real vel_0 = 227.0;

//Real rho0_s = 2700.0; /**< Reference density. */
//Real poisson = 0.3;  /**< Poisson ratio. */
//Real Youngs_modulus = 78.2e9;
//Real yield_stress = 0.29e9;
//Real vel_0 = 373.0;

Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));

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
    explicit Column(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(0.0, 0.0, 0.6 * PW);
        Vecd halfsize_column(0.003, 0.003, 0.5 * PW);
        add<TriangleMeshShapeBrick>(halfsize_column, resolution, translation_column);
    }
};
//class Column : public ComplexShape
//{
//public:
//    explicit Column(const std::string& shape_name) : ComplexShape(shape_name)
//    {
//        Vecd translation_column(0.0, 0.0, 0.6 * PW);
//        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), inner_circle_radius,
//            0.5 * PW, resolution, translation_column);
//    }
//};

/**
 * application dependent initial condition
 */

class InitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][2] = -vel_0;
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

//class ColumnInitialCondition
//    : public fluid_dynamics::FluidInitialCondition
//{
//public:
//    explicit ColumnInitialCondition(RealBody& sph_body)
//        : fluid_dynamics::FluidInitialCondition(sph_body) {};
//
//protected:
//    void update(size_t index_i, Real dt)
//    {
//        /** initial velocity profile */
//        Real x = pos_[index_i][0] / PL;
//        if (x > 0.0)
//        {
//            vel_[index_i][1] = vf * c0 *
//                (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
//        }
//    };
//};