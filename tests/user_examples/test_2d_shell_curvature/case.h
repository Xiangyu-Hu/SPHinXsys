/**
 * @file 	case.h
 * @brief 	This is the case study the curvature of shell.
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real radius_x = 2;
const Real radius_y = 1;
const int particle_number_mid_surface = 100;
const Real resolution_ref = M_PI * (radius_x + radius_y) / Real(particle_number_mid_surface);
const Real shell_thickness = resolution_ref;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-radius_x, -radius_y), Vec2d(radius_x, radius_y));
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** Particle generator and constraint boundary for shell baffle. */
// x=R*cos(theta), y=R*sin(theta)
class ShellParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        Real theta = 0;
        while (theta < 2 * Pi)
        {
            Real x = radius_x * cos(theta);
            Real y = radius_y * sin(theta);
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_ref);
            Vec2d normal_direction = Vec2d(x / radius_x / radius_x, y / radius_y / radius_y);
            normal_direction /= normal_direction.norm();
            initializeSurfaceProperties(normal_direction, shell_thickness);
            Real ds = radius_x * radius_y * sqrt(pow(x / radius_x / radius_x, 2) + pow(y / radius_y / radius_y, 2));
            theta += resolution_ref / ds;
        }
    }
};