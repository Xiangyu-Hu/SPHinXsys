/**
 * @file 	case.h
 * @brief 	This is the case study the curvature of shell.
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

// //----------------------------------------------------------------------
// //	Basic geometry parameters and numerical setup.
// //----------------------------------------------------------------------
// const Real radius = 1;
// const Real length = 2;
// const Real resolution_ref = 0.1;
// const int particle_number_mid_surface = int(2.0 * radius * Pi / resolution_ref);
// const int particle_number_height = int(length / resolution_ref);
// const Real shell_thickness = resolution_ref;
// /** Domain bounds of the system. */
// BoundingBox system_domain_bounds(Vec3d(-radius, -radius, 0), Vec3d(radius, radius, length));
// //----------------------------------------------------------------------
// //	Define the body shape.
// //----------------------------------------------------------------------
// class ShellParticleGenerator : public SurfaceParticleGenerator
// {
//   public:
//     explicit ShellParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
//     virtual void initializeGeometricVariables() override
//     {
//         for (int i = 0; i < particle_number_mid_surface; i++)
//         {
//             for (int j = 0; j < particle_number_height; j++)
//             {
//                 Real theta = (i + 0.5) * 2 * Pi / (Real)particle_number_mid_surface;
//                 Real x = radius * cos(theta);
//                 Real y = radius * sin(theta);
//                 Real z = length * j / (Real)particle_number_height + 0.5 * resolution_ref;
//                 initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
//                                                        resolution_ref * resolution_ref);
//                 Vec3d n_0 = Vec3d(x / radius, y / radius, 0.0);
//                 initializeSurfaceProperties(n_0, shell_thickness);
//             }
//         }
//     }
// };

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real radius = 1.0;
const Real resolution_ref = 0.1;
const Real shell_thickness = resolution_ref;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius, -radius, -radius), Vec3d(radius, radius, radius));
//----------------------------------------------------------------------
//	Material Parameters
//----------------------------------------------------------------------
Real rho0_s = 10.0;
Real Youngs_modulus = 1e3;
Real poisson = 0.3;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ShellParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        Real dphi = resolution_ref / radius;
        int particle_number_phi = int(2 * Pi / dphi);
        for (int i = 0; i < particle_number_phi; i++)
        {
            Real phi = dphi * (i + 0.5);
            Real radius_phi = radius * sin(phi);
            Real z = radius * cos(phi);
            int particle_number_theta = int(2 * Pi * radius_phi / resolution_ref);
            for (int j = 0; j < particle_number_theta; j++)
            {
                Real theta = (j + 0.5) * 2 * Pi / (Real)particle_number_theta;
                Real x = radius_phi * cos(theta);
                Real y = radius_phi * sin(theta);
                initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
                                                       resolution_ref * resolution_ref);
                Vec3d n_0 = Vec3d(x / radius, y / radius, z / radius);
                initializeSurfaceProperties(n_0, shell_thickness);
            }
        }
    }
};
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ControlDisplacement : public LocalDynamics, public thin_structure_dynamics::ShellDataSimple
{
  public:
    ControlDisplacement(SPHBody &sph_body, Real moving_v)
        : LocalDynamics(sph_body),
          thin_structure_dynamics::ShellDataSimple(sph_body),
          n_(particles_->n_), vel_(particles_->vel_), moving_v_(moving_v){};

  protected:
    StdLargeVec<Vecd> &n_;
    StdLargeVec<Vecd> &vel_;
    Real moving_v_;

    inline void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = moving_v_ * n_[index_i];
    };
};
