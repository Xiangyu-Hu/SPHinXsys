/**
 * @file twisting_column_chi.cpp
 * @brief Rewrites the 3D twisting column test using the high-level SPHSimulation
 *        facade.  The column is a NeoHookean solid with one end kinematically
 *        fixed (the holder).  A sinusoidal angular-velocity profile is applied
 *        as an initial condition, and the tip position is written to file.
 * @author Chi Zhang and Xiangyu Hu
 * @ref DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"
#include "sph_simulation.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real total_physical_time = 0.5; /**< TOTAL SIMULATION TIME*/
Real PL = 6.0;                  /**< X-direction size. */
Real PH = 1.0;                  /**< Y-direction size. */
Real PW = 1.0;                  /**< Z-direction size. */
Real particle_spacing_ref = PH / 10.0;
Real BW = particle_spacing_ref * 0.0; /**< no wall boundary in this case. */
Real SL = particle_spacing_ref * 1.0; /**< Length of the holder is one layer particle. */
Vecd halfsize_column(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_column(0.5 * (PL - SL), 0.0, 0.0);
Vecd halfsize_holder(0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
Vecd translation_holder(-0.5 * (SL + BW), 0.0, 0.0);
Vec3d domain_lower_bound(-SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
Vec3d domain_upper_bound(PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 1100.0; /**< Reference density. */
Real poisson = 0.45;  /**< Poisson ratio. */
Real Youngs_modulus = 1.7e7;
Real angular_0 = -300.0; /**< Initial angular velocity amplitude. */
//------------------------------------------------------------------------------
//	Sinusoidal twisting velocity profile applied as an initial condition.
//------------------------------------------------------------------------------
class VelocityProfile : public ReturnFunction<Vecd>
{
    Real angular_0_ = angular_0;
    Real PL_ = PL;

  public:
    VelocityProfile() {}

    Vecd operator()(const Vec3d &position)
    {
        Real x = position[0];
        Real y = position[1];
        Real z = position[2];
        Real angular_velocity = angular_0_ * sin((M_PI * x) / (2.0 * PL_));
        Real local_radius = sqrt(pow(y, 2) + pow(z, 2));
        Real angular = atan2(y, z);

        if (x > 0.0)
        {
            return Vecd(0.0,
                        angular_velocity * local_radius * cos(angular),
                        -angular_velocity * local_radius * sin(angular));
        }
        return Vecd::Zero();
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build simulation with the high-level SPHSimulation facade.
    //----------------------------------------------------------------------
    SPHSimulation sim;
    sim.defineDomain(domain_lower_bound, domain_upper_bound, particle_spacing_ref);

    //----------------------------------------------------------------------
    //	Configure the solid block: geometry, material, constraint, velocity.
    //----------------------------------------------------------------------
    sim.addSolidBlock("Column")
        .addBox(halfsize_column, translation_column)
        .addBox(halfsize_holder, translation_holder)
        .materialNeoHookean(rho0_s, Youngs_modulus, poisson)
        .constrainBox(halfsize_holder, translation_holder)
        .initialVelocity(VelocityProfile());

    //----------------------------------------------------------------------
    //	Observer at the column tip.
    //----------------------------------------------------------------------
    sim.addObserver("MyObserver", Vecd(PL, 0.0, 0.0));

    //----------------------------------------------------------------------
    //	Run the simulation.
    //----------------------------------------------------------------------
    sim.run(total_physical_time);

    return 0;
}
