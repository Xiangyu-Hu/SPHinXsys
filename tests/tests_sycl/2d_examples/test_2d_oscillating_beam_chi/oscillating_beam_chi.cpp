/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D oscillating beam — chi (high-level API) version   *
 * ----------------------------------------------------------------------------*
 * Rewrites the standard 2D oscillating beam test using the high-level         *
 * SPHSimulation facade.  The beam is a SaintVenantKirchhoff solid clamped     *
 * at its left insert.  An analytical initial velocity profile (linear-theory   *
 * first mode) is applied, and the tip displacement is written to file.        *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "sph_simulation.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real total_physical_time = 1.0; /**< TOTAL SIMULATION TIME*/
Real PL = 0.2;                  // beam length
Real PH = 0.02;                 // for thick plate; =0.01 for thin plate
Real SL = 0.06;                 // depth of the insert
Real global_resolution = PH / 10.0;
Real BW = global_resolution * 4; // boundary width, at least three particles
//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;       // Poisson ratio
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
//----------------------------------------------------------------------
//	Application-dependent initial velocity profile (first beam mode).
//----------------------------------------------------------------------
class LinearProfile : public ReturnFunction<Vecd>
{
    Real vf_ = vf;
    Real c0_;
    Real kl_ = kl;
    Real M_ = M;
    Real N_ = N;
    Real Q_ = Q;
    Real PL_ = PL;

  public:
    explicit LinearProfile(Real c0) : c0_(c0) {}

    Vecd operator()(const Vec2d &position)
    {
        Real x = position[0] / PL_;
        Vecd result = Vec2d::Zero();
        if (x > 0.0)
        {
            result[1] = vf_ * c0_ *
                        (M_ * (cos(kl_ * x) - cosh(kl_ * x)) -
                         N_ * (sin(kl_ * x) - sinh(kl_ * x))) /
                        Q_;
        }
        return result;
    }
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Geometry: beam base box and beam column box (in halfsize + translation).
    //----------------------------------------------------------------------
    Vec2d base_halfsize(0.5 * (SL + BW), 0.5 * PH + BW);
    Vec2d base_translation(-0.5 * (SL + BW), 0.0);
    Vec2d column_halfsize(0.5 * (PL + SL), 0.5 * PH);
    Vec2d column_translation(0.5 * (PL - SL), 0.0);

    //----------------------------------------------------------------------
    //	Build simulation with the high-level SPHSimulation facade.
    //----------------------------------------------------------------------
    SPHSimulation sim;
    sim.defineDomain(Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0),
                     global_resolution);

    //----------------------------------------------------------------------
    //	Configure the solid block: geometry, material, constraint, damping.
    //----------------------------------------------------------------------
    auto &beam = sim.addSolidBlock("BeamBody")
                     .addBox(base_halfsize, base_translation)
                     .addBox(column_halfsize, column_translation)
                     .materialSVK(rho0_s, Youngs_modulus, poisson)
                     .constrainBox(base_halfsize, base_translation)
                     .subtractFromConstraint(column_halfsize, column_translation)
                     .withNumericalDamping();

    // Retrieve the reference sound speed after the material has been configured
    // and use it to scale the initial velocity profile.
    beam.initialVelocity(LinearProfile(beam.getReferenceSoundSpeed()));

    //----------------------------------------------------------------------
    //	Observer at the beam tip.
    //----------------------------------------------------------------------
    sim.addObserver("BeamObserver", Vec2d(PL, 0.0));

    //----------------------------------------------------------------------
    //	Run the simulation.
    //----------------------------------------------------------------------
    sim.run(total_physical_time);

    return 0;
}
