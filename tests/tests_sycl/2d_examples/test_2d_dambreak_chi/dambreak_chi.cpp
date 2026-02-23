/**
 * @file dambreak_chi.cpp
 * @brief 2D dambreak example using the high-level SPHSimulation user API.
 * @details This file demonstrates the simplified user entry point for SPH simulations.
 *          All internal SPH setup (bodies, relations, dynamics, I/O, time stepping)
 *          is handled by the SPHSimulation facade.
 * @author Xiangyu Hu
 */
#include "sphinxsys.h"
#include "sph_simulation.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                        /**< Reference density of fluid. */
Real gravity_g = 1.0;                     /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH);  /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                  /**< Artificial sound speed. */
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    SPHSimulation sim;

    sim.createDomain(DL, DH, particle_spacing_ref);

    sim.addFluidBlock("WaterBody")
        .rectangle(LL, LH)
        .material(rho0_f, c_f);

    sim.addWall("WallBoundary")
        .hollowBox(DL, DH, BW);

    sim.enableGravity(0.0, -gravity_g);

    sim.addObserver("FluidObserver", Vec2d(DL, 0.2));

    sim.useSolver()
        .dualTimeStepping()
        .freeSurfaceCorrection();

    sim.run(20.0);

    return 0;
}
