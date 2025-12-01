/**
 * @file play_simbody.cpp
 * @brief This is a simple example for play with simbody.
 */

#include "sphinxsys.h"

#include "UdfMotion.h"

int main(int ac, char *av[])
{
    SPH::BoundingBoxd system_domain_bounds(SPH::Vec3d(0, 0, 0), SPH::Vec3d(1, 1, 1));
    SPH::SPHSystem sph_system(system_domain_bounds, 1.0);
    sph_system.handleCommandlineOptions(ac, av);
    double Pi = SPH::Pi;
    /** Create the multi_body_system, with subsystems for the bodies and some forces. */
    SimTK::MultibodySystem multi_body_system;
    SimTK::SimbodyMatterSubsystem matter(multi_body_system);
    SimTK::GeneralForceSubsystem forces(multi_body_system);
    SPH::SimbodyStateEngine state_engine(sph_system, multi_body_system);

    // Force::UniformGravity gravity(forces, matter, SimTK::Vec3(10, Real(-9.8), 3));
    /** Create the body and some artwork for it. */
    SimTK::Body::Rigid pendulumBody(SimTK::MassProperties(1.0, SimTK::Vec3(0), SimTK::Inertia(1)));
    pendulumBody.addDecoration(SimTK::Transform(),
                               SimTK::DecorativeSphere(0.1).setColor(SimTK::Red));
                               /** Add an instance of the body to the multibody system by connecting
      it to Ground via a pin mobilizer.
     */
    SimTK::MobilizedBody::Pin pendulum1(matter.updGround(),
                                        SimTK::Transform(/*x45,*/ SimTK::Vec3(0, -1, 0)),
                                        pendulumBody,
                                        SimTK::Transform(SimTK::Vec3(0, 1, 0)));
    MyMotion motion1(pendulum1, Pi / 40.0, 10.0, 2.0 * Pi, 0.0);
    SimTK::MobilizedBody::Pin pendulum1b(pendulum1,
                                         SimTK::Transform(/*x45,*/ SimTK::Vec3(0, -1, 0)),
                                         pendulumBody,
                                         SimTK::Transform(SimTK::Vec3(0, 1, 0)));
    MyMotion motion2(pendulum1b, Pi / 30.0, 10.0, 2.0 * Pi, -0.5 * Pi);
    SimTK::MobilizedBody::Pin pendulum1c(pendulum1b,
                                         SimTK::Transform(/*x45,*/ SimTK::Vec3(0, -1, 0)),
                                         pendulumBody,
                                         SimTK::Transform(SimTK::Vec3(0, 1, 0)));
    MyMotion motion3(pendulum1c, Pi / 20.0, 10.0, 2.0 * Pi, -0.8 * Pi);
    multi_body_system.realizeTopology();
    SimTK::State state = multi_body_system.getDefaultState();
    pendulum1.setOneQ(state, 0, 0.0);
    pendulum1b.setOneQ(state, 0, 0.0);
    pendulum1c.setOneQ(state, 0, 0.0);
    /** Restart the system if needed. */
    size_t restart = sph_system.RestartStep();
    size_t step = restart; // Current step
    if (restart != 0)
    {
        state_engine.readStateFromXml(restart, state);
        state.setTime(double(restart));
    }
    multi_body_system.realize(state);
    SimTK::RungeKuttaMersonIntegrator integ(multi_body_system);
    SimTK::TimeStepper ts(multi_body_system, integ);
    ts.initialize(state);
    size_t num_steps = 10; // Number of steps to simulate
    while (step < num_steps)
    {
        step++;
        ts.stepTo(double(step));
        state_engine.writeStateToXml(step, integ);
    }
    return 0;
}