 /**
 * @file play_simbody.cpp
 * @brief This is a simple example for play with simbody.
 */

#include "sphinxsys.h"
#include "UdfMotion.h"

using namespace SimTK;
using std::cout; 
using std::endl;

int main() {
  try {   
    /** Create the system, with subsystems for the bodies and some forces. */
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
	SPH::StateEngine state_engine(system);
    //Force::UniformGravity gravity(forces, matter, Vec3(10, Real(-9.8), 3));
    /** Create the body and some artwork for it. */
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), 
                               DecorativeSphere(Real(0.1)).setColor(Red));
    /** Add an instance of the body to the multibody system by connecting
      it to Ground via a pin mobilizer.
     */
    MobilizedBody::Pin pendulum1(matter.updGround(), 
                                Transform(/*x45,*/Vec3(0,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));
    MyMotion 			motion1(pendulum1, Pi / 40.0, 10.0, 2.0 * Pi, 0.0);
    MobilizedBody::Pin pendulum1b(pendulum1, 
                                Transform(/*x45,*/Vec3(0,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));
    MyMotion 			motion2(pendulum1b, Pi / 30.0, 10.0, 2.0 * Pi, -0.5 * Pi);
    MobilizedBody::Pin pendulum1c(pendulum1b, 
                                Transform(/*x45,*/Vec3(0,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));
    MyMotion 			motion3(pendulum1c, Pi / 20.0, 10.0, 2.0 * Pi, -0.8 * Pi);
    /** Visualize with default options; ask for a report every 1/30 of a second
      to match the Visualizer's default 30 frames per second rate.
     */
    // Visualizer not included in dependency-free version
    /*
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, Real(1./30)));
    */
    /** Initialize the system and state. */
    system.realizeTopology();
    State state = system.getDefaultState();
    Real restart = 0.0;
    if(restart == 0.0)
    {
        pendulum1.setOneQ(state, 0, 0.0);
        pendulum1b.setOneQ(state, 0,0.0);
        pendulum1c.setOneQ(state, 0,0.0);
        system.realize(state);
    }else
    {
        state = state_engine.readAndSetStateInfoFromXml(int(restart), system);
	    state.setTime(restart);
    }

    // Visualizer not included in dependency-free version
    /*
    viz.report(state);
    */
    /** Simulate it. */
    cout << "Hit ENTER to run a short simulation ...";
    getchar();
    
	RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    if (restart == 0.0)
    {
        ts.stepTo(5.0);
	    const SimTK::State& state_for_output = integ.getState();
	    state_engine.writeStateInfoToXml(5, state_for_output);
    }
	ts.stepTo(10.0);
  } catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }

    return 0;
}