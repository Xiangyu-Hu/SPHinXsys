/* ----------------------------------------------------------------------------*
 *                    SPHinXsys: 2D oscillation beam                           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases for understanding SPH method for    *
 * solid simulation based on updated Lagrangian method.                        *
 * A generalized hourglass control method is used here.                         *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu                              *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;  // beam length
Real PH = 0.02; // for thick plate
Real SL = 0.06; // depth of the insert
Real global_resolution = PH / 10;
Real BW = global_resolution * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
                                 Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;       // Poisson ratio
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
Real U_ref = vf * c0 * (M * (cos(kl) - cosh(kl)) - N * (sin(kl) - sinh(kl))) / Q;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a beam base shape
std::vector<Vecd> beam_base_shape{
    Vecd(-SL - BW, -PH / 2 - BW), Vecd(-SL - BW, PH / 2 + BW), Vecd(0.0, PH / 2 + BW),
    Vecd(0.0, -PH / 2 - BW), Vecd(-SL - BW, -PH / 2 - BW)};
// a beam shape
std::vector<Vecd> beam_shape{
    Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(-SL, -PH / 2)};
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class BeamInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit BeamInitialCondition(RealBody &beam_column)
        : fluid_dynamics::FluidInitialCondition(beam_column) {};

  protected:
    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        Real x = pos_[index_i][0] / PL;
        if (x > 0.0)
        {
            vel_[index_i][1] = vf * c0 *
                               (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
        }
    };
};
//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::sub);
    return multi_polygon;
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody beam_body(sph_system, makeShared<Beam>("BeamBody"));
    beam_body.defineMaterial<GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody beam_observer(sph_system, "BeamObserver");
    beam_observer.getSPHAdaptation().resetAdaptationRatios(1.15, 2.0);
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    InnerRelation beam_body_inner(beam_body);
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    /** initial condition */
    SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correction_matrix(beam_body_inner);
    Dynamics1Level<continuum_dynamics::Integration1stHalf> beam_pressure_relaxation(beam_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerDissipativeRiemann> beam_density_relaxation(beam_body_inner);
    InteractionWithUpdate<continuum_dynamics::ShearStressRelaxationHourglassControl1stHalf> beam_shear_stress(beam_body_inner);
    InteractionDynamics<continuum_dynamics::ShearStressRelaxationHourglassControl2ndHalf> beam_shear_acceleration(beam_body_inner);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> beam_volume_update(beam_body);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> advection_time_step(beam_body, U_ref, 0.2);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> acoustic_time_step(beam_body, 0.4);
    // clamping a solid body part.
    BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
    SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_beam_states(beam_body);
    write_beam_states.addToWrite<Real>(beam_body, "Density");
    write_beam_states.addToWrite<Real>(beam_body, "Pressure");
    SimpleDynamics<continuum_dynamics::VonMisesStress> beam_von_mises_stress(beam_body);
    write_beam_states.addToWrite<Real>(beam_body, "VonMisesStress");
    ObservedQuantityRecording<Vecd> write_beam_tip_displacement("Position", beam_observer_contact);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_beam_kinetic_energy(beam_body);
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    beam_initial_velocity.exec();
    correction_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup computing time-step controls.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 1.0;
    Real End_Time = T0;
    Real D_Time = End_Time / 100; /**< Time period for data observing */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    write_beam_states.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);
    write_beam_kinetic_energy.writeToFile(0);
    // computation loop starts
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < D_Time)
        {
            Real relaxation_time = 0.0;
            Real advection_dt = advection_time_step.exec();
            beam_volume_update.exec();
            while (relaxation_time < advection_dt)
            {
                Real acoustic_dt = acoustic_time_step.exec();
                beam_pressure_relaxation.exec(acoustic_dt);
                constraint_beam_base.exec();
                beam_shear_stress.exec(acoustic_dt);
                beam_shear_acceleration.exec(acoustic_dt);
                beam_density_relaxation.exec(acoustic_dt);
                ite++;
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                physical_time += acoustic_dt;
                if (ite % 500 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	advection_dt: "
                              << advection_dt << "	acoustic_dt: "
                              << acoustic_dt << "\n";
                }
            }
            beam_body.updateCellLinkedList();
            beam_body_inner.updateConfiguration();
            correction_matrix.exec();
        }
        beam_von_mises_stress.exec();
        write_beam_tip_displacement.writeToFile(ite);
        write_beam_kinetic_energy.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_beam_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_beam_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_beam_kinetic_energy.testResult();
    }
    return 0;
}