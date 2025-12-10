/* ----------------------------------------------------------------------------*
 *                    SPHinXsys: 2D spinning plate                             *
 * ----------------------------------------------------------------------------*
 * This case is designed to test tensile instability and                       *
 * angular momentum conservation properties.                                   *
 * @author Shuaihao Zhang and Xiangyu Hu                                       *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 1.0; // square length
Real PH = PL;
Real resolution_ref = PH / 20;
Real BW = resolution_ref * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-PL, -PH), Vec2d(PL, PH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.1e3;         // reference density
Real Youngs_modulus = 1.7e7; // reference Youngs modulus
Real poisson = 0.45;         // Poisson ratio
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
Real angular_0 = -50.0; // rad/s
Real U_ref = angular_0 * 0.5 * sqrt(2);
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a square shape
std::vector<Vecd> square_shape{
    Vecd(-PL / 2, -PH / 2), Vecd(-PL / 2, PH / 2), Vecd(PL / 2, PH / 2), Vecd(PL / 2, -PH / 2), Vecd(-PL / 2, -PH / 2)};
// Square observer location
StdVec<Vecd> observation_location = {Vecd(PL / 2, PH / 2)};
//----------------------------------------------------------------------
//	Define the square body
//----------------------------------------------------------------------
class Square : public MultiPolygonShape
{
  public:
    explicit Square(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(square_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class SquareInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit SquareInitialCondition(RealBody &square_column)
        : fluid_dynamics::FluidInitialCondition(square_column) {};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real x = pos_[index_i][0];
        Real y = pos_[index_i][1];
        Real local_radius = sqrt(pow(x, 2) + pow(y, 2));
        Real angle = atan2(x, y);
        vel_[index_i][0] = angular_0 * local_radius * cos(angle);
        vel_[index_i][1] = -angular_0 * local_radius * sin(angle);
    };
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody square_body(sph_system, makeShared<Square>("SquareBody"));
    square_body.defineMaterial<GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
    square_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody square_observer(sph_system, "SquareObserver");
    square_observer.getSPHAdaptation().resetAdaptationRatios(1.15, 2.0);
    square_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    ContactRelation square_observer_contact(square_observer, {&square_body});
    InnerRelation square_body_inner(square_body);
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    /** initial condition */
    SimpleDynamics<SquareInitialCondition> square_initial_velocity(square_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correction_matrix(square_body_inner);
    Dynamics1Level<continuum_dynamics::Integration1stHalf> square_pressure_relaxation(square_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerDissipativeRiemann> square_density_relaxation(square_body_inner);
    InteractionWithUpdate<continuum_dynamics::ShearStressRelaxationHourglassControl1stHalf> square_shear_stress(square_body_inner);
    InteractionDynamics<continuum_dynamics::ShearStressRelaxationHourglassControl2ndHalf> square_shear_acceleration(square_body_inner);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> square_volume_update(square_body);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> acoustic_time_step(square_body, 0.1);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_square_states(square_body);
    write_square_states.addToWrite<Real>(square_body, "Density");
    write_square_states.addToWrite<Real>(square_body, "Pressure");
    SimpleDynamics<continuum_dynamics::VonMisesStress> square_von_mises_stress(square_body);
    write_square_states.addToWrite<Real>(square_body, "VonMisesStress");
    ObservedQuantityRecording<Vecd> write_square_tip_displacement("Position", square_observer_contact);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_square_kinetic_energy(square_body);
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    square_initial_velocity.exec();
    correction_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup computing time-step controls.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 0.5;
    Real End_Time = T0;
    Real D_Time = End_Time / 50; /**< Time period for data observing */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    write_square_states.writeToFile(0);
    write_square_tip_displacement.writeToFile(0);
    write_square_kinetic_energy.writeToFile(0);
    // computation loop starts
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < D_Time)
        {
            square_volume_update.exec();
            Real acoustic_dt = acoustic_time_step.exec();
            square_pressure_relaxation.exec(acoustic_dt);
            square_shear_stress.exec(acoustic_dt);
            square_shear_acceleration.exec(acoustic_dt);
            square_density_relaxation.exec(acoustic_dt);
            ite++;
            integration_time += acoustic_dt;
            physical_time += acoustic_dt;
            if (ite % 500 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	advection_dt: "
                          << acoustic_dt << "\n";
            }

            square_body.updateCellLinkedList();
            square_body_inner.updateConfiguration();
            correction_matrix.exec();
        }
        square_von_mises_stress.exec();
        write_square_tip_displacement.writeToFile(ite);
        write_square_kinetic_energy.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_square_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_square_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_square_kinetic_energy.testResult();
    }
    return 0;
}