/**
 * @file 	oscillating_plate_UL.cpp
 * @brief 	This is the test case for the hourglass manuscript.
 * @details  We consider vibration deformation of a square plate under initial vertical velocity field.
 * @author 	Shuaihao Zhang, Dong Wu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real PL = 0.4;                                          /** Length of the square plate. */
Real PH = 0.4;                                          /** Width of the square plate. */
Real PT = 0.01;                                         /** Thickness of the square plate. */
int particle_number = 3;                                /** Particle number in the direction of the thickness. */
Real particle_spacing_ref = PT / (Real)particle_number; /** Initial reference particle spacing. */
int particle_number_PL = PL / particle_spacing_ref;
int particle_number_PH = PH / particle_spacing_ref;
int BWD = 1;
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -PT / 2), Vec3d(PL + BW, PH + BW, PT / 2));
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.5 * PH, 0.0), Vecd(-BW, -BW, 0.0)};
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_s = 1000.0;          /** Normalized density. */
Real Youngs_modulus = 100.0e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;            /** Poisson ratio. */
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
Real gravity_g = 0.0;

Real governing_vibration_integer_x = 2.0;
Real governing_vibration_integer_y = 2.0;
Real U_ref = 1.0; // Maximum velocity

namespace SPH
{
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class Plate;
template <>
class ParticleGenerator<BaseParticles, Plate> : public ParticleGenerator<BaseParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles){};
    virtual void prepareGeometricData() override
    {
        for (int k = 0; k < particle_number; k++)
        {
            for (int i = 0; i < particle_number_PL + 2 * BWD; i++)
            {
                for (int j = 0; j < particle_number_PH + 2 * BWD; j++)
                {
                    Real x = particle_spacing_ref * i - BW + particle_spacing_ref * 0.5;
                    Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
                    Real z = particle_spacing_ref * (k - ((particle_number - 1.0) / 2.0));
                    addPositionAndVolumetricMeasure(
                        Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref * particle_spacing_ref);
                }
            }
        }
    }
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry(){};

  private:
    void tagManually(size_t index_i)
    {
        if ((base_particles_.ParticlePositions()[index_i][2] < 0.5 * particle_spacing_ref) &&
            (base_particles_.ParticlePositions()[index_i][2] > -0.5 * particle_spacing_ref) &&
            (base_particles_.ParticlePositions()[index_i][0] < 0.0 ||
             base_particles_.ParticlePositions()[index_i][0] > PL ||
             base_particles_.ParticlePositions()[index_i][1] < 0.0 ||
             base_particles_.ParticlePositions()[index_i][1] > PH))
        {
            body_part_particles_.push_back(index_i);
        }
    };
};

/** Define the initial condition. */
class BeamInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit BeamInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        vel_[index_i][2] = sin(governing_vibration_integer_x * Pi * pos_[index_i][0] / PL) *
                           sin(governing_vibration_integer_y * Pi * pos_[index_i][1] / PH);
    };
};
} // namespace SPH

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    SolidBody plate_body(sph_system, makeShared<DefaultShape>("PlateBody"));
    plate_body.defineMaterial<GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
    plate_body.generateParticles<BaseParticles, Plate>();

    ObserverBody plate_observer(sph_system, "PlateObserver");
    plate_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation plate_body_inner(plate_body);
    ContactRelation plate_observer_contact(plate_observer, {&plate_body});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxillary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<BeamInitialCondition> initial_velocity(plate_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(plate_body_inner);

    InteractionDynamics<continuum_dynamics::ShearAccelerationRelaxation> plate_shear_acceleration(plate_body_inner);
    Dynamics1Level<continuum_dynamics::Integration1stHalf> plate_pressure_relaxation(plate_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerDissipativeRiemann> plate_density_relaxation(plate_body_inner);
    Dynamics1Level<continuum_dynamics::ShearStressRelaxation> plate_shear_stress_relaxation(plate_body_inner);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(plate_body, U_ref, 0.2);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(plate_body, 0.4);
    BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
    SimpleDynamics<FixedInAxisDirection> constrain_holder(boundary_geometry, Vecd(1.0, 1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center(plate_body, Vecd(1.0, 1.0, 0.0));
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(plate_body, "VonMisesStress");
    write_states.addToWrite<Real>(plate_body, "VonMisesStress");
    write_states.addToWrite<Real>(plate_body, "Pressure");
    write_states.addToWrite<Real>(plate_body, "Density");
    RestartIO restart_io(sph_system);
    ObservedQuantityRecording<Vecd> write_plate_displacement("Position", plate_observer_contact);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_kinetic_energy(plate_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_velocity.exec();
    constrain_holder.exec();
    corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
    }
    GlobalStaticVariables::physical_time_ = 0.0;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 500;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 0.02;
    Real output_period = end_time / 50.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(0);
    write_plate_displacement.writeToFile(0);
    write_kinetic_energy.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            Real relaxation_time = 0.0;
            Real advection_dt = fluid_advection_time_step.exec();

            while (relaxation_time < advection_dt)
            {
                Real acoustic_dt = fluid_acoustic_time_step.exec();
                plate_shear_stress_relaxation.exec(acoustic_dt);
                plate_pressure_relaxation.exec(acoustic_dt);
                constrain_holder.exec(acoustic_dt);
                plate_density_relaxation.exec(acoustic_dt);
                plate_shear_acceleration.exec(acoustic_dt);
                number_of_iterations++;
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                GlobalStaticVariables::physical_time_ += acoustic_dt;
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << "N=" << number_of_iterations << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	advection_dt: "
                              << advection_dt << "	acoustic_dt: "
                              << acoustic_dt << "\n";
                    if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != sph_system.RestartStep())
                        restart_io.writeToFile(Real(number_of_iterations));
                }
            }
            plate_body.updateCellLinkedList();
            plate_body_inner.updateConfiguration();
            corrected_configuration.exec();
        }
        write_plate_displacement.writeToFile(number_of_iterations);
        write_kinetic_energy.writeToFile(number_of_iterations);
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TickCount::interval_t tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_kinetic_energy.testResult();
    }
    return 0;
}
