/**
 * @file 	test_3d_dynamic_plate.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider deformation of a square plate under stepping load.
 * @author 	Dong Wu and Xiangyu Hu
 * @ref		A New Continuum-Based Thick Shell Finite Element for Soft Biological Tissues in Dynamics:
 *			Part 1 - Preliminary Benchmarking Using Classic Verification Experiments
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 0.254;                                        /** Length of the square plate. */
Real PH = 0.254;                                        /** Width of the square plate. */
Real PT = 0.0127;                                       /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                       /** Pseudo-normal. */
int particle_number = 40;                               /** Particle number in the direction of the length */
Real particle_spacing_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -0.5 * particle_spacing_ref),
                                 Vec3d(PL + BW, PH + BW, 0.5 * particle_spacing_ref));
// Observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.5 * PH, 0.0), Vecd(-BW, -BW, 0.0)};
/** For material properties of the solid. */
Real rho0_s = 1.0;                /** Normalized density. */
Real Youngs_modulus = 68.94e9;    /** Normalized Youngs Modulus. */
Real poisson = 0.3;               /** Poisson ratio. */
Real physical_viscosity = 2000.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 2068427.0; /** Total distributed load. */
Real time_to_full_external_force = 0.0;

Real gravitational_acceleration = 0.0;

/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit PlateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        // the plate and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {
            for (int j = 0; j < (particle_number + 2 * BWD); j++)
            {
                Real x = particle_spacing_ref * i - BW + particle_spacing_ref * 0.5;
                Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
                initializePositionAndVolumetricMeasure(Vecd(x, y, 0), particle_spacing_ref * particle_spacing_ref);
                initializeSurfaceProperties(n_0, PT);
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
        if (base_particles_.pos_[index_i][0] < 0.0 || base_particles_.pos_[index_i][1] < 0.0 ||
            base_particles_.pos_[index_i][0] > PL || base_particles_.pos_[index_i][1] > PH)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class TimeStepPartInitialization1 : public TimeStepInitialization
{
  public:
    TimeStepPartInitialization1(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
        : TimeStepInitialization(sph_body, gravity_ptr){};
    virtual ~TimeStepPartInitialization1(){};

  protected:
    SolidParticles *solid_particles = dynamic_cast<SolidParticles *>(particles_);

    virtual void update(size_t index_i, Real dt = 0.0)
    {
        if (solid_particles->pos0_[index_i][0] > 0.0 && solid_particles->pos0_[index_i][1] > 0.0 &&
            solid_particles->pos0_[index_i][0] < PL && solid_particles->pos0_[index_i][1] < PH)
        {
            acc_prior_[index_i] = gravity_->InducedAcceleration(pos_[index_i]);
        }
    }
};
class TimeStepPartInitialization2 : public TimeStepInitialization
{
  public:
    TimeStepPartInitialization2(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
        : TimeStepInitialization(sph_body, gravity_ptr){};
    virtual ~TimeStepPartInitialization2(){};

  protected:
    SolidParticles *solid_particles = dynamic_cast<SolidParticles *>(particles_);
    virtual void update(size_t index_i, Real dt = 0.0)
    {
        if (solid_particles->pos0_[index_i][0] < 0.0 || solid_particles->pos0_[index_i][1] < 0.0 ||
            solid_particles->pos0_[index_i][0] > PL || solid_particles->pos0_[index_i][1] > PH)
        {
            acc_prior_[index_i] = gravity_->InducedAcceleration(pos_[index_i]);
        }
    }
};
/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
  public:
    explicit TimeDependentExternalForce(Vecd external_force)
        : Gravity(external_force) {}
    virtual Vecd InducedAcceleration(Vecd &position) override
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        return current_time < time_to_full_external_force
                   ? current_time * global_acceleration_ / time_to_full_external_force
                   : global_acceleration_;
    }
};
/**
 *  The main program
 */
int main()
{
    /** Setup the system. */
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    system.generate_regression_data_ = false;
    /** create a plate body. */
    SolidBody plate_body(system, makeShared<DefaultShape>("PlateBody"));
    plate_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    plate_body.generateParticles<PlateParticleGenerator>();
    plate_body.addBodyStateForRecording<Vec3d>("PriorAcceleration");

    /** Define Observer. */
    ObserverBody plate_observer(system, "PlateObserver");
    plate_observer.defineParticlesAndMaterial();
    plate_observer.generateParticles<ObserverParticleGenerator>(observation_location);

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation plate_body_inner(plate_body);
    ContactRelation plate_observer_contact(plate_observer, {&plate_body});

    /** Common particle dynamics. */
    SimpleDynamics<TimeStepInitialization> initialize_external_force(plate_body,
                                                                     makeShared<TimeDependentExternalForce>(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration)));

    /**
     * This section define all numerical methods will be used in this case.
     */
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
        corrected_configuration(plate_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
    /** active-passive stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
        stress_relaxation_first_half(plate_body_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
        stress_relaxation_second_half(plate_body_inner);
    /** Constrain the Boundary. */
    BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constrain_holder(boundary_geometry);
    /** Output */
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_plate_max_displacement("Position", io_environment, plate_observer_contact);

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();

    /**
     * From here the time stepping begines.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    write_states.writeToFile(0);
    write_plate_max_displacement.writeToFile(0);

    /** Setup physical parameters. */
    int ite = 0;
    Real end_time = 4.0e-5;
    // Real end_time = 0.8;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integeral_time = 0.0;
        while (integeral_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            initialize_external_force.exec(dt);
            stress_relaxation_first_half.exec(dt);
            constrain_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integeral_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        write_plate_max_displacement.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_plate_max_displacement.generateDataBase(0.005);
    }
    else
    {
        write_plate_max_displacement.testResult();
    }

    return 0;
}
