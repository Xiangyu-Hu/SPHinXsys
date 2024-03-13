/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: Spinning Cube                         *
 *-----------------------------------------------------------------------------*
 * This is a basic test to examine if the viscous operator is capable of       *
 * accurately capturing angular momentum conservation                          *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup properties
Real particle_spacing = 0.1;
Real gravity_g = 0.0;
Real end_time = 13.1;
int nmbr_of_outputs = 100;
bool linearized_iteration = false;

// material properties
Real rho = 1.0;                          // reference density
Real omega = -1.0;                       // lid velocity
Real SOS = 10.0 * (Real)std::abs(omega); // numerical speed of sound

// non-Newtonian properties
Real K = 100;   // consistency index
Real n = 1;     // power index
Real tau_y = 0; // yield stress

Real min_shear_rate = 1e-2; // cutoff low shear rate
Real max_shear_rate = 1e+2; // cutoff high shear rate

// geometry data
Real length = 1;

//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class FluidFilling : public ComplexShape
{
  public:
    explicit FluidFilling(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd cube(0.5 * length, 0.5 * length, 0.5 * length);
        Vecd translation(0.0, 0.0, 0.0);
        Transform translate_to_point(translation);
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_point), cube);
    }
};

class InitialVelocity
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    InitialVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          fluid_particles_(dynamic_cast<BaseParticles *>(&sph_body.getBaseParticles())){};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        Matd rotation;
        rotation << 0, 1, 0,
            -1, 0, 0,
            0, 0, 0;
        vel_[index_i] = omega * rotation * pos_[index_i];
    }

  protected:
    BaseParticles *fluid_particles_;
};

void output_setup()
{
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "XXXXXXXX Spinning Cube Case XXXXXXXX" << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    std::cout << "     particle_spacing= " << particle_spacing << std::endl;
    std::cout << "     end_time= " << end_time << std::endl;
    std::cout << "     K= " << K << std::endl;
    std::cout << "     n= " << n << std::endl;
    std::cout << "     tau_y= " << tau_y << std::endl;
    std::cout << "     min_shear_rate= " << min_shear_rate << std::endl;
    std::cout << "     max_shear_rate= " << max_shear_rate << std::endl;
    std::cout << "     rho= " << rho << std::endl;
    std::cout << "     SOS= " << SOS << std::endl;

    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
}

int main(int ac, char *av[])
{
    output_setup();
    //	Build up an SPHSystem
    BoundingBox system_domain_bounds(Vecd(-1, -1, -1), Vecd(1, 1, 1));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<FluidFilling>("FluidBody"));
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, min_shear_rate, max_shear_rate, K, n, tau_y);
    fluid.generateParticles<ParticleGeneratorLattice>();

    //	Define body relation map
    InnerRelation fluid_inner(fluid);

    //	Define the numerical methods used in the simulation
    SimpleDynamics<InitialVelocity> initial_condition(fluid);

    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(fluid_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> density_relaxation(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(fluid_inner);

    InteractionDynamics<fluid_dynamics::VelocityGradient<Inner<>>> vel_grad_calculation(fluid_inner);
    InteractionDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_rate_calculation(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::GeneralizedNewtonianViscousForce<Inner<>>> viscous_acceleration(fluid_inner);

    // InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<ZerothInconsistencyLimiter, BulkParticles>> transport_velocity_correction(fluid_inner);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid, 10.0 * omega);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_acoustic_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    //	Define the methods for I/O operations, observations
    fluid.addBodyStateForRecording<Real>("Pressure");
    fluid.addBodyStateForRecording<Real>("Density");
    fluid.addBodyStateForRecording<Real>("Mass");
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);

    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();

    //	Setup for time-stepping control
    // size_t number_of_iterations = sph_system.RestartStep();
    Real output_interval = end_time / nmbr_of_outputs;
    Real dt = 0;
    Real Dt = 0;
    Real Dt_visc = 0;
    Real Dt_adv = 0;
    Real Dt_aco = 0;
    int iteration = 0;
    int output_counter = 1;

    //	First output before the main loop.
    write_fluid_states.writeToFile(0);
    TickCount t1 = TickCount::now();
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        iteration++;
        TimeInterval tt;
        TickCount t2 = TickCount::now();
        tt = t2 - t1;
        Dt_adv = get_fluid_advection_time_step_size.exec();
        Dt_visc = get_viscous_time_step_size.exec();

        if (linearized_iteration == true && Dt_visc < Dt_adv)
        {
            Real viscous_time = 0.0;
            update_density_by_summation.exec(Dt_adv);
            vel_grad_calculation.exec(Dt_adv);
            shear_rate_calculation.exec(Dt_adv);

            while (viscous_time < Dt_adv)
            {
                viscous_acceleration.exec(Dt_visc);
                viscous_time += Dt_visc;
                if (viscous_time + Dt_visc > Dt_adv)
                {
                    Dt_visc = Dt_adv - viscous_time;
                }
            }
        }
        else
        {
            Dt = SMIN(Dt_visc, Dt_adv);
            update_density_by_summation.exec(Dt);
            vel_grad_calculation.exec(Dt);
            shear_rate_calculation.exec(Dt);
            viscous_acceleration.exec(Dt);
        }

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            Dt_aco = get_acoustic_time_step_size.exec();
            dt = SMIN(Dt_aco, Dt);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            relaxation_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }

        if (iteration < 100 || output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            std::cout << std::fixed << std::setprecision(2) << std::scientific << "Iteration: " << iteration << " | sim time in %: " << GlobalStaticVariables::physical_time_ / end_time * 100 << " | physical time in s: " << GlobalStaticVariables::physical_time_ << " | computation time in s: " << tt.seconds() << " | dt_adv: " << Dt_adv << " | dt_visc: " << Dt_visc << " | dt_aco: " << Dt_aco << "\r" << std::flush;
        }

        if (output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            write_fluid_states.writeToFile();
            output_counter++;
        }
        fluid.updateCellLinkedListWithParticleSort(100);
        fluid_inner.updateConfiguration();
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    std::cout << "Dt " << Dt << " | Dt_adv " << Dt_adv << " | Dt_visc " << Dt_visc << " | dt " << dt << std::endl;
    return 0;
}