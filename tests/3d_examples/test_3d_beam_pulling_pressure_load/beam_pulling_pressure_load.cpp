/**
 * @file 	beam_pulling_pressure_load.cpp
 * @brief 	This is the test for comparing SPH with ABAQUS.
 * @author 	Anyong Zhang, Huiqiang Yue
 */

#include "sphinxsys.h"
/** Name space. */
using namespace SPH;

/** Geometry parameters. */
Real resolution_ref = 0.005;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-0.026, -0.026, -0.021), Vecd(0.026, 0.026, 0.101));
StdVec<Vecd> observation_location = {Vecd(0.0, 0.0, 0.04)};

/** Physical parameters */
Real rho = 1265; // kg/m^3
Real poisson_ratio = 0.45;
Real Youngs_modulus = 5e4; // Pa
Real physical_viscosity = 500;

/** Load Parameters */
// Real load_total_force = 12.5; // N
//  Don't be confused with the name of force, here force means pressure.
Real load_total_force = 5000; // pa

/**
 * @brief define the beam body
 */
class Beam : public ComplexShape
{
  public:
    Beam(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        std::string fname_ = "./input/beam.stl";
        Vecd translation(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(fname_, translation, 0.001);
    }
};

/* define load*/
class LoadForce : public BaseLocalDynamics<BodyPartByParticle>, public solid_dynamics::ElasticSolidDataSimple
{
  public:
    LoadForce(BodyPartByParticle &body_part, StdVec<std::array<Real, 2>> f_arr)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          solid_dynamics::ElasticSolidDataSimple(sph_body_),
          acc_prior(particles_->acc_prior_),
          mass_n_(particles_->mass_),
          Vol_(particles_->Vol_),
          F_(particles_->F_),
          force_arr_(f_arr),
          particles_num_(body_part.body_part_particles_.size())
    {
        area_0_.resize(particles_->total_real_particles_);
        for (size_t i = 0; i < particles_->total_real_particles_; ++i)
            area_0_[i] = pow(particles_->Vol_[i], 2.0 / 3.0);
    }

    void update(size_t index_i, Real time = 0.0)
    {
        // pulling direction, i.e. positive z direction
        Vecd normal(0, 0, 1);
        // compute the new normal direction
        const Vecd current_normal = F_[index_i].inverse().transpose() * normal;
        const Real current_normal_norm = current_normal.norm();

        Real J = F_[index_i].determinant();
        // using Nanson’s relation to compute the new area of the surface particle.
        // current_area * current_normal = det(F) * trans(inverse(F)) * area_0 * normal	   =>
        // current_area = J * area_0 * norm(trans(inverse(F)) * normal)   =>
        // current_area = J * area_0 * current_normal_norm
        Real mean_force_ = getForce(time) * J * area_0_[index_i] * current_normal_norm;

        acc_prior[index_i] += (mean_force_ / mass_n_[index_i]) * normal;
    }

  protected:
    StdLargeVec<Vecd> &acc_prior;
    StdLargeVec<Real> &mass_n_;
    StdLargeVec<Real> area_0_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Matd> &F_;

    StdVec<std::array<Real, 2>> force_arr_;
    size_t particles_num_;

  protected:
    virtual Real getForce(Real time)
    {
        for (size_t i = 1; i < force_arr_.size(); i++)
        {
            if (time >= force_arr_[i - 1][0] && time < force_arr_[i][0])
            {
                Real slope = (force_arr_[i][1] - force_arr_[i - 1][1]) / (force_arr_[i][0] - force_arr_[i - 1][0]);
                Real vel = (time - force_arr_[i - 1][0]) * slope + force_arr_[i - 1][1];
                return vel;
            }
            else if (time > force_arr_.back()[0])
                return force_arr_.back()[1];
        }
        return 0.0;
    }
};

/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(system);

    /** Import a beam body, with corresponding material and particles. */
    SolidBody beam_body(system, makeShared<Beam>("beam"));
    beam_body.defineParticlesAndMaterial<ElasticSolidParticles, LinearElasticSolid>(rho, Youngs_modulus, poisson_ratio);
    beam_body.generateParticles<ParticleGeneratorLattice>();

    // Define Observer
    ObserverBody beam_observer(system, "BeamObserver");
    beam_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    /** topology */
    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    /** initialize a time step */
    SimpleDynamics<TimeStepInitialization> beam_initialize(beam_body);

    /** Corrected configuration. */
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(beam_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> update_beam_normal(beam_body);

    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);

    /** specify end-time for defining the force-time profile */
    Real end_time = 1;

    /** === define load === */
    /** create a brick to tag the surface */
    Vecd half_size_0(0.03, 0.03, resolution_ref);
    BodyRegionByParticle load_surface(beam_body, makeShared<TriangleMeshShapeBrick>(half_size_0, 1, Vecd(0.00, 0.00, 0.1)));
    StdVec<std::array<Real, 2>> force_over_time = {
        {Real(0), Real(0)},
        {Real(0.1) * end_time, Real(0.1) * load_total_force},
        {Real(0.4) * end_time, load_total_force},
        {Real(end_time), Real(load_total_force)}};
    SimpleDynamics<LoadForce> pull_force(load_surface, force_over_time);
    std::cout << "load surface particle number: " << load_surface.body_part_particles_.size() << std::endl;

    //=== define constraint ===
    /* create a brick to tag the region */
    Vecd half_size_1(0.03, 0.03, 0.02);
    BodyRegionByParticle holder(beam_body, makeShared<TriangleMeshShapeBrick>(half_size_1, 1, Vecd(0.0, 0.0, -0.02)));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);

    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        beam_damping(0.1, beam_body_inner, "Velocity", physical_viscosity);

    /** Output */
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    RegressionTestTimeAverage<ObservedQuantityRecording<Real>>
        write_beam_stress("VonMisesStress", io_environment, beam_observer_contact);
    /* time step begins */
    GlobalStaticVariables::physical_time_ = 0.0;
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    /** apply initial condition */
    corrected_configuration.exec();
    write_states.writeToFile(0);
    write_beam_stress.writeToFile(0);
    /** Setup physical parameters. */
    int ite = 0;
    Real output_period = end_time / 200.0;
    Real dt = 0.0;

    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }

            beam_initialize.exec();
            pull_force.exec(GlobalStaticVariables::physical_time_);

            /** Stress relaxation and damping. */
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            beam_damping.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = system.getSmallestTimeStepAmongSolidBodies();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        TickCount t2 = TickCount::now();
        write_beam_stress.writeToFile(ite);
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_beam_stress.generateDataBase(0.01, 0.01);
    }
    else
    {
        write_beam_stress.testResult();
    }

    return 0;
}
