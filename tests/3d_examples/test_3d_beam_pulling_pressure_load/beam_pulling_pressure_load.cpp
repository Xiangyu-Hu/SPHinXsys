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
BoundingBoxd system_domain_bounds(Vecd(-0.026, -0.026, -0.021), Vecd(0.026, 0.026, 0.101));
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
class PullingForce : public solid_dynamics::BaseLoadingForce<BodyPartByParticle>
{
  public:
    PullingForce(BodyPartByParticle &body_part, StdVec<std::array<Real, 2>> f_arr)
        : solid_dynamics::BaseLoadingForce<BodyPartByParticle>(body_part, "PullingForce"),
          mass_n_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          F_(particles_->getVariableDataByName<Matd>("DeformationGradient")),
          force_arr_(f_arr),
          particles_num_(body_part.body_part_particles_.size())
    {
        area_0_.resize(particles_->TotalRealParticles());
        for (size_t i = 0; i < particles_->TotalRealParticles(); ++i)
            area_0_[i] = pow(Vol_[i], 2.0 / 3.0);
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

        loading_force_[index_i] = mean_force_ * normal;
        solid_dynamics::BaseLoadingForce<BodyPartByParticle>::update(index_i, time);
    }

  protected:
    Real *mass_n_;
    StdVec<Real> area_0_;
    Real *Vol_;
    Matd *F_;

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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif

    /** Import a beam body, with corresponding material and particles. */
    SolidBody beam_body(sph_system, makeShared<Beam>("beam"));
    beam_body.defineMaterial<LinearElasticSolid>(rho, Youngs_modulus, poisson_ratio);
    beam_body.generateParticles<BaseParticles, Lattice>();

    // Define Observer
    ObserverBody beam_observer(sph_system, "BeamObserver");
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    /** topology */
    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});

    /** Corrected configuration. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(beam_body_inner);

    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(beam_body);

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
    SimpleDynamics<PullingForce> pull_force(load_surface, force_over_time);
    std::cout << "load surface particle number: " << load_surface.body_part_particles_.size() << std::endl;

    //=== define constraint ===
    /* create a brick to tag the region */
    Vecd half_size_1(0.03, 0.03, 0.02);
    BodyRegionByParticle holder(beam_body, makeShared<TriangleMeshShapeBrick>(half_size_1, 1, Vecd(0.0, 0.0, -0.02)));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);

    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        beam_damping(0.1, beam_body_inner, "Velocity", physical_viscosity);

    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(beam_body);
    RegressionTestTimeAverage<ObservedQuantityRecording<Real>>
        write_beam_stress("VonMisesStress", beam_observer_contact);
    /* time step begins */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    /** apply initial condition */
    corrected_configuration.exec();
    write_states.writeToFile(0);
    write_beam_stress.writeToFile(0);
    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real output_period = end_time / 200.0;
    Real dt = 0.0;

    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }

            pull_force.exec(physical_time);

            /** Stress relaxation and damping. */
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            beam_damping.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = sph_system.getSmallestTimeStepAmongSolidBodies();
            integration_time += dt;
            physical_time += dt;
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

    if (sph_system.GenerateRegressionData())
    {
        write_beam_stress.generateDataBase(0.01, 0.01);
    }
    else
    {
        write_beam_stress.testResult();
    }

    return 0;
}
