/**
 * @file 	stretching.cpp
 * @brief   Plane strain necking of a bar.
 * @author 	Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.cma.2013.09.024
 */
#include "sphinxsys.h"

using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.05334;  // beam length
Real PH = 0.012826; // for thick plate;

// reference particle spacing
Real global_resolution = PH / 30;
Real BW = global_resolution * 4.0; // boundary width, at least three particles

/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-PL / 2.0, -PL / 2.0),
                                 Vec2d(2.0 * PL, PL / 2.0));
// two dimensional should be circle smooth between two parts.
//----------------------------------------------------------------------
Real rho0_s = 7850.0;           /**< Reference density. */
Real Shear_modulus = 80.1938e9; /**< Poisson ratio. */
Real Bulk_modulus = 164.21e9;

Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /**< Poisson ratio. */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);

Real yield_stress = 0.45e9;
Real hardening_modulus = 1.2924e8;
Real saturation_flow_stress = 7.15e8;
Real saturation_exponent = 16.93;

Real physical_viscosity = 1.0e4;
Real refer_energy = 0.5 * 8000 * 0.01; // 40

Vecd norm_(1.0, 0.0);
Vecd upper_face_point_(0.02 + 3.0 * global_resolution, 0.0);
Vecd lower_face_point_(0.02, 0.0);

Vecd norm_4(1.0, 0.0);
Vecd upper_face_point_4(0.04 + 3.0 * global_resolution, 0.0);
Vecd lower_face_point_4(0.04, 0.0);

// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL / 2.0, PH / 2.0 - PH * 0.01)};

//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------

std::vector<Vecd> beam_left_stretch_shape{
    Vecd(-BW, -PH / 2), Vecd(-BW, PH / 2), Vecd(0.0, PH / 2),
    Vecd(0.0, -PH / 2), Vecd(-BW, -PH / 2)};

// a beam shape
std::vector<Vecd> beam_shape{
    Vecd(0.0, -PH / 2), Vecd(0.0, PH / 2),
    Vecd(PL / 2.0, PH / 2 - PH * 0.01),
    Vecd(PL, PH / 2), Vecd(PL, -PH / 2),
    Vecd(PL / 2.0, -PH / 2 + PH * 0.01),
    Vecd(0.0, -PH / 2)};

std::vector<Vecd> beam_right_stretch_shape{
    Vecd(PL, -PH / 2), Vecd(PL, PH / 2), Vecd(PL + BW, PH / 2),
    Vecd(PL + BW, -PH / 2), Vecd(PL, -PH / 2)};

//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(beam_right_stretch_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(beam_left_stretch_shape, ShapeBooleanOps::add);
    }
};

class LeftStretchSolidBodyRegion : public BodyPartMotionConstraint
{
  public:
    // TODO: use only body part as argment since body can be referred from it already
    LeftStretchSolidBodyRegion(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")) {};

    virtual ~LeftStretchSolidBodyRegion() {};

  protected:
    Vecd *vel_;
    Vecd *pos_;
    virtual void update(size_t index_i, Real Dt = 0.0)
    {
        pos_[index_i][0] -= 0.5e-4 * Dt;
    };
};

class RightStretchSolidBodyRegion : public BodyPartMotionConstraint
{
  public:
    // TODO: use only body part as argment since body can be referred from it already
    RightStretchSolidBodyRegion(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")) {};

    virtual ~RightStretchSolidBodyRegion() {};

  protected:
    Vecd *vel_;
    Vecd *pos_;
    virtual void update(size_t index_i, Real Dt = 0.0)
    {
        pos_[index_i][0] += 0.5e-4 * Dt;
    };
};

MultiPolygon createBeamRightStretchShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(beam_right_stretch_shape, ShapeBooleanOps::add);
    return multi_polygon;
};

MultiPolygon createBeamLeftStretchShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(beam_left_stretch_shape, ShapeBooleanOps::add);
    return multi_polygon;
};

MultiPolygon createConstrainBeamShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(beam_left_stretch_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(beam_right_stretch_shape, ShapeBooleanOps::add);

    return multi_polygon;
};

class ConstrainXVelocity : public BodyPartMotionConstraint
{
  public:
    // TODO: use only body part as argment since body can be referred from it already
    ConstrainXVelocity(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")), pos_(particles_->getVariableDataByName<Vecd>("Position")) {};

    virtual ~ConstrainXVelocity() {};

  protected:
    Vecd *vel_;
    Vecd *pos_;
    virtual void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = Vecd(0.0, vel_[index_i][1]);
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
    SPHSystem system(system_domain_bounds, global_resolution);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    system.setReloadParticles(true);
    system.setGenerateRegressionData(false);
    system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody beam_body(system, makeShared<Beam>("StretchingBody"));
    beam_body.defineBodyLevelSetShape();
    beam_body.defineMaterial<NonLinearHardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus, saturation_flow_stress, saturation_exponent);

    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? beam_body.generateParticles<BaseParticles, Reload>(beam_body.getName())
        : beam_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody beam_observer(system, "BeamObserver");
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation beam_body_inner(beam_body);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> beam_body_random_particles(beam_body);
        RelaxationStepInner beam_body_relaxation_step_inner(beam_body_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_ball_state(system);
        ReloadParticleIO write_particle_reload_files(beam_body);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        beam_body_random_particles.exec(0.25);
        write_ball_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            beam_body_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_ball_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }

    InnerRelation beam_body_inner(beam_body);
    ContactRelation beam_observer_contact(beam_observer, {&beam_body});
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> beam_corrected_configuration(beam_body_inner);

    // stress relaxation for the beam
    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf> stress_relaxation_first_half(beam_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);
    ReduceDynamics<TotalKineticEnergy> get_kinetic_energy(beam_body);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(beam_body);
    BodyRegionByParticle beam_left_stretch(beam_body, makeShared<MultiPolygonShape>(createBeamLeftStretchShape()));
    SimpleDynamics<LeftStretchSolidBodyRegion> stretch_beam_left_end(beam_left_stretch);
    BodyRegionByParticle beam_right_stretch(beam_body, makeShared<MultiPolygonShape>(createBeamRightStretchShape()));
    SimpleDynamics<RightStretchSolidBodyRegion> stretch_beam_right_end(beam_right_stretch);
    BodyRegionByParticle beam_constrain(beam_body, makeShared<MultiPolygonShape>(createConstrainBeamShape()));
    SimpleDynamics<ConstrainXVelocity> constrain_beam_end(beam_constrain);
    InteractionDynamics<solid_dynamics::DeformationGradientBySummation> beam_deformation_gradient_tensor(beam_body_inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> damping(0.5, beam_body_inner, "Velocity", physical_viscosity);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_beam_states(system);
    ReducedQuantityRecording<TotalKineticEnergy> write_kinetic_mechanical_energy(beam_body);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_displacement("Position", beam_observer_contact);
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    beam_corrected_configuration.exec();

    //----------------------------------------------------------------------
    //	Setup computing time-step controls.
    //----------------------------------------------------------------------
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int Dt_ite = 0;
    Real End_Time = 100.0;

    Real D_Time = End_Time / 100.0; // time step size for output file
    Real Dt = End_Time / 10000.0;   /**< Time period for stretching */
    Real dt = 0.0;                  // default acoustic time step sizes
    int observation_sample_interval = 1000.0;

    // statistics for computing time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    write_beam_states.writeToFile(0);
    write_kinetic_mechanical_energy.writeToFile(0);
    write_displacement.writeToFile(0);

    // computation loop starts
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < D_Time)
        {
            Real relaxation_time = 0.0;
            stretch_beam_left_end.exec(Dt);
            stretch_beam_right_end.exec(Dt);
            beam_deformation_gradient_tensor.exec(Dt);
            Dt_ite++;

            int stress_ite = 0;
            Real refer_total_kinetic_energy = 10000.0;
            while (relaxation_time < Dt)
            {
                if (refer_total_kinetic_energy > 0.005)
                {
                    stress_relaxation_first_half.exec(dt);
                    constrain_beam_end.exec(dt);
                    damping.exec(Dt);
                    constrain_beam_end.exec(dt);
                    stress_relaxation_second_half.exec(dt);

                    refer_total_kinetic_energy = get_kinetic_energy.exec() / refer_energy;

                    ite++;
                    stress_ite++;

                    dt = computing_time_step_size.exec();
                    if (ite % 500 == 0)
                    {
                        std::cout << "N=" << ite << " Time: "
                                  << physical_time
                                  << "	Dt: " << Dt << "	dt: " << dt
                                  << "	Dt:dt = " << Dt / dt << "\n";
                    }

                    if (ite != 0 && ite % observation_sample_interval == 0)
                    {
                        write_kinetic_mechanical_energy.writeToFile(ite);
                        write_displacement.writeToFile(ite);
                    }
                }
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
        }

        TickCount t2 = TickCount::now();
        write_beam_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds."
              << "  Iterations:  " << ite << std::endl;
    std::cout << "Total iterations computation:  " << physical_time / dt
              << std::endl;

    if (system.GenerateRegressionData())
    {
        write_displacement.generateDataBase(0.05);
    }
    else
    {
        write_displacement.testResult();
    }

    return 0;
}
