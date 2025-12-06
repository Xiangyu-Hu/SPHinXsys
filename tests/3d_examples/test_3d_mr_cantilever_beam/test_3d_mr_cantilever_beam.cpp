/**
 * @file 	test_3d_mr_cantilever_beam.cpp
 * @brief 	This is the case file for the test of multi-resolution solid dynamics.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

struct return_data
{
    Real deflection = std::numeric_limits<Real>::max();
    Real error = std::numeric_limits<Real>::max(); // % error
    size_t number_of_particles = std::numeric_limits<size_t>::max();
    Real dp_max = std::numeric_limits<Real>::max();       // the maximum particle spacing
    Real dp_min = std::numeric_limits<Real>::max();       // the minimum particle spacing
    Real run_time = std::numeric_limits<Real>::max();     // the total simulation time
    Real damping_time = std::numeric_limits<Real>::max(); // time spent on damping
};

return_data beam_multi_resolution(Real dp_factor, bool damping_on, int refinement_level = 1);

TEST(mr_solid, dp_2)
{
    auto print_result = [&](const return_data &data)
    {
        std::cout << data.dp_max << "\t"
                  << data.dp_min << "\t"
                  << data.number_of_particles << "\t"
                  << data.run_time << "\t"
                  << data.damping_time << "\t"
                  << data.deflection << "\t"
                  << data.error << "\n";
    };
    auto data_single = beam_multi_resolution(4, true, 0);
    auto data_mr = beam_multi_resolution(2, true, 1);
    std::cout << "Resolution max [m]" << "\t"
              << "" << "Resolution min [m]" << "\t"
              << "Particles [-]" << "\t"
              << "Run time [s]" << "\t"
              << "Damping time [s]" << "\t"
              << "Deflection [m]" << "\t"
              << "Error [%]" << "\n";
    print_result(data_single);
    print_result(data_mr);

    ASSERT_LE(data_mr.error, 15); // 15% error
    ASSERT_LE(data_mr.run_time, data_single.run_time);
}

int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

//------------------------------------------------------------------------------
void relax_solid(RealBody &body, BaseInnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    body.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
        }
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
}

struct solid_algs
{
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    SimpleDynamics<NormalDirectionFromBodyShape> normal_direction;

    explicit solid_algs(BaseInnerRelation &inner_relation)
        : corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          normal_direction(inner_relation.getSPHBody()) {};

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void normal_update() { normal_direction.exec(); }
};

class FixPart : public BodyPartByParticle
{
  public:
    FixPart(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    bool tagManually(size_t index_i)
    {
        return pos_[index_i].x() < 0;
    };
};

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for Beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}
//------------------------------------------------------------------------------
struct beam_parameters
{
    const Real length = 0.1;
    const Real height = 0.04;
    const Real width = height;

    const Real rho = 1265;
    const Real youngs_modulus = 5e4;
    const Real poisson_ratio = 0.45;
    const Real physical_viscosity = get_physical_viscosity_general(rho, youngs_modulus, width);
    const Real gravity = 9.8;
};
//------------------------------------------------------------------------------
return_data beam_multi_resolution(Real dp_factor, bool damping_on, int refinement_level)
{
    const beam_parameters params;
    return_data data;

    // resolution
    const Real dp_ref = params.height / 4.0;
    const Real dp = dp_ref / dp_factor;
    const Real extension_length = 4 * dp_ref;

    // load
    Gravity gravity(-params.gravity * Vec3d::UnitY());

    // time
    const Real end_time = 1.5;

    // Material models
    auto material = makeShared<NeoHookeanSolid>(params.rho, params.youngs_modulus, params.poisson_ratio);

    // Import meshes
    const Vec3d halfsize = 0.5 * (Vec3d(params.length, params.height, params.width) + extension_length * Vec3d::UnitX());
    const Vec3d translation = 0.5 * (params.length - extension_length) * Vec3d::UnitX();
    auto mesh = makeShared<GeometricShapeBox>(Transform(translation), halfsize, "beam");

    // refinement region
    const Real refinement_region_length = 0.5 * params.length;
    const Vec3d refinement_halfsize = 0.5 * Vec3d(refinement_region_length, params.height, params.width);
    const Vec3d refinement_translation = (params.length - 0.5 * refinement_region_length) * Vec3d::UnitX();
    auto refinement_region = makeShared<GeometricShapeBox>(Transform(refinement_translation), refinement_halfsize);

    // System bounding box
    BoundingBoxd bb_system = mesh->getBounds();

    // System
    SPHSystem system(bb_system, dp);

    // Create objects
    SolidBody beam_body(system, mesh);
    if (refinement_level > 0)
        beam_body.defineAdaptation<AdaptiveWithinShape>(1.15, 1.0, refinement_level);
    beam_body.defineBodyLevelSetShape();
    beam_body.defineMaterial<NeoHookeanSolid>(*material.get());
    if (refinement_level > 0)
        beam_body.generateParticles<BaseParticles, Lattice, AdaptiveByShape>(*refinement_region);
    else
        beam_body.generateParticles<BaseParticles, Lattice>();

    data.number_of_particles = beam_body.getBaseParticles().TotalRealParticles();
    data.dp_max = beam_body.getSPHBodyResolutionRef();
    data.dp_min = beam_body.getSPHAdaptation().MinimumSpacing();

    // Inner relation
    std::unique_ptr<BaseInnerRelation> inner;
    if (refinement_level > 0)
        inner = std::make_unique<AdaptiveInnerRelation>(beam_body);
    else
        inner = std::make_unique<InnerRelation>(beam_body);

    // relax solid
    relax_solid(beam_body, *inner);

    // Methods
    solid_algs algs(*inner);

    // damping
    std::unique_ptr<AdaptiveSplittingInnerRelation> adaptive_inner_split = nullptr;
    if (refinement_level > 0)
        adaptive_inner_split = std::make_unique<AdaptiveSplittingInnerRelation>(beam_body);
    std::unique_ptr<DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>> damping_single_res;
    std::unique_ptr<DampingWithRandomChoice<InteractionAdaptiveSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>> damping_mr_res;
    if (refinement_level > 0) // for multi-resolution, the damping physical viscosity is doubled as a pair of particles only see each other once
        damping_mr_res = std::make_unique<DampingWithRandomChoice<InteractionAdaptiveSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>>(0.2, *adaptive_inner_split, "Velocity", 2 * params.physical_viscosity);
    else
        damping_single_res = std::make_unique<DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>>(0.2, *inner, "Velocity", params.physical_viscosity);
    auto damping_exec = [&](Real dt)
    {
        if (refinement_level > 0)
            damping_mr_res->exec(dt);
        else
            damping_single_res->exec(dt);
    };

    // Boundary conditions
    FixPart fix_bc_part(beam_body);
    SimpleDynamics<FixBodyPartConstraint> fix_bc(fix_bc_part);

    // gravity
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(beam_body, gravity);

    // Output

    beam_body.getBaseParticles().addVariableToWrite<Vec3d>("GravityForce");
    beam_body.getBaseParticles().addVariableToWrite<Vec3d>("NormalDirection");
    beam_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    if (refinement_level > 0)
    {
        beam_body.getBaseParticles().addVariableToWrite<Real>("SmoothingLengthRatio");
        beam_body.getBaseParticles().addVariableToWrite<int>("ParticleMeshLevel");
    }
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(beam_body);
    vtp_output.writeToFile(0);

    // Observer
    const StdVec<Vecd> observation_locations = {params.length * Vec3d::UnitX()};
    ObserverBody observer(system, "MidObserver");
    observer.generateParticles<ObserverParticles>(observation_locations);
    AdaptiveContactRelation observer_contact(observer, {&beam_body});
    ObservedQuantityRecording<Vecd> write_position("Displacement", observer_contact);
    write_position.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs.corrected_config();

    constant_gravity.exec();

    // Simulation
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    TimeInterval time_damping;
    const Real dt_ref = system.getSmallestTimeStepAmongSolidBodies();
    std::cout << "dt_ref: " << dt_ref << std::endl;

    auto run_simulation = [&]()
    {
        while (physical_time < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                dt = system.getSmallestTimeStepAmongSolidBodies();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                algs.stress_relaxation_first(dt);
                fix_bc.exec();
                if (damping_on)
                {
                    TickCount t_damping_0 = TickCount::now();
                    damping_exec(dt);
                    time_damping += (TickCount::now() - t_damping_0);
                    fix_bc.exec();
                }
                algs.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            write_position.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        data.run_time = tt.seconds();
        data.damping_time = time_damping.seconds();
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        vtp_output.writeToFile(ite_output + 1);
        exit(0);
    }

    // error
    // Reference: A dynamic relaxation method with operator splitting and random-choice strategy for SPH
    // https://doi.org/10.1016/j.jcp.2022.111105
    // deflection with the conventional TLSPH method with h/dp=24
    const Real deflection_ref = -0.02 - 5.86e-3;
    data.deflection = write_position.getObservedQuantity()[0].y();
    data.error = std::abs(data.deflection - deflection_ref) / std::abs(deflection_ref) * 100.0;

    return data;
}
