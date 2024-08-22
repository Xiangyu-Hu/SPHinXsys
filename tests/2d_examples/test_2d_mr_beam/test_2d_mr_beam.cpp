/**
 * @file 	Three ring impact.cpp
 * @brief 	This is the case file for the test of dynamic contacts between shell and solid.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"
using namespace SPH;

void beam_single_resolution(int dp_factor);
void beam_multi_resolution(int dp_factor, int refinement_level = 1);

int main(int ac, char *av[])
{
    beam_multi_resolution(1, 2);
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
    SimpleDynamics<UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(body);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    update_smoothing_length_ratio.exec();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        update_smoothing_length_ratio.exec();
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 200 == 0)
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
}

struct solid_algs
{
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size;
    SimpleDynamics<NormalDirectionFromBodyShape> normal_direction;

    explicit solid_algs(BaseInnerRelation &inner_relation)
        : corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          computing_time_step_size(inner_relation.getSPHBody()),
          normal_direction(inner_relation.getSPHBody()){};

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void normal_update() { normal_direction.exec(); }
    Real time_step_size() { return computing_time_step_size.exec(); }
};

class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name, const Vec2d &center, Real length, Real height, Real extension_length) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(center), 0.5 * Vec2d(length + extension_length, height), ShapeBooleanOps::add);
    }
};

class FixPart : public BodyPartByParticle
{
  private:
    const Vec2d &center_;
    Real length_;

  public:
    FixPart(SPHBody &body, const std::string &body_part_name, const Vec2d &center, Real length)
        : BodyPartByParticle(body, body_part_name),
          center_(center), length_(length)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (pos_[index_i].x() > center_.x() + 0.5 * length_ || pos_[index_i].x() < center_.x() - 0.5 * length_)
            body_part_particles_.push_back(index_i);
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
    const Real length = 2;
    const Real height = 0.05;
    const Real rho = 2700;
    const Real youngs_modulus = 6.75e6;
    const Real poisson_ratio = 0.34;
    const Real physical_viscosity = get_physical_viscosity_general(rho, youngs_modulus, height);
    const Real gravity = 50;
};
//------------------------------------------------------------------------------
void beam_single_resolution(int dp_factor)
{
    const beam_parameters params;

    // resolution
    const Real dp = params.height / (4.0 * dp_factor);
    const Real extension_length = 4 * dp;

    // load
    Gravity gravity(-params.gravity * Vec2d::UnitY());

    // time
    const Real end_time = 5;

    // Material models
    auto material = makeShared<SaintVenantKirchhoffSolid>(params.rho, params.youngs_modulus, params.poisson_ratio);

    // Import meshes
    const auto translation = Vec2d::Zero();
    auto mesh = std::make_shared<Beam>("Beam", Vec2d::Zero(), params.length, params.height, extension_length);

    // System bounding box
    BoundingBox bb_system = mesh->getBounds();

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody beam_body(system, mesh);
    beam_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    beam_body.assignMaterial(material.get());
    beam_body.generateParticles<BaseParticles, Lattice>();

    // Inner relation
    InnerRelation inner(beam_body);

    // relax solid
    relax_solid(beam_body, inner);

    // Methods
    solid_algs algs(inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> damping(0.5, inner, "Velocity", params.physical_viscosity);

    // Boundary conditions
    FixPart fix_bc_part(beam_body, "ClampingPart", translation, params.length);
    SimpleDynamics<FixBodyPartConstraint> fix_bc(fix_bc_part);

    // gravity
    SimpleDynamics<GravityForce> constant_gravity(beam_body, gravity);

    // Check
    auto check_nan = [&](BaseParticles &particles)
    {
        const auto &pos_ = particles.ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("GravityForce");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("NormalDirection");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("Position");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Observer
    const StdVec<Vecd> observation_locations = {Vec2d::Zero()};
    ObserverBody observer(system, "MidObserver");
    observer.generateParticles<ObserverParticles>(observation_locations);
    ContactRelation observer_contact(observer, {&beam_body});
    ObservedQuantityRecording<Vecd> write_position("Position", observer_contact);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs.corrected_config();

    constant_gravity.exec();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = algs.time_step_size();
    auto run_simulation = [&]()
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                dt = algs.time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                algs.stress_relaxation_first(dt);
                fix_bc.exec();
                damping.exec(dt);
                fix_bc.exec();
                algs.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                // checking if any position has become nan
                check_nan(beam_body.getBaseParticles());
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            write_position.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
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
}
//------------------------------------------------------------------------------
void beam_multi_resolution(int dp_factor, int refinement_level)
{
    const beam_parameters params;

    // resolution
    const Real dp = params.height / (4.0 * dp_factor);
    const Real extension_length = 4 * dp;

    // load
    Gravity gravity(-params.gravity * Vec2d::UnitY());

    // time
    const Real end_time = 5;

    // Material models
    auto material = makeShared<SaintVenantKirchhoffSolid>(params.rho, params.youngs_modulus, params.poisson_ratio);

    // Import meshes
    const auto translation = Vec2d::Zero();
    auto mesh = std::make_shared<Beam>("Beam", Vec2d::Zero(), params.length, params.height, extension_length);

    // refinement region
    MultiPolygonShape refinement_region = [&, refinement_region_length = 0.25 * params.length]()
    {
        MultiPolygon shape;
        shape.addABox(Transform(translation), 0.5 * Vec2d(refinement_region_length, params.height), ShapeBooleanOps::add);
        return MultiPolygonShape(shape, "RefinementRegion");
    }();

    // System bounding box
    BoundingBox bb_system = mesh->getBounds();

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody beam_body(system, mesh);
    beam_body.defineAdaptation<ParticleRefinementWithinShape>(1.15, 1.0, refinement_level);
    beam_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    beam_body.assignMaterial(material.get());
    beam_body.generateParticles<BaseParticles, Lattice, Adaptive>(refinement_region);

    // Inner relation
    AdaptiveInnerRelation inner(beam_body);

    // relax solid
    relax_solid(beam_body, inner);

    // Methods
    solid_algs algs(inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> damping(0.5, inner, "Velocity", params.physical_viscosity);

    // Boundary conditions
    FixPart fix_bc_part(beam_body, "ClampingPart", translation, params.length);
    SimpleDynamics<FixBodyPartConstraint> fix_bc(fix_bc_part);

    // gravity
    SimpleDynamics<GravityForce> constant_gravity(beam_body, gravity);

    // Check
    auto check_nan = [&](BaseParticles &particles)
    {
        const auto &pos_ = particles.ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("GravityForce");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("NormalDirection");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("Position");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Observer
    const StdVec<Vecd> observation_locations = {Vec2d::Zero()};
    ObserverBody observer(system, "MidObserver");
    observer.generateParticles<ObserverParticles>(observation_locations);
    AdaptiveContactRelation observer_contact(observer, {&beam_body});
    ObservedQuantityRecording<Vecd> write_position("Position", observer_contact);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs.corrected_config();

    constant_gravity.exec();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = algs.time_step_size();
    auto run_simulation = [&]()
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                dt = algs.time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                algs.stress_relaxation_first(dt);
                fix_bc.exec();
                damping.exec(dt);
                fix_bc.exec();
                algs.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                // checking if any position has become nan
                check_nan(beam_body.getBaseParticles());
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            write_position.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
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
}