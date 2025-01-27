/**
 * @file twisting_rigid_elastic_bar.cpp
 * @brief This is an example of rigid-flexible coupling simulation.
 * @author Weiyi Kong
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

void run_rigid_elastic_coupling(int res_factor = 1);

int main(int ac, char *av[])
{
    run_rigid_elastic_coupling(4);
}

class FixPart : public BodyPartByParticle
{
  private:
    std::function<bool(Vec3d &)> contains_;

  public:
    FixPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec3d &)> contains)
        : BodyPartByParticle(body, body_part_name),
          contains_(std::move(contains))
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (contains_(pos_[index_i]))
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

Vec3d get_central_position(const BoundingBox &bbox) { return 0.5 * (bbox.first_ + bbox.second_); };

struct solid_algs
{
    InnerRelation inner_relation;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> time_step;

    explicit solid_algs(RealBody &body, Real physical_damping)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          damping(0.5, inner_relation, "Velocity", physical_damping),
          time_step(body)
    {
        // Initial normal
        SimpleDynamics<NormalDirectionFromBodyShape>{body}.exec();
    }

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void damping_exec(Real dt) { damping.exec(dt); }
    Real get_time_step() { return time_step.exec(); }
};

struct rigid_algs
{
    // multi-body system
    SimTK::MultibodySystem MBsystem{};
    SimTK::SimbodyMatterSubsystem matter;
    SimTK::GeneralForceSubsystem forces;

    // rigid body
    SolidBodyPartForSimbody rigid_multibody;
    SimTK::Body::Rigid rigid_info;
    SimTK::MobilizedBody::Free rigidMBody;
    SimTK::Force::DiscreteForces force_on_bodies;

    // state
    SimTK::State state;
    SimTK::RungeKuttaMersonIntegrator integ;

    // coupling between SimBody and SPH
    std::unique_ptr<ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>> force_on_rigid_body;
    std::unique_ptr<SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>> constraint_rigid_body;

    rigid_algs(RealBody &body, SharedPtr<Shape> rigid_shape)
        : matter(MBsystem),
          forces(MBsystem),
          rigid_multibody(body, rigid_shape),
          rigid_info(*rigid_multibody.body_part_mass_properties_),
          rigidMBody(matter.Ground(), SimTK::Transform(EigenToSimTK(get_central_position(rigid_shape->getBounds()))), rigid_info, SimTK::Transform(SimTKVec3(0))),
          force_on_bodies(forces, matter),
          state(MBsystem.realizeTopology()),
          integ(MBsystem)
    {
        integ.setAccuracy(1e-3);
        integ.setAllowInterpolation(false);
        integ.initialize(state);

        force_on_rigid_body = std::make_unique<ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>>(
            rigid_multibody, MBsystem, rigidMBody, integ);
        constraint_rigid_body = std::make_unique<SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>>(
            rigid_multibody, MBsystem, rigidMBody, integ);
    }

    void exec(const SimTK::SpatialVec &ext_force, Real dt)
    {
        SimTK::State &state_for_update = integ.updAdvancedState();
        force_on_bodies.clearAllBodyForces(state_for_update);
        auto coupling_force = force_on_rigid_body->exec();
        force_on_bodies.setOneBodyForce(state_for_update, rigidMBody, ext_force + coupling_force);
        integ.stepBy(dt);
        constraint_rigid_body->exec();
    }

    void apply_constraint() { constraint_rigid_body->exec(); }
};

inline std::vector<Vec3d>
read_ref_data(const fs::path &ref_name)
{
    std::vector<Vec3d> variable_ref;

    std::vector<std::string> row;
    std::string line;
    std::string word;

    std::fstream file(ref_name, std::ios::in);
    if (!file.is_open())
        throw std::runtime_error("Could not open the file: " + ref_name.string());

    double variable_x;
    double variable_y;
    double variable_z;
    while (std::getline(file, line))
    {
        row.clear();
        std::stringstream str(line);
        while (std::getline(str, word, ','))
        {
            row.emplace_back(word);
        }
        // row[0] is node id, skip
        std::stringstream(row[1]) >> variable_x;
        std::stringstream(row[2]) >> variable_y;
        std::stringstream(row[3]) >> variable_z;
        // Relative position y/radius
        variable_ref.emplace_back(variable_x, variable_y, variable_z);
    }

    return variable_ref;
}

void run_rigid_elastic_coupling(int res_factor)
{
    // unit transformation
    constexpr Real unit_mm = 1e-3;

    //   elastic   rigid
    //|-----------|----|
    //|-----------|----|
    // geometry
    const Real elastic_length = 4; // length of the elastic part
    const Real rigid_length = 1;   // length of the rigid part
    const Real total_length = elastic_length + rigid_length;
    const Real height = 1; // height of the bar
    const Real width = 1;  // width of the bar
    const Real x0 = -2;    // the x-coordinate of the left end of the bar

    // resolution
    const Real dp = width / (4.0 * res_factor);
    const Real constraint_length = dp;
    const Real min_x_pos = x0 - dp;

    // material properties
    const Real rho = 1000 * pow(unit_mm, 2); // density
    const Real youngs_modulus = 5;           // MPa
    const Real poisson_ratio = 0.45;         // Poisson ratio
    auto material = makeShared<NeoHookeanSolid>(rho, youngs_modulus, poisson_ratio);

    // load and end time
    const Real end_time = 5.0;
    auto get_force_p = [&](Real t)
    {
        return t < 1.0 ? 0.05 * t : 0.05;
    };
    auto get_rigid_force = [&](Real time)
    {
        Real force_p = get_force_p(time);
        SimTKVec3 force_vec(0, -force_p, 0);           // downward force
        SimTKVec3 torque_vec(-force_p * height, 0, 0); // torque
        return SimTK::SpatialVec(torque_vec, force_vec);
    };

    // Import meshes
    // mesh of the total bar
    const Vec3d halfsize = 0.5 * Vec3d(total_length + constraint_length, height, width);
    const Vec3d translation = (min_x_pos + halfsize.x()) * Vec3d::UnitX() + 0.5 * width * Vec3d::UnitZ();
    auto mesh = makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "bar");

    // mesh of the rigid body part
    const Vec3d rigid_halfsize = 0.5 * Vec3d(rigid_length, height, width);
    const Vec3d rigid_translation = (x0 + elastic_length + 0.5 * rigid_length) * Vec3d::UnitX() + 0.5 * width * Vec3d::UnitZ();
    auto mesh_rigid = makeShared<TransformShape<GeometricShapeBox>>(Transform(rigid_translation), rigid_halfsize, "rigid_bar");

    // System bounding box
    auto bbox = mesh->getBounds();

    // System
    SPHSystem system(bbox, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody body(system, mesh);
    body.defineMaterial<NeoHookeanSolid>(*material.get());
    body.generateParticles<BaseParticles, Lattice>();

    // Methods
    solid_algs algs(body, get_physical_viscosity_general(rho, youngs_modulus, total_length));
    rigid_algs rigid_algs(body, mesh_rigid);

    // Boundary conditions
    FixPart fix_elastic_part(body, "ClampingPart", [&](Vec3d &pos)
                             { return pos.x() < x0; });
    SimpleDynamics<FixBodyPartConstraint> fix_elastic_bc(fix_elastic_part);

    // record rigid particle id for debug
    body.getBaseParticles().registerStateVariable<int>("isRigid", [&](size_t i)
                                                       { return mesh_rigid->checkContain(body.getBaseParticles().ParticlePositions()[i]); });

    // output
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addToWrite<Vec3d>(body, "Velocity");
    vtp_output.addToWrite<Mat3d>(body, "DeformationGradient");
    vtp_output.addToWrite<Mat3d>(body, "StressPK1OnParticle");
    vtp_output.addToWrite<Vec3d>(body, "ForcePrior");
    vtp_output.addToWrite<Vec3d>(body, "Force");
    vtp_output.addToWrite<int>(body, "isRigid");
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(body);
    vtp_output.writeToFile(0);

    // Observer
    const auto observation_locations = read_ref_data("./input/initial_position");
    ObserverBody observer(system, "InterfaceObserver");
    observer.generateParticles<ObserverParticles>(observation_locations);
    ContactRelation contact_observer(observer, {&body});
    InteractionDynamics<CorrectInterpolationKernelWeights>{contact_observer}.exec();
    ObservedQuantityRecording<Vecd> write_disp("Displacement", contact_observer);
    write_disp.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    algs.corrected_config();

    // Simulation
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    const Real dt_ref = algs.get_time_step();
    Real dt = 0.0;
    TickCount t1 = TickCount::now();

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

                dt = algs.get_time_step();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // first half integration of elastic body
                algs.stress_relaxation_first(dt);
                fix_elastic_bc.exec();
                algs.damping_exec(dt);
                fix_elastic_bc.exec();

                // rigid_algs.apply_constraint();

                // update rigid body state based on the inner force computed in the previous step
                auto rigid_force = get_rigid_force(physical_time);
                rigid_algs.exec(rigid_force, dt);

                // second half integration of elastic body
                algs.stress_relaxation_second(dt);

                // rigid_algs.apply_constraint();

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            write_disp.writeToFile(ite_output);
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

    // post-processing
    auto *disp = write_disp.getObservedQuantity();
    const auto disp_ref = read_ref_data("./input/displacement");
    std::cout << "Error: " << std::endl;
    for (size_t i = 0; i < observation_locations.size(); i++)
    {
        std::cout << (disp[i] - disp_ref[i]).norm() / disp_ref[i].norm() * 100 << "%" << std::endl;
    }
}