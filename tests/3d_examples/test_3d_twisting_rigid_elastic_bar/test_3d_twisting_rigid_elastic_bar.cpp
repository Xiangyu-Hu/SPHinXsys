/**
 * @file twisting_rigid_elastic_bar.cpp
 * @brief This is an example of rigid-flexible coupling simulation.
 * @author Weiyi Kong
 */
#include "sphinxsys.h"
using namespace SPH;

void run_rigid_elastic_coupling(int res_factor = 1);

int main(int ac, char *av[])
{
    run_rigid_elastic_coupling(1);
}

struct solid_algs
{
    InnerRelation inner_relation;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> time_step;

    explicit solid_algs(RealBody &body)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          time_step(body){};

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    double get_time_step() { return time_step.exec(); }
};

BoundingBox union_bounding_box(const BoundingBox &a, const BoundingBox &b)
{
    BoundingBox out = a;
    out.first_[0] = std::min(a.first_[0], b.first_[0]);
    out.first_[1] = std::min(a.first_[1], b.first_[1]);
    out.first_[2] = std::min(a.first_[2], b.first_[2]);
    out.second_[0] = std::max(a.second_[0], b.second_[0]);
    out.second_[1] = std::max(a.second_[1], b.second_[1]);
    out.second_[2] = std::max(a.second_[2], b.second_[2]);
    return out;
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

void run_rigid_elastic_coupling(int res_factor)
{
    // unit transformation
    constexpr Real unit_mm = 1e-3;

    //   elastic   rigid
    //|-----------|----|
    //|           |    |
    //|           |    |
    //|           |    |
    //|-----------|----|
    // geometry
    const Real elastic_length = 6; // length of the elastic part
    const Real rigid_length = 2;   // length of the rigid part
    const Real height = 5;         // height of the bar
    const Real width = 1;          // width of the bar

    // resolution
    const Real dp = width / (4.0 * res_factor);

    // material properties
    const Real rho = 1000 * pow(unit_mm, 2); // density
    const Real youngs_modulus = 5;           // MPa
    const Real poisson_ratio = 0.45;         // Poisson ratio
    auto material = makeShared<SaintVenantKirchhoffSolid>(rho, youngs_modulus, poisson_ratio);

    // load and end time
    const Real end_time = 1.0;
    auto get_force_p = [&](Real t)
    {
        if (t < 0.1)
            return 0.5 * t;
        else if (t < 0.2)
            return 0.5 * (0.2 - t);
        else
            return 0.0;
    };

    // Import meshes
    const Vec3d elastic_halfsize = 0.5 * Vec3d(elastic_length, height, width);
    const Vec3d elastic_translation = elastic_halfsize.x() * Vec3d::UnitX();
    auto mesh_elastic = makeShared<TransformShape<GeometricShapeBox>>(Transform(elastic_translation), elastic_halfsize, "elastic_bar");

    const Vec3d rigid_halfsize = 0.5 * Vec3d(rigid_length, height, width);
    const Vec3d rigid_translation = (elastic_length + 0.5 * rigid_length) * Vec3d::UnitX();
    auto mesh_rigid = makeShared<TransformShape<GeometricShapeBox>>(Transform(rigid_translation), rigid_halfsize, "rigid_bar");

    // System bounding box
    auto bbox = union_bounding_box(mesh_elastic->getBounds(), mesh_rigid->getBounds());

    // System
    SPHSystem system(bbox, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody elastic_body(system, mesh_elastic);
    elastic_body.defineMaterial<SaintVenantKirchhoffSolid>(*material.get());
    elastic_body.generateParticles<BaseParticles, Lattice>();

    SolidBody rigid_body(system, mesh_rigid);
    rigid_body.defineMaterial<Solid>(rho);
    rigid_body.generateParticles<BaseParticles, Lattice>();

    // Initial normal
    SimpleDynamics<NormalDirectionFromBodyShape>{elastic_body}.exec();
    SimpleDynamics<NormalDirectionFromBodyShape>{rigid_body}.exec();

    // Methods
    solid_algs elastic_algs(elastic_body);

    // Boundary conditions
    FixPart fix_elastic_part(elastic_body, "ClampingPart", [&, x0 = mesh_elastic->getBounds().first_.x()](Vec3d &pos)
                             { return pos.x() < x0 + 0.5 * dp; });
    SimpleDynamics<FixBodyPartConstraint> fix_elastic_bc(fix_elastic_part);

    // Multibody system
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the rigid solid box. */
    SolidBodyPartForSimbody rigid_multibody(rigid_body, mesh_rigid);
    SimTK::Body::Rigid rigid_info(*rigid_multibody.body_part_mass_properties_);
    SimTK::MobilizedBody::Free rigidMBody(matter.Ground(), SimTK::Transform(SimTKVec3(0)), rigid_info, SimTK::Transform(SimTKVec3(0)));

    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);

    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    /** Coupling between SimBody and SPH.*/
    // ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_rigid_bar(rigid_multibody, MBsystem, rigidMBody, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_rigid_bar(rigid_multibody, MBsystem, rigidMBody, integ);

    // output
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addToWrite<Vec3d>(elastic_body, "Velocity");
    vtp_output.addToWrite<Vec3d>(rigid_body, "Velocity");
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(elastic_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(rigid_body);
    vtp_output.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    elastic_algs.corrected_config();

    // Simulation
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    const Real dt_ref = elastic_algs.get_time_step();
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

                dt = elastic_algs.get_time_step();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                elastic_algs.stress_relaxation_first(dt);
                fix_elastic_bc.exec();
                elastic_algs.stress_relaxation_second(dt);

                {
                    SimTK::State &state_for_update = integ.updAdvancedState();
                    force_on_bodies.clearAllBodyForces(state_for_update);

                    Real force_p = get_force_p(physical_time);
                    SimTKVec3 force_vec(0, -force_p, 0);           // downward force
                    SimTKVec3 torque_vec(-force_p * height, 0, 0); // torque
                    SimTK::SpatialVec force_on_rigid_body(torque_vec, force_vec);

                    force_on_bodies.setOneBodyForce(state_for_update, rigidMBody, force_on_rigid_body);
                    //   force_on_bodies.setOneBodyForce(state_for_update, rigidMBody, force_on_rigid_bar.exec());
                    integ.stepBy(dt);
                    constraint_rigid_bar.exec();
                }

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
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