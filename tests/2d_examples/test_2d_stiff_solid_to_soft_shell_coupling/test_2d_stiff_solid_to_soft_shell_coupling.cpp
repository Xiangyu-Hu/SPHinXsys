#include "sphinxsys.h"
using namespace SPH;

void run_solid_to_shell_coupling(size_t res_factor_solid, size_t res_factor_shell, Real stiffness_ratio, bool run_relax);
void run_solid(size_t res_factor, Real stiffness_ratio, bool run_relax);

int main(int ac, char *av[])
{
    // run_solid_to_shell_coupling(1, 1, 1.0, true);
    run_solid(1, 1.0, false);
}

//-------Relaxation for the solid body-------------------------------------------------
inline void relax_solid(BaseInnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(inner.getSPHBody());
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    inner.real_body_->updateCellLinkedList();
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

class SolidMaterialInitialization : public MaterialIdInitialization
{
  private:
    std::function<bool(Vec2d &)> contain_;

  public:
    SolidMaterialInitialization(SolidBody &solid_body, std::function<bool(Vec2d &)> contain)
        : MaterialIdInitialization(solid_body),
          contain_(std::move(contain)) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        if (contain_(pos_[index_i]))
            material_id_[index_i] = 1;
        else
            material_id_[index_i] = 0;
    };
};

struct solid_algs
{
    InnerRelation inner_relation;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    SimpleDynamics<NormalDirectionFromBodyShape> initial_normal_direction;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> normal_direction;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size;

    explicit solid_algs(SolidBody &body, Real physical_viscosity)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          initial_normal_direction(body),
          normal_direction(body),
          damping(0.5, inner_relation, "Velocity", physical_viscosity),
          computing_time_step_size(body)
    {
        initial_normal_direction.exec();
    };

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void normal_update() { normal_direction.exec(); }
    void damping_exec(Real dt) { damping.exec(dt); }
    Real time_step_size() { return computing_time_step_size.exec(); }
};

struct shell_algs
{
    InnerRelation inner_relation;
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half;
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> normal_direction;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> velocity_damping;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> rotation_damping;
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size;

    explicit shell_algs(SolidBody &body, Real physical_viscosity)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          normal_direction(inner_relation.getSPHBody()),
          velocity_damping(0.5, inner_relation, "Velocity", physical_viscosity),
          rotation_damping(0.5, inner_relation, "AngularVelocity", physical_viscosity),
          computing_time_step_size(body) {};

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void normal_update() { normal_direction.exec(); }
    void damping_exec(Real dt)
    {
        velocity_damping.exec(dt);
        rotation_damping.exec(dt);
    }
    Real time_step_size() { return computing_time_step_size.exec(); }
};

class SolidBodyPart : public BodyPartByParticle
{
  private:
    std::function<bool(Vec2d &)> contain_;

  public:
    SolidBodyPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec2d &)> contain)
        : BodyPartByParticle(body, body_part_name),
          contain_(std::move(contain))
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&SolidBodyPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (contain_(pos_[index_i]))
            body_part_particles_.push_back(index_i);
    };
};

namespace SPH
{
class ShellDirectGenerator;
template <>
class ParticleGenerator<SurfaceParticles, ShellDirectGenerator> : public ParticleGenerator<SurfaceParticles>
{
  private:
    Real dp_;
    Real thickness_;
    StdLargeVec<Vecd> pos_;
    StdLargeVec<Vecd> n_;

  public:
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, const StdLargeVec<Vecd> &pos, const StdLargeVec<Vecd> &n, Real dp, Real thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          dp_(dp),
          thickness_(thickness),
          pos_(pos),
          n_(n)
    {
        if (pos_.size() != n_.size())
        {
            std::cout << "\n Error: the size of position and normal vector are not consistent!" << std::endl;
            exit(0);
        }
    };

    void prepareGeometricData() override
    {
        for (size_t i = 0; i < pos_.size(); i++)
        {
            addPositionAndVolumetricMeasure(pos_[i], pow(dp_, Dimensions - 1));
            addSurfaceProperties(n_[i], thickness_);
        }
    }
};
} // namespace SPH

class CouplingPart : public BodyPartByParticle
{
  private:
    int *is_coupled_;
    std::function<bool(Vec2d &)> contain_;

  public:
    CouplingPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec2d &)> contain)
        : BodyPartByParticle(body, body_part_name),
          is_coupled_(base_particles_.registerStateVariable<int>("IsCoupled", 0)),
          contain_(std::move(contain))
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&CouplingPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (contain_(pos_[index_i]))
        {
            is_coupled_[index_i] = 1;
            body_part_particles_.push_back(index_i);
        }
    }
};

struct solid_coupling_algs
{
    NearestNeighborContactRelation contact_relation;
    CouplingPart part;
    InteractionWithUpdate<solid_dynamics::InterpolationForceConstraint> force_bc;

    solid_coupling_algs(SolidBody &body,
                        RealBodyVector contact_bodies,
                        std::function<bool(Vec2d &)> contain,
                        const std::vector<Real> &factors = {})
        : contact_relation(body, std::move(contact_bodies), factors),
          part(body, "CouplingPart", std::move(contain)),
          force_bc(part, contact_relation) {}

    void exec() { force_bc.exec(); }
};

struct shell_coupling_algs
{
    NearestNeighborContactRelation contact_relation;
    CouplingPart part;
    SimpleDynamics<solid_dynamics::TotalWeightComputation> total_weight;
    SimpleDynamics<solid_dynamics::InterpolationVelocityConstraint> vel_bc;

    shell_coupling_algs(SolidBody &body,
                        RealBodyVector contact_bodies,
                        std::function<bool(Vec2d &)> contain,
                        const std::vector<Real> &factors = {})
        : contact_relation(body, std::move(contact_bodies), factors),
          part(body, "CouplingPart", std::move(contain)),
          total_weight(part, contact_relation),
          vel_bc(part, contact_relation) {}

    void initialize_total_weight() { total_weight.exec(); }
    void exec() { vel_bc.exec(); }
};

inline Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

void run_solid_to_shell_coupling(size_t res_factor_solid, size_t res_factor_shell, Real stiffness_ratio, bool run_relax)
{
    // unit
    constexpr Real unit_mm = 1e-3;

    // geometry
    const Real cube_length = 1.0;
    const Real shell_thickness = 0.1;
    const Real shell_length = 5.0;

    const Real dp_ref = cube_length / 10.0;
    const Real dp_solid = dp_ref / Real(res_factor_solid);
    const Real dp_shell = dp_ref / Real(res_factor_shell);

    // import model
    MultiPolygon cube_shape{};
    cube_shape.addABox(Transform((0.5 * cube_length + dp_shell) * Vec2d::UnitY()), 0.5 * cube_length * Vec2d::Ones(), SPH::ShapeBooleanOps::add);
    MultiPolygonShape cube_mesh(cube_shape, "Cube");

    // shell positions
    StdLargeVec<Vecd> shell_pos;
    StdLargeVec<Vecd> shell_n;
    {
        Real x = -0.5 * shell_length + 0.5 * dp_shell;
        Real y = 0.5 * dp_solid;
        while (x < 0.5 * shell_length)
        {
            shell_pos.emplace_back(x, y);
            shell_n.emplace_back(0, 1);
            x += dp_shell;
        }
    }

    // Material
    Real rho = 1000 * pow(unit_mm, 2);
    Real youngs_modulus_shell = 3; // 3 MPa
    Real youngs_modulus_solid = youngs_modulus_shell * stiffness_ratio;
    Real poisson_ratio = 0.45;
    Real physical_viscosity = get_physical_viscosity_general(rho, youngs_modulus_solid, cube_length);
    auto material_solid = makeShared<NeoHookeanSolid>(rho, youngs_modulus_solid, poisson_ratio);
    auto material_shell = makeShared<NeoHookeanSolid>(rho, youngs_modulus_shell, poisson_ratio);

    // System bounding box
    BoundingBox bb_system = cube_mesh.getBounds();

    // System
    SPHSystem system(bb_system, dp_solid);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody cube_body(system, cube_mesh, "Cube");
    cube_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    cube_body.defineMaterial<NeoHookeanSolid>(rho, youngs_modulus_solid, poisson_ratio);
    cube_body.generateParticles<BaseParticles, Lattice>();

    SolidBody shell_body(system, makeShared<DefaultShape>("shell"));
    shell_body.defineAdaptationRatios(1.15, dp_solid / dp_shell);
    shell_body.defineMaterial<NeoHookeanSolid>(rho, youngs_modulus_shell, poisson_ratio);
    shell_body.generateParticles<SurfaceParticles, ShellDirectGenerator>(shell_pos, shell_n, dp_shell, shell_thickness);

    // Algorithm
    solid_algs algs_solid(cube_body, physical_viscosity);
    shell_algs algs_shell(shell_body, physical_viscosity);

    // Relax
    if (run_relax)
        relax_solid(algs_solid.inner_relation);

    // Boundary conditions
    SolidBodyPart shell_fixed_part(shell_body, "ShellFixedPart", [&](Vec2d &pos)
                                   { return pos.x() < -0.5 * shell_length + 0.7 * dp_shell || pos.x() > 0.5 * shell_length - 0.7 * dp_shell; });
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> fix_shell_bc(shell_fixed_part);

    // Gravity
    Gravity gravity(-1 * Vec2d::UnitY());
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(cube_body, gravity);

    // Coupling
    shell_coupling_algs shell_coupling(shell_body, {&cube_body}, [&](Vec2d &pos)
                                       { return pos.x() > bb_system.first_.x() && pos.x() < bb_system.second_.x(); }, {2.3});
    solid_coupling_algs solid_coupling(cube_body, {&shell_body}, [&](Vec2d &pos)
                                       { return true; }, {2.3});

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    algs_solid.corrected_config();
    algs_shell.corrected_config();
    shell_coupling.initialize_total_weight();
    constant_gravity.exec();

    // Output
    cube_body.getBaseParticles().addVariableToWrite<int>("IsCoupled");
    shell_body.getBaseParticles().addVariableToWrite<int>("IsCoupled");
    cube_body.getBaseParticles().addVariableToWrite<Vecd>("CouplingForce");
    cube_body.getBaseParticles().addVariableToWrite<Vecd>("ForcePrior");
    shell_body.getBaseParticles().addVariableToWrite<Real>("TotalWeight");
    cube_body.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    shell_body.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(cube_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell_body);
    vtp_output.writeToFile(0);

    // Output
    body.getBaseParticles().addVariableToWrite<int>("MaterialID");
    body.getBaseParticles().addVariableToWrite<Vecd>("ForcePrior");
    body.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(body);
    vtp_output.writeToFile(0);

    // Observer
    StdVec<Vecd> observation_locations{Vec2d(0, 0.5 * shell_thickness)};
    ObserverBody observer_body(system, "Observer");
    observer_body.generateParticles<ObserverParticles>(observation_locations);
    ContactRelation obs_contact(observer_body, {&shell_body});
    obs_contact.updateConfiguration();
    InteractionDynamics<CorrectInterpolationKernelWeights>{obs_contact}.exec();
    ObservedQuantityRecording<Vecd> disp_recorder("Displacement", obs_contact);
    ObservedQuantityRecording<Real> stress_recorder("VonMisesStress", obs_contact);
    disp_recorder.writeToFile(0);
    stress_recorder.writeToFile(0);

    // Simulation
    const Real end_time = 10.0;
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min(algs_solid.time_step_size(), algs_shell.time_step_size());
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

                dt = std::min(algs_solid.time_step_size(), algs_shell.time_step_size());
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // shell 1st half
                algs_shell.stress_relaxation_first(dt);

                // compute solid coupling force
                solid_coupling.exec();

                // update solid
                algs_solid.stress_relaxation_first(dt);
                algs_solid.damping_exec(dt);
                algs_solid.stress_relaxation_second(dt);

                // update shell kinematic constraint and 2nd half
                shell_coupling.exec();
                fix_shell_bc.exec();
                algs_shell.damping_exec(dt);
                shell_coupling.exec();
                fix_shell_bc.exec();
                algs_shell.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            disp_recorder.writeToFile(ite_output);
            stress_recorder.writeToFile(ite_output);
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

    // Output results
    std::cout << "Deflection: " << disp_recorder.getObservedQuantity()[0].y() << std::endl;
    std::cout << "Von Mises Stress: " << stress_recorder.getObservedQuantity()[0] << std::endl;
}

class ForcePartByParticle : public BaseForcePrior<BodyPartByParticle>
{
  private:
    Vecd acc_;
    Real *mass_;

  public:
    ForcePartByParticle(BodyPartByParticle &body_part, const std::string &force_name, const Vecd &acc)
        : BaseForcePrior<BodyPartByParticle>(body_part, force_name),
          acc_(acc),
          mass_(particles_->registerStateVariable<Real>("Mass")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        current_force_[index_i] = mass_[index_i] * acc_;
        BaseForcePrior<BodyPartByParticle>::update(index_i, dt);
    };
};

void run_solid(size_t res_factor, Real stiffness_ratio, bool run_relax)
{
    // unit
    constexpr Real unit_mm = 1e-3;

    // geometry
    const Real cube_length = 1.0;
    const Real shell_thickness = 0.1;
    const Real shell_length = 5.0;

    const Real dp = shell_thickness / (4.0 * res_factor);

    // import model
    MultiPolygon shape{};
    shape.addABox(Transform((0.5 * cube_length + shell_thickness) * Vec2d::UnitY()), 0.5 * cube_length * Vec2d::Ones(), SPH::ShapeBooleanOps::add);
    shape.addABox(Transform(0.5 * shell_thickness * Vec2d::UnitY()), 0.5 * Vec2d(shell_length, shell_thickness), SPH::ShapeBooleanOps::add);
    MultiPolygonShape mesh(shape, "Solid");

    // Material
    Real rho = 1000 * pow(unit_mm, 2);
    Real youngs_modulus_shell = 3; // 3 MPa
    Real youngs_modulus_solid = youngs_modulus_shell * stiffness_ratio;
    Real poisson_ratio = 0.45;
    Real physical_viscosity = get_physical_viscosity_general(rho, youngs_modulus_solid, cube_length);

    // System bounding box
    BoundingBox bb_system = mesh.getBounds();

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody body(system, mesh, "solid");
    body.defineBodyLevelSetShape()->cleanLevelSet(0);
    body.defineMaterial<CompositeSolid>(rho);
    body.generateParticles<BaseParticles, Lattice>();
    auto &material = *dynamic_cast<CompositeSolid *>(&body.getBaseMaterial());
    material.add<NeoHookeanSolid>(rho, youngs_modulus_solid, poisson_ratio);
    material.add<NeoHookeanSolid>(rho, youngs_modulus_shell, poisson_ratio);

    // Algorithm
    solid_algs algs_solid(body, physical_viscosity);

    // Relax
    if (run_relax)
        relax_solid(algs_solid.inner_relation);

    // assign material ids
    SimpleDynamics<SolidMaterialInitialization> material_id_initialization(body, [&](Vec2d &pos) -> bool
                                                                           { return pos.y() < shell_thickness; });
    material_id_initialization.exec();

    // Boundary conditions
    SolidBodyPart shell_fixed_part(body, "ShellFixedPart", [&](Vec2d &pos)
                                   { return pos.x() < -0.5 * shell_length + 0.7 * dp || pos.x() > 0.5 * shell_length - 0.7 * dp; });
    SimpleDynamics<FixBodyPartConstraint> fix_shell_bc(shell_fixed_part);

    // Gravity
    auto gravity(-1 * Vec2d::UnitY());
    SolidBodyPart cube_part(body, "CubePart", [&](Vec2d &pos)
                            { return pos.y() > shell_thickness; });
    SimpleDynamics<ForcePartByParticle> constant_gravity(cube_part, "Gravity", gravity);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    algs_solid.corrected_config();
    constant_gravity.exec();

    // Output
    body.getBaseParticles().addVariableToWrite<int>("MaterialID");
    body.getBaseParticles().addVariableToWrite<Vecd>("ForcePrior");
    body.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(body);
    vtp_output.writeToFile(0);

    // Observer
    StdVec<Vecd> observation_locations{Vec2d(0, 0.5 * shell_thickness)};
    ObserverBody observer_body(system, "Observer");
    observer_body.generateParticles<ObserverParticles>(observation_locations);
    ContactRelation obs_contact(observer_body, {&body});
    obs_contact.updateConfiguration();
    InteractionDynamics<CorrectInterpolationKernelWeights>{obs_contact}.exec();
    ObservedQuantityRecording<Vecd> disp_recorder("Displacement", obs_contact);
    ObservedQuantityRecording<Real> stress_recorder("VonMisesStress", obs_contact);
    disp_recorder.writeToFile(0);
    stress_recorder.writeToFile(0);

    // Simulation
    const Real end_time = 10.0;
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = algs_solid.time_step_size();
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

                dt = algs_solid.time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                algs_solid.stress_relaxation_first(dt);
                fix_shell_bc.exec();
                algs_solid.damping_exec(dt);
                fix_shell_bc.exec();
                algs_solid.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            disp_recorder.writeToFile(ite_output);
            stress_recorder.writeToFile(ite_output);
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

    // Output results
    std::cout << "Deflection: " << disp_recorder.getObservedQuantity()[0].y() << std::endl;
    std::cout << "Von Mises Stress: " << stress_recorder.getObservedQuantity()[0] << std::endl;
}