/**
 * @file 	test_3d_shell_to_solid_coupling.cpp
 * @brief 	This is the case file for the test of solid-shell coupling.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"

using namespace SPH;

void run_solid(int res_factor, Real stiffness_ratio, int load_type = 0);
void run_shell_to_solid_coupling(int res_factor, Real stiffness_ratio, int load_type = 0);

int main(int ac, char *av[])
{
    run_shell_to_solid_coupling(1, 0.1, 1);
    // run_solid(1, 0.1, 1);
}

struct solid_algs
{
    InnerRelation inner_relation;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    SimpleDynamics<NormalDirectionFromBodyShape> initial_normal_direction;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> normal_direction;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size;

    explicit solid_algs(SolidBody &body, Real physical_viscosity)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          initial_normal_direction(body),
          normal_direction(body),
          damping(0.2, inner_relation, "Velocity", physical_viscosity),
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
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> velocity_damping;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> rotation_damping;
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size;

    explicit shell_algs(SolidBody &body, Real physical_viscosity)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          normal_direction(inner_relation.getSPHBody()),
          velocity_damping(0.2, inner_relation, "Velocity", physical_viscosity),
          rotation_damping(0.2, inner_relation, "AngularVelocity", physical_viscosity),
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

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

namespace SPH
{
class ShellPlate;
template <>
class ParticleGenerator<SurfaceParticles, ShellPlate> : public ParticleGenerator<SurfaceParticles>
{
  private:
    SharedPtr<Shape> mesh_;
    const Real thickness_;
    Vec3d normal_;

  public:
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, SharedPtr<Shape> mesh, Real thickness, bool normal_flip = false)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          mesh_(mesh),
          thickness_(thickness),
          normal_((normal_flip ? -1 : 1) * Vec3d::UnitZ()) {};
    void prepareGeometricData() override
    {
        Real dp = (mesh_->getBounds().second_.z() - mesh_->getBounds().first_.z());
        Real z0 = 0.5 * (mesh_->getBounds().first_.z() + mesh_->getBounds().second_.z());
        Real x0 = mesh_->getBounds().first_.x();
        Real y0 = mesh_->getBounds().first_.y();
        Real x1 = mesh_->getBounds().second_.x();
        Real y1 = mesh_->getBounds().second_.y();
        Real x = x0 + 0.5 * dp;
        while (x < x1)
        {
            Real y = y0 + 0.5 * dp;
            while (y < y1)
            {
                addPositionAndVolumetricMeasure(Vec3d(x, y, z0), dp);
                addSurfaceProperties(normal_, thickness_);
                y += dp;
            }
            x += dp;
        }
    }
};
} // namespace SPH

class SolidBodyPart : public BodyPartByParticle
{
  private:
    std::function<bool(Vec3d &)> contain_;

  public:
    SolidBodyPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec3d &)> contain)
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

class VelocityBoundaryCondition : public MotionConstraint<BodyPartByParticle>
{
  private:
    Real &physical_time;
    std::function<Vec3d(Real time)> vel_func_;

  public:
    VelocityBoundaryCondition(BodyPartByParticle &body_part, std::function<Vec3d(Real time)> vel_func)
        : MotionConstraint<BodyPartByParticle>(body_part),
          physical_time(*body_part.getSPHSystem().getSystemVariableDataByName<Real>("PhysicalTime")),
          vel_func_(std::move(vel_func)) {};
    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = vel_func_(physical_time);
    }
};

class CouplingConstraint : public BaseLocalDynamics<BodyPartByParticle>,
                           public DataDelegateContact
{
  private:
    int *coupling_ids_;

  public:
    explicit CouplingConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          DataDelegateContact(contact_relation),
          coupling_ids_(particles_->registerStateVariable<int>("CouplingID"))
    {
        if (contact_configuration_.size() != 1)
            throw std::runtime_error("CouplingConstraint currently only supports one contact configuration.");
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        const auto &contact_neighbors = (*contact_configuration_[0])[index_i];
        int coupling_id = [&contact_neighbors]()
        {
            int id = -1;
            Real r_min = std::numeric_limits<Real>::max();
            for (size_t n = 0; n < contact_neighbors.current_size_; n++)
            {
                if (contact_neighbors.r_ij_[n] < r_min)
                {
                    r_min = contact_neighbors.r_ij_[n];
                    id = contact_neighbors.j_[n];
                }
            }
            return id;
        }();
        coupling_ids_[index_i] = coupling_id;
    }
};

class SolidVelocityConstraint : public MotionConstraint<BodyPartByParticle>,
                                public DataDelegateContact
{
  private:
    int *coupling_ids_;
    StdVec<Vecd *> contact_vel_;

  public:
    SolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
        : MotionConstraint<BodyPartByParticle>(body_part),
          DataDelegateContact(contact_relation),
          coupling_ids_(particles_->getVariableDataByName<int>("CouplingID"))
    {
        for (auto *contact_particle : contact_particles_)
            contact_vel_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Velocity"));
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        int coupling_id = coupling_ids_[index_i];
        vel_[index_i] = contact_vel_[0][coupling_id];
    }
};

class ShellForceConstraint : public BaseForcePrior<BodyPartByParticle>,
                             public DataDelegateInner,
                             public DataDelegateContact
{
  private:
    int *coupling_ids_;
    Real *mass_;
    Real *Vol_;
    Matd *B_;
    Matd *transformation_matrix_;
    StdVec<Matd *> contact_F_;
    StdVec<ElasticSolid *> contact_materials_;
    Real inv_rho0_;

  public:
    ShellForceConstraint(BodyPartByParticle &body_part, BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : BaseForcePrior<BodyPartByParticle>(body_part, "CouplingForce"),
          DataDelegateInner(inner_relation),
          DataDelegateContact(contact_relation),
          coupling_ids_(particles_->getVariableDataByName<int>("CouplingID")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          transformation_matrix_(particles_->getVariableDataByName<Matd>("TransformationMatrix")),
          inv_rho0_(1.0 / sph_body_.getBaseMaterial().ReferenceDensity())
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_F_.emplace_back(contact_particles_[k]->getVariableDataByName<Matd>("DeformationGradient"));
            // contact_B_.emplace_back(contact_particles_[k]->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"));
            contact_materials_.emplace_back(&DynamicCast<ElasticSolid>(this, contact_bodies_[k]->getBaseMaterial()));
        }
    };

    void interaction(size_t index_i, Real dt = 0.0)
    {
        int coupling_id_i = coupling_ids_[index_i];
        Vecd force = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Matd stress_i = contact_materials_[0]->StressPK1(contact_F_[0][coupling_id_i], coupling_id_i);
        // Matd global_B_i = transformation_matrix_[index_i].transpose() * B_[index_i] * transformation_matrix_[index_i];
        // Matd stress_B_i = stress_i * global_B_i.transpose();
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            int coupling_id_j = coupling_ids_[index_j];
            Matd stress_j = contact_materials_[0]->StressPK1(contact_F_[0][coupling_id_j], coupling_id_j);
            // Matd global_B_j = transformation_matrix_[index_j].transpose() * B_[index_j] * transformation_matrix_[index_j];
            // Matd stress_B_j = stress_j * global_B_j.transpose();
            Vecd e_ij = inner_neighborhood.e_ij_[n];
            force += mass_[index_i] * inv_rho0_ * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                     (stress_i + stress_j) * e_ij;
        }
        current_force_[index_i] = force;
    }
};

struct coupling_algs
{
    SolidBodyPart part;
    ContactRelation contact_relation;
    SimpleDynamics<CouplingConstraint> coupling;

    coupling_algs(SolidBody &body, SolidBody &contact_body, std::string part_name, std::function<bool(Vec3d &)> contain)
        : part(body, std::move(part_name), std::move(contain)),
          contact_relation(body, {&contact_body}),
          coupling(part, contact_relation) {}

    void initialize() { coupling.exec(); }
};

struct solid_coupling_algs
{
    coupling_algs base_algs;
    SimpleDynamics<SolidVelocityConstraint> vel_bc;

    solid_coupling_algs(SolidBody &body, SolidBody &contact_body, std::string part_name, std::function<bool(Vec3d &)> contain)
        : base_algs(body, contact_body, std::move(part_name), std::move(contain)),
          vel_bc(base_algs.part, base_algs.contact_relation) {}

    void initialize() { base_algs.coupling.exec(); }
    void exec() { vel_bc.exec(); }
};

struct shell_coupling_algs
{
    coupling_algs base_algs;
    InteractionWithUpdate<ShellForceConstraint> force_bc;

    shell_coupling_algs(BaseInnerRelation &inner, SolidBody &contact_body, std::string part_name, std::function<bool(Vec3d &)> contain)
        : base_algs(DynamicCast<SolidBody>(this, inner.getSPHBody()), contact_body, std::move(part_name), std::move(contain)),
          force_bc(base_algs.part, inner, base_algs.contact_relation) {}

    void exec() { force_bc.exec(); }
};

namespace plate
{
constexpr Real unit_mm = 1e-3;
constexpr Real scale = 1.0;

struct plate_parameters
{
    // Geometry
    const Real length = 50 * scale;
    const Real width = length;
    const Real height = 20 * scale;
    const Real thickness_shell = 4 * scale;

    // Material properties
    const Real rho = 1000 * pow(unit_mm, 2);
    const Real youngs_modulus_shell = 1e9 * pow(unit_mm, 2); // 10 GPa
    const Real poisson_ratio = 0.3;

    inline auto get_solid_extended_mesh(Real dp)
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, height + 2 * dp);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    }

    inline auto get_whole_mesh()
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, height + 2 * thickness_shell);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    }

    inline auto get_upper_shell_mesh(Real dp)
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, dp);
        const Vec3d translation = 0.5 * (height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "upper_shell");
    }

    inline auto get_lower_shell_mesh(Real dp)
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, dp);
        const Vec3d translation = -0.5 * (height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "lower_shell");
    }
};

struct boundary_condition
{
    Real max_elongation = 10 * scale;
    Real speed = 10.0; // mm/ms
    Real elongation_time = 0.5 * max_elongation / speed;

    Vec3d direction_left;

    SolidBodyPart part_left;
    SolidBodyPart part_right;
    SimpleDynamics<VelocityBoundaryCondition> left_bc;
    SimpleDynamics<VelocityBoundaryCondition> right_bc;

    Vec3d get_velocity(Real time, const Vec3d &direction) const
    {
        Real vel = time < elongation_time ? speed : 0;
        return vel * direction;
    }

    boundary_condition(SolidBody &body, int problem_type)
        : direction_left(problem_type == 0 ? -Vec3d::UnitX() : -Vec3d::UnitZ()),
          part_left(body, "LeftPart", [x0 = body.getInitialShape().getBounds().first_.x(), dp = body.getSPHBodyResolutionRef()](Vec3d &pos)
                    { return pos.x() < x0 + dp; }),
          part_right(body, "RightPart", [x1 = body.getInitialShape().getBounds().second_.x(), dp = body.getSPHBodyResolutionRef()](Vec3d &pos)
                     { return pos.x() > x1 - dp; }),
          left_bc(part_left, [this](Real time)
                  { return get_velocity(time, direction_left); }),
          right_bc(part_right,
                   [this](Real time)
                   { return get_velocity(time, -direction_left); }) {};

    void exec()
    {
        left_bc.exec();
        right_bc.exec();
    }
};
}; // namespace plate

void run_shell_to_solid_coupling(int res_factor, Real stiffness_ratio, int load_type)
{
    // parameters
    plate::plate_parameters params;
    Real dp = params.thickness_shell / (4.0 * res_factor);
    Real youngs_modulus_solid = params.youngs_modulus_shell * stiffness_ratio;

    // Import meshes
    auto solid_mesh = params.get_solid_extended_mesh(dp);
    auto lower_shell_mesh = params.get_lower_shell_mesh(dp);
    auto upper_shell_mesh = params.get_upper_shell_mesh(dp);

    // System bounding box
    BoundingBox bb_system = solid_mesh->getBounds();

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Body
    SolidBody solid_body(system, solid_mesh);
    solid_body.defineMaterial<LinearElasticSolid>(params.rho, youngs_modulus_solid, params.poisson_ratio);
    solid_body.generateParticles<BaseParticles, Lattice>();

    SolidBody upper_shell_body(system, upper_shell_mesh);
    upper_shell_body.defineMaterial<LinearElasticSolid>(params.rho, params.youngs_modulus_shell, params.poisson_ratio);
    upper_shell_body.generateParticles<SurfaceParticles, ShellPlate>(upper_shell_mesh, params.thickness_shell, true);

    SolidBody lower_shell_body(system, lower_shell_mesh);
    lower_shell_body.defineMaterial<LinearElasticSolid>(params.rho, params.youngs_modulus_shell, params.poisson_ratio);
    lower_shell_body.generateParticles<SurfaceParticles, ShellPlate>(lower_shell_mesh, params.thickness_shell, false);

    auto &solid_material = *dynamic_cast<ElasticSolid *>(&solid_body.getBaseMaterial());
    std::cout << "solid sound speed: " << solid_material.ReferenceSoundSpeed() << std::endl;

    // algorithms
    Real eta = get_physical_viscosity_general(params.rho, params.youngs_modulus_shell, params.height + 2 * params.thickness_shell);
    solid_algs algs_solid(solid_body, eta);
    shell_algs algs_shell_upper(upper_shell_body, eta);
    shell_algs algs_shell_lower(lower_shell_body, eta);

    // Boundary condition
    plate::boundary_condition solid_bc(solid_body, load_type);
    plate::boundary_condition lower_shell_bc(lower_shell_body, load_type);
    plate::boundary_condition upper_shell_bc(upper_shell_body, load_type);

    // Coupling conditions
    solid_coupling_algs solid_coupling_algs_lower(solid_body, lower_shell_body, "SolidLowerCouplingPart",
                                                  [z0 = solid_mesh->getBounds().first_.z(), dp](Vec3d &pos)
                                                  { return pos.z() < z0 + dp; });
    solid_coupling_algs solid_coupling_algs_upper(solid_body, upper_shell_body, "SolidUpperCouplingPart",
                                                  [z1 = solid_mesh->getBounds().second_.z(), dp](Vec3d &pos)
                                                  { return pos.z() > z1 - dp; });
    shell_coupling_algs shell_coupling_algs_lower(algs_shell_lower.inner_relation, solid_body, "LowerShellCouplingPart", [z1 = lower_shell_mesh->getBounds().second_.z(), dp](Vec3d &pos)
                                                  { return pos.z() > z1 - dp; });
    shell_coupling_algs shell_coupling_algs_upper(algs_shell_upper.inner_relation, solid_body, "UpperShellCouplingPart",
                                                  [z0 = upper_shell_mesh->getBounds().first_.z(), dp](Vec3d &pos)
                                                  { return pos.z() < z0 + dp; });

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs_solid.corrected_config();
    algs_shell_lower.corrected_config();
    algs_shell_upper.corrected_config();

    solid_coupling_algs_upper.initialize();
    solid_coupling_algs_lower.initialize();
    shell_coupling_algs_upper.base_algs.coupling.exec();
    shell_coupling_algs_lower.base_algs.coupling.exec();

    // Output
    solid_body.getBaseParticles().addVariableToWrite<int>("CouplingID");
    upper_shell_body.getBaseParticles().addVariableToWrite<int>("CouplingID");
    lower_shell_body.getBaseParticles().addVariableToWrite<int>("CouplingID");
    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    upper_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    lower_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    solid_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    upper_shell_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    lower_shell_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    upper_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    lower_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    upper_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("CouplingForce");
    lower_shell_body.getBaseParticles().addVariableToWrite<Vec3d>("CouplingForce");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(upper_shell_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(lower_shell_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(upper_shell_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(lower_shell_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(upper_shell_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(lower_shell_body);
    vtp_output.writeToFile(0);

    // Simulation
    const Real end_time = 4.0;
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min({algs_shell_upper.time_step_size(), algs_shell_lower.time_step_size(), algs_solid.time_step_size()});
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

                dt = std::min({algs_shell_upper.time_step_size(), algs_shell_lower.time_step_size(), algs_solid.time_step_size()});
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // update solid first half
                algs_solid.stress_relaxation_first(dt);

                // compute coupling force on shell
                shell_coupling_algs_lower.exec();
                shell_coupling_algs_upper.exec();

                // update shell
                algs_shell_upper.stress_relaxation_first(dt);
                upper_shell_bc.exec();
                algs_shell_upper.damping_exec(dt);
                upper_shell_bc.exec();
                algs_shell_upper.stress_relaxation_second(dt);

                algs_shell_lower.stress_relaxation_first(dt);
                lower_shell_bc.exec();
                algs_shell_lower.damping_exec(dt);
                lower_shell_bc.exec();
                algs_shell_lower.stress_relaxation_second(dt);

                // update solid second half with shell kinematic constraints
                solid_coupling_algs_lower.exec();
                solid_coupling_algs_upper.exec();
                solid_bc.exec();
                algs_solid.damping_exec(dt);
                solid_coupling_algs_lower.exec();
                solid_coupling_algs_upper.exec();
                solid_bc.exec();
                algs_solid.stress_relaxation_second(dt);

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

void run_solid(int res_factor, Real stiffness_ratio, int load_type)
{
    // parameters
    plate::plate_parameters params;
    Real dp = params.thickness_shell / (4.0 * res_factor);
    Real youngs_modulus_solid = params.youngs_modulus_shell * stiffness_ratio;

    // Import meshes
    auto solid_part_mesh = params.get_solid_extended_mesh(0);
    auto whole_mesh = params.get_whole_mesh();

    // System bounding box
    BoundingBox bb_system = whole_mesh->getBounds();

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Body
    SolidBody solid_body(system, whole_mesh);
    solid_body.defineMaterial<CompositeSolid>(params.rho);
    solid_body.generateParticles<BaseParticles, Lattice>();
    { // assign material ids
        auto *composite_solid = dynamic_cast<CompositeSolid *>(&solid_body.getBaseMaterial());
        composite_solid->add<LinearElasticSolid>(params.rho, youngs_modulus_solid, params.poisson_ratio);
        composite_solid->add<LinearElasticSolid>(params.rho, params.youngs_modulus_shell, params.poisson_ratio);
        int *material_id = solid_body.getBaseParticles().getVariableDataByName<int>("MaterialID");
        Vec3d *pos = solid_body.getBaseParticles().getVariableDataByName<Vec3d>("Position");
        for (size_t index_i = 0; index_i < solid_body.getBaseParticles().TotalRealParticles(); index_i++)
        {
            if (solid_part_mesh->checkContain(pos[index_i]))
                material_id[index_i] = 0;
            else
                material_id[index_i] = 1;
        }
    }

    // algorithms
    Real eta = get_physical_viscosity_general(params.rho, params.youngs_modulus_shell, params.height + 2 * params.thickness_shell);
    solid_algs algs_solid(solid_body, eta);

    // Boundary condition
    plate::boundary_condition solid_bc(solid_body, load_type);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs_solid.corrected_config();

    // Output

    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    solid_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    solid_body.getBaseParticles().addVariableToWrite<int>("MaterialID");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(solid_body);
    vtp_output.writeToFile(0);

    // Simulation
    const Real end_time = 4.0;
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

                // compute coupling force on shell
                dt = algs_solid.time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // update solid
                algs_solid.stress_relaxation_first(dt);
                solid_bc.exec();
                algs_solid.damping_exec(dt);
                solid_bc.exec();
                algs_solid.stress_relaxation_second(dt);

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