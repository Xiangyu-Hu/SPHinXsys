/**
 * @file 	test_3d_shell_to_solid_coupling.cpp
 * @brief 	This is the case file for the test of solid-shell coupling.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"

using namespace SPH;

void test_shell_to_solid_coupling();

int main(int ac, char *av[])
{
    test_shell_to_solid_coupling();
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

namespace SPH
{
class ShellPlate;
template <>
class ParticleGenerator<SurfaceParticles, ShellPlate> : public ParticleGenerator<SurfaceParticles>
{
    SharedPtr<Shape> mesh_;
    const Real thickness_;

  public:
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, SharedPtr<Shape> mesh, Real thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          mesh_(mesh),
          thickness_(thickness) {};
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
                addSurfaceProperties(Vec3d::UnitZ(), thickness_);
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

class ExternalForce : public BaseForcePrior<BodyPartByParticle>
{
  private:
    Vec3d force_;

  public:
    ExternalForce(BodyPartByParticle &body_part, const Vec3d &force)
        : BaseForcePrior<BodyPartByParticle>(body_part, "ExternalForce"),
          force_(force) {};
    void update(size_t index_i, Real dt = 0.0)
    {
        current_force_[index_i] = force_;
        BaseForcePrior<BodyPartByParticle>::update(index_i, dt);
    }
};

class CouplingConstraint : public BaseLocalDynamics<BodyPartByParticle>,
                           public DataDelegateContact
{
  protected:
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
    Vecd *vel_;
    StdVec<Vecd *> contact_vel_;

  public:
    SolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
        : MotionConstraint<BodyPartByParticle>(body_part),
          DataDelegateContact(contact_relation),
          coupling_ids_(particles_->getVariableDataByName<int>("CouplingID")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
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
                             public DataDelegateContact
{
  private:
    int *coupling_ids_;
    Vecd *n0_;
    StdVec<Matd *> contact_F_;
    StdVec<Matd *> contact_B_;
    StdVec<Real> contact_A0_;
    StdVec<ElasticSolid *> contact_materials_;

  public:
    ShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
        : BaseForcePrior<BodyPartByParticle>(body_part, "CouplingForce"),
          DataDelegateContact(contact_relation),
          coupling_ids_(particles_->getVariableDataByName<int>("CouplingID")),
          n0_(particles_->getVariableDataByName<Vecd>("InitialNormalDirection"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_A0_.emplace_back(pow(contact_bodies_[k]->getSPHBodyResolutionRef(), Dimensions - 1)); // initial area
            contact_F_.emplace_back(contact_particles_[k]->getVariableDataByName<Matd>("DeformationGradient"));
            contact_B_.emplace_back(contact_particles_[k]->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"));
            contact_materials_.emplace_back(&DynamicCast<ElasticSolid>(this, contact_bodies_[k]->getBaseMaterial()));
        }
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        int coupling_id = coupling_ids_[index_i];
        const Matd stress_PK1_B = contact_materials_[0]->StressPK1(contact_F_[0][coupling_id], coupling_id) * contact_B_[0][coupling_id];
        const Vecd force = stress_PK1_B * n0_[index_i] * contact_A0_[0];

        current_force_[index_i] = force;
        BaseForcePrior<BodyPartByParticle>::update(index_i, dt);
    }
};

void test_shell_to_solid_coupling()
{
    constexpr Real unit_mm = 1e-3;
    constexpr Real scale = 1.0;

    // Geometry
    const Real length = 10 * scale;
    const Real width = 5 * scale;
    const Real thickness_solid = 2 * scale;
    const Real dp = thickness_solid / 4.0;
    const Real thickness_shell = dp;

    // Material properties
    const Real rho = 1000 * pow(unit_mm, 2);
    const Real youngs_modulus_shell = 100; // 100 MPa
    const Real youngs_modulus_solid = 10;  // 1 MPa
    const Real poisson_ratio = 0.3;

    // Import meshes
    const Vec3d solid_halfsize = 0.5 * Vec3d(length, width, thickness_solid + dp);
    const Vec3d solid_translation = (solid_halfsize.z() - dp) * Vec3d::UnitZ();
    auto solid_mesh = makeShared<TransformShape<GeometricShapeBox>>(Transform(solid_translation), solid_halfsize, "solid");

    const Vec3d shell_halfsize = 0.5 * Vec3d(length, width, dp);
    const Vec3d shell_translation = -0.5 * dp * Vec3d::UnitZ();
    auto shell_mesh = makeShared<TransformShape<GeometricShapeBox>>(Transform(shell_translation), shell_halfsize, "shell");

    // System bounding box
    BoundingBox bb_system = union_bounding_box(solid_mesh->getBounds(), shell_mesh->getBounds());

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Body
    SolidBody solid_body(system, solid_mesh);
    solid_body.defineMaterial<NeoHookeanSolid>(rho, youngs_modulus_solid, poisson_ratio);
    solid_body.generateParticles<BaseParticles, Lattice>();

    SolidBody shell_body(system, shell_mesh);
    shell_body.defineMaterial<NeoHookeanSolid>(rho, youngs_modulus_shell, poisson_ratio);
    shell_body.generateParticles<SurfaceParticles, ShellPlate>(shell_mesh, thickness_shell);

    // algorithms
    solid_algs algs_solid(solid_body, get_physical_viscosity_general(rho, youngs_modulus_solid, thickness_solid));
    shell_algs algs_shell(shell_body, get_physical_viscosity_general(rho, youngs_modulus_shell, thickness_shell));

    // Boundary condition
    // fix the left end of the solid
    SolidBodyPart solid_fix_part(solid_body, "SolidFixedPart",
                                 [x0 = solid_mesh->getBounds().first_.x(), dp](Vec3d &pos)
                                 { return pos.x() < x0 + dp; });
    SimpleDynamics<FixBodyPartConstraint> fix_bc_solid(solid_fix_part);

    // fix the left end of shell
    SolidBodyPart shell_fix_part(shell_body, "ShellFixedPart",
                                 [x0 = shell_mesh->getBounds().first_.x(), dp](Vec3d &pos)
                                 { return pos.x() < x0 + dp; });
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> fix_bc_shell(shell_fix_part);

    // pressure on the upper surface of solid
    Vec3d pressure = -5e-2 * Vec3d::UnitZ(); // 2e-2 N/mm^2
    Vec3d force = pressure * dp * dp;
    SolidBodyPart solid_upper_surface(solid_body, "SolidUpperSurface",
                                      [z0 = solid_mesh->getBounds().second_.z(), dp](Vec3d &pos)
                                      { return pos.z() > z0 - dp; });
    SimpleDynamics<ExternalForce> pressure_bc(solid_upper_surface, force);

    // Coupling conditions
    SolidBodyPart solid_coupling_part(solid_body, "SolidCouplingPart",
                                      [z0 = solid_mesh->getBounds().first_.z(), dp](Vec3d &pos)
                                      { return pos.z() < z0 + dp; });
    std::cout << "Number of coupled solid ids: " << solid_coupling_part.body_part_particles_.size() << std::endl;
    ContactRelation solid_shell_contact(solid_body, {&shell_body});
    SimpleDynamics<CouplingConstraint> solid_coupling(solid_coupling_part, solid_shell_contact);
    SimpleDynamics<SolidVelocityConstraint> solid_coupling_vel_bc(solid_coupling_part, solid_shell_contact);

    SolidBodyPart shell_coupling_part(shell_body, "ShellCouplingPart",
                                      [z1 = shell_mesh->getBounds().second_.z(), dp](Vec3d &pos)
                                      { return pos.z() > z1 - dp; });
    std::cout << "Number of coupled shell ids: " << shell_coupling_part.body_part_particles_.size() << std::endl;
    ContactRelation shell_solid_contact(shell_body, {&solid_body});
    SimpleDynamics<CouplingConstraint> shell_coupling(shell_coupling_part, shell_solid_contact);
    SimpleDynamics<ShellForceConstraint> shell_coupling_force_bc(shell_coupling_part, shell_solid_contact);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs_solid.corrected_config();
    algs_shell.corrected_config();

    solid_coupling.exec();
    shell_coupling.exec();

    pressure_bc.exec();

    // Output
    solid_body.getBaseParticles().addVariableToWrite<int>("CouplingID");
    shell_body.getBaseParticles().addVariableToWrite<int>("CouplingID");
    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    shell_body.getBaseParticles().addVariableToWrite<Vec3d>("Velocity");
    solid_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    shell_body.getBaseParticles().addVariableToWrite<Mat3d>("DeformationGradient");
    solid_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    shell_body.getBaseParticles().addVariableToWrite<Vec3d>("InitialNormalDirection");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell_body);
    vtp_output.writeToFile(0);

    // Simulation
    const Real end_time = 2.0;
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min(algs_shell.time_step_size(), algs_solid.time_step_size());
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
                shell_coupling_force_bc.exec();

                dt = std::min(algs_shell.time_step_size(), algs_solid.time_step_size());
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // update shell first
                algs_shell.stress_relaxation_first(dt);
                fix_bc_shell.exec();
                algs_shell.damping_exec(dt);
                fix_bc_shell.exec();
                algs_shell.stress_relaxation_second(dt);

                // update solid
                algs_solid.stress_relaxation_first(dt);
                fix_bc_solid.exec();
                solid_coupling_vel_bc.exec();
                algs_solid.damping_exec(dt);
                fix_bc_solid.exec();
                solid_coupling_vel_bc.exec();
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