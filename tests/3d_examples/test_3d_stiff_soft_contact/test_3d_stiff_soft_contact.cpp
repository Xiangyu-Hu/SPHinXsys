#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <numeric>

using namespace SPH;

static constexpr Real unit_mm = 1e-3; // mm, g, ms

struct material_parameter
{
    Real density = 1 * std::pow(unit_mm, 2);          // 1kg/m^3
    Real Youngs_modulus = 5e5 * std::pow(unit_mm, 2); // 5e5 Pa
    Real poisson_ratio = 0.49;
};

// velocity ratio: velocity / soft_c
void test_stiff_soft(material_parameter &params_stiff, material_parameter &params_soft, double velocity_ratio);

inline double get_K_from_E_nu(double E, double nu) { return E / 3.0 / (1.0 - 2.0 * nu); }
inline double get_E_from_K_nu(double K, double nu) { return K * 3 * (1 - 2 * nu); }

TEST(stiff_soft_contact, ratio_x2)
{
    material_parameter params_stiff;
    material_parameter params_soft;
    params_soft.Youngs_modulus = 5e1 * std::pow(unit_mm, 2);
    test_stiff_soft(params_stiff, params_soft, 2);
}

TEST(stiff_soft_contact, ratio_x1)
{
    material_parameter params_stiff;
    material_parameter params_soft;
    params_soft.Youngs_modulus = 5e1 * std::pow(unit_mm, 2);
    test_stiff_soft(params_stiff, params_soft, 1);
}

TEST(stiff_soft_contact, ratio_x0_5)
{
    material_parameter params_stiff;
    material_parameter params_soft;
    params_soft.Youngs_modulus = 5e1 * std::pow(unit_mm, 2);
    test_stiff_soft(params_stiff, params_soft, 0.5);
}

TEST(stiff_soft_contact, same_K_ratio_x1)
{
    material_parameter params_stiff;
    material_parameter params_soft;
    params_soft.Youngs_modulus = 5e1 * std::pow(unit_mm, 2);
    params_stiff.poisson_ratio = 0.3;
    const double bulk_soft = get_K_from_E_nu(params_soft.Youngs_modulus, params_soft.poisson_ratio);
    const double E_stiff = get_K_from_E_nu(bulk_soft, params_stiff.poisson_ratio);
    params_stiff.Youngs_modulus = E_stiff;
    std::cout << "Youngs modulus of stiff material: " << E_stiff / std::pow(unit_mm, 2) << "Pa" << std::endl;
    std::cout << "Stiff / soft: " << E_stiff / params_soft.Youngs_modulus << std::endl;
    test_stiff_soft(params_stiff, params_soft, 1);
}

TEST(stiff_soft_contact, same_K_ratio_x0_1)
{
    material_parameter params_stiff;
    material_parameter params_soft;
    params_soft.Youngs_modulus = 5e1 * std::pow(unit_mm, 2);
    params_stiff.poisson_ratio = 0.3;
    const double bulk_soft = get_K_from_E_nu(params_soft.Youngs_modulus, params_soft.poisson_ratio);
    const double E_stiff = get_K_from_E_nu(bulk_soft, params_stiff.poisson_ratio);
    params_stiff.Youngs_modulus = E_stiff;
    std::cout << "Youngs modulus of stiff material: " << E_stiff / std::pow(unit_mm, 2) << "Pa" << std::endl;
    std::cout << "Stiff / soft: " << E_stiff / params_soft.Youngs_modulus << std::endl;
    test_stiff_soft(params_stiff, params_soft, 0.1);
}

int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    testing::GTEST_FLAG(filter) = "stiff_soft_contact.same_K_ratio_x0_1";
    return RUN_ALL_TESTS();
}

namespace SPH
{
class ShellHalfSphere;
template <>
class ParticleGenerator<SurfaceParticles, ShellHalfSphere> : public ParticleGenerator<SurfaceParticles>
{
    const Vec3d center_;
    const Real radius_;
    const Real dp_;
    const Real thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               const Vec3d &center, Real radius, Real dp, Real thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          center_(center),
          radius_(radius),
          dp_(dp),
          thickness_(thickness){};
    void prepareGeometricData() override
    {
        int number_of_particles_phi = int(Pi * radius_ / dp_) + 1;
        const Real dphi = Pi / Real(number_of_particles_phi - 1);
        for (int i = 0; i < number_of_particles_phi; i++)
        {
            const Real phi = i * dphi;
            const Real r = radius_ * sin(phi);
            int number_of_particles_theta = int(Pi * r / dp_) + 1;
            if (number_of_particles_theta == 1)
            {
                Vec3d pos(0, radius_ * cos(phi), 0);
                addPositionAndVolumetricMeasure(pos + center_, dp_ * dp_);
                addSurfaceProperties((pos - center_).normalized(), thickness_);
            }
            else
            {
                const Real dtheta = Pi / Real(number_of_particles_theta - 1);
                for (int j = 0; j < number_of_particles_theta; j++)
                {
                    const Real theta = j * dtheta;
                    Vec3d pos(r * sin(theta), radius_ * cos(phi), r * cos(theta));
                    addPositionAndVolumetricMeasure(pos + center_, dp_ * dp_);
                    addSurfaceProperties(pos.normalized(), thickness_);
                }
            }
        }
    }
};
} // namespace SPH

class BoundaryGeometry : public BodyPartByParticle
{
  private:
    Vec3d center_;
    Real length_;
    Real bc_length_;

  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name, const Vec3d &center, Real length, Real bc_length)
        : BodyPartByParticle(body, body_part_name),
          center_(center), length_(length), bc_length_(bc_length)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (pos_[index_i].x() > (center_.x() + 0.5 * length_ - bc_length_))
            body_part_particles_.push_back(index_i);
    };
};

class VelocityBoundaryCondition : public BaseLocalDynamics<BodyPartByParticle>, public DataDelegateSimple
{
  private:
    StdLargeVec<Vec3d> *vel_;
    Vec3d velocity_;

  public:
    VelocityBoundaryCondition(BodyPartByParticle &body_part, Vec3d velocity)
        : BaseLocalDynamics<BodyPartByParticle>(body_part), DataDelegateSimple(body_part.getSPHBody()),
          vel_(this->particles_->template registerSharedVariable<Vec3d>("Velocity")),
          velocity_(std::move(velocity)){};
    inline void update(size_t index_i, [[maybe_unused]] Real dt = 0.0)
    {
        (*vel_)[index_i] = velocity_;
    }
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

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

struct solid_algs
{
    InnerRelation inner_relation;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> velocity_damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size;

    solid_algs(RealBody &body, double physical_viscosity, double damping_ratio = 0.2)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          velocity_damping(damping_ratio, inner_relation, "Velocity", physical_viscosity),
          computing_time_step_size(body){};

    void config_update() { inner_relation.updateConfiguration(); }
    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void damping(Real dt) { velocity_damping.exec(dt); }
    Real time_step_size() { return computing_time_step_size.exec(); }
};

struct shell_algs
{
    InnerRelation inner_relation;
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> velocity_damping;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> rotation_damping;
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_s;
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size_s;

    shell_algs(RealBody &body, double physical_viscosity, double damping_ratio = 0.2)
        : inner_relation(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation, 3, true),
          stress_relaxation_second_half(inner_relation),
          velocity_damping(damping_ratio, inner_relation, "Velocity", physical_viscosity),
          rotation_damping(damping_ratio, inner_relation, "AngularVelocity", physical_viscosity),
          update_normal_s(body),
          computing_time_step_size_s(body){};

    void config_update() { inner_relation.updateConfiguration(); }
    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void damping(Real dt)
    {
        velocity_damping.exec(dt);
        rotation_damping.exec(dt);
    }
    void update_normal() { update_normal_s.exec(); }
    Real time_step_size() { return computing_time_step_size_s.exec(); }
};

struct shell_curvature_algs
{
    ShellInnerRelationWithContactKernel curvature_inner;
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> average_curvature;

    shell_curvature_algs(RealBody &sph_body, RealBody &contact_body)
        : curvature_inner(sph_body, contact_body),
          average_curvature(curvature_inner){};

    void config_update() { curvature_inner.updateConfiguration(); }
    void curvature_update() { average_curvature.exec(); }
};

struct contact_algs
{
    SurfaceContactRelation contact_relation;
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_density;
    std::unique_ptr<InteractionWithUpdate<solid_dynamics::ContactForce>> contact_force;

    contact_algs(RealBody &body, RealBodyVector contact_bodies, StdVec<bool> normal_corrections = {})
        : contact_relation(body, std::move(contact_bodies), std::move(normal_corrections)),
          contact_density(contact_relation)
    {
        for (auto *contact_body : contact_relation.contact_bodies_)
        {
            contact_body->getBaseParticles().registerSharedVariable<Real>("RepulsionFactor");
        }
        contact_force = std::make_unique<InteractionWithUpdate<solid_dynamics::ContactForce>>(contact_relation);
    }

    void config_update() { contact_relation.updateConfiguration(); }
    void density_update() { contact_density.exec(); }
    void force_update() { contact_force->exec(); }
};

void test_stiff_soft(material_parameter &params_stiff, material_parameter &params_soft, double velocity_ratio)
{
    const Real displacement_max = 5;

    // geometry
    const Real brick_length = 10;
    const Real brick_height = 5;
    const Real sphere_radius = 5;

    // resolutions
    const Real dp = brick_height / 10;
    const Real sphere_thickness = dp;

    // Material models
    auto material_brick = makeShared<SaintVenantKirchhoffSolid>(params_stiff.density, params_stiff.Youngs_modulus, params_stiff.poisson_ratio);
    auto material_sphere = makeShared<SaintVenantKirchhoffSolid>(params_soft.density, params_soft.Youngs_modulus, params_soft.poisson_ratio);

    // velocity
    const Real velocity = material_sphere->ReferenceSoundSpeed() * velocity_ratio;
    const Real end_time = displacement_max / velocity;

    std::cout << "Sound speed of the stiff material: " << material_brick->ReferenceSoundSpeed() << std::endl;
    std::cout << "Sound speed of the soft material: " << material_sphere->ReferenceSoundSpeed() << std::endl;
    std::cout << "Impact speed: " << velocity << std::endl;

    const Real physical_viscosity_brick = get_physical_viscosity_general(params_stiff.density, params_stiff.Youngs_modulus, dp);
    const Real physical_viscosity_sphere = get_physical_viscosity_general(params_soft.density, params_soft.Youngs_modulus, sphere_thickness);

    // Import meshes
    const auto brick_translation = Vec3d(0.5 * brick_length + 1.15 * dp, 0.0, 0.0);
    const auto sphere_translation = Vec3d(-sphere_radius - 1.15 * dp, 0.0, 0.0);
    auto mesh_brick = std::make_shared<TriangleMeshShapeBrick>(0.5 * Vec3d(brick_length, brick_height, brick_height), 5, brick_translation, "Brick");
    auto mesh_sphere = std::make_shared<TriangleMeshShapeSphere>(sphere_radius, 5, sphere_translation, "Sphere");

    // System bounding box
    BoundingBox bb_system = union_bounding_box(mesh_brick->getBounds(), mesh_sphere->getBounds());

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody brick_body(system, mesh_brick);
    brick_body.assignMaterial(material_brick.get());
    brick_body.generateParticles<BaseParticles, Lattice>();

    SolidBody sphere_body(system, mesh_sphere);
    sphere_body.assignMaterial(material_sphere.get());
    sphere_body.generateParticles<SurfaceParticles, ShellHalfSphere>(sphere_translation, sphere_radius, dp, sphere_thickness);

    // Methods
    solid_algs brick_algs(brick_body, physical_viscosity_brick, 0.2);
    shell_algs sphere_algs(sphere_body, physical_viscosity_sphere, 0.2);

    // Shell average curvature
    shell_curvature_algs sphere_curvature_algs(sphere_body, brick_body);

    // Contact
    contact_algs contact_brick_sphere(brick_body, {&sphere_body}, {true});
    contact_algs contact_sphere_brick(sphere_body, {&brick_body}, {true});

    // Boundary conditions
    BoundaryGeometry brick_vel_bc_part(brick_body, "BrickVelocityBC", brick_translation, brick_length, 0.25 * brick_length);
    SimpleDynamics<VelocityBoundaryCondition> brick_vel_bc(brick_vel_bc_part, Vec3d(-velocity, 0.0, 0.0));

    // Check
    auto check_nan = [&](BaseParticles &particles)
    {
        const auto &pos_ = particles.ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    sphere_body.getBaseParticles().addVariableToWrite<Vec3d>("NormalDirection");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    brick_algs.corrected_config();
    sphere_algs.corrected_config();

    sphere_curvature_algs.curvature_update();
    contact_brick_sphere.config_update();
    contact_sphere_brick.config_update();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min(brick_algs.time_step_size(), sphere_algs.time_step_size());
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

                contact_brick_sphere.density_update();
                contact_sphere_brick.density_update();
                contact_brick_sphere.force_update();
                contact_sphere_brick.force_update();

                dt = std::min(brick_algs.time_step_size(), sphere_algs.time_step_size());
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                brick_algs.stress_relaxation_first(dt);
                sphere_algs.stress_relaxation_first(dt);

                brick_vel_bc.exec();

                brick_algs.damping(dt);
                sphere_algs.damping(dt);

                brick_vel_bc.exec();

                brick_algs.stress_relaxation_second(dt);
                sphere_algs.stress_relaxation_second(dt);

                brick_body.updateCellLinkedList();
                sphere_body.updateCellLinkedList();

                sphere_algs.update_normal();

                sphere_curvature_algs.config_update();
                sphere_curvature_algs.curvature_update();

                contact_brick_sphere.config_update();
                contact_sphere_brick.config_update();

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                { // checking if any position has become nan
                    check_nan(brick_body.getBaseParticles());
                    check_nan(sphere_body.getBaseParticles());
                }
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