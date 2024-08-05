/**
 * @file 	Three ring impact.cpp
 * @brief 	This is the case file for the test of dynamic contacts between shell and solid.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"
using namespace SPH;

static constexpr Real unit_mm = 1e-3; // mm, g, ms

void ball_beam_contact(int dp_factor, bool use_pseudo_stiffness);

int main(int ac, char *av[])
{
    ball_beam_contact(1, true);
}

//------------------------------------------------------------------------------
void relax_solid(RealBody &body, InnerRelation &inner)
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
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
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
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> velocity_damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size;

    solid_algs(InnerRelation &inner_relation, double physical_viscosity, double damping_ratio = 0.2)
        : corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          velocity_damping(damping_ratio, inner_relation, "Velocity", physical_viscosity),
          computing_time_step_size(inner_relation.getSPHBody()){};

    void corrected_config() { corrected_configuration.exec(); }
    void stress_relaxation_first(Real dt) { stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second(Real dt) { stress_relaxation_second_half.exec(dt); }
    void damping(Real dt) { velocity_damping.exec(dt); }
    Real time_step_size() { return computing_time_step_size.exec(); }
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

class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name, const double length, const double height, const Vec2d &center) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(center), 0.5 * Vec2d(length, height), ShapeBooleanOps::add);
    }
};

class Sphere : public MultiPolygonShape
{
  public:
    explicit Sphere(const std::string &shape_name, const double radius, const Vec2d &center) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, radius, 100, ShapeBooleanOps::add);
    }
};

class VelocityBoundaryCondition : public BaseLocalDynamics<BodyPartByParticle>, public DataDelegateSimple
{
  private:
    StdLargeVec<Vec2d> *vel_;
    Vec2d velocity_;

  public:
    VelocityBoundaryCondition(BodyPartByParticle &body_part, Vec2d velocity)
        : BaseLocalDynamics<BodyPartByParticle>(body_part), DataDelegateSimple(body_part.getSPHBody()),
          vel_(this->particles_->template registerSharedVariable<Vec2d>("Velocity")),
          velocity_(std::move(velocity)){};
    inline void update(size_t index_i, [[maybe_unused]] Real dt = 0.0)
    {
        (*vel_)[index_i] = velocity_;
    }
};

class SphereBoundaryGeometry : public BodyPartByParticle
{
  private:
    Real radius_;
    Real dp_;
    Vec2d center_;

  public:
    SphereBoundaryGeometry(SPHBody &body, const std::string &body_part_name, Real radius, Real dp, const Vec2d &center)
        : BodyPartByParticle(body, body_part_name),
          radius_(radius), dp_(dp), center_(center)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&SphereBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (pos_[index_i].y() < center_.y() - radius_ + dp_)
            body_part_particles_.push_back(index_i);
    };
};

class BeamBoundaryGeometry : public BodyPartByParticle
{
  private:
    Real half_length_;
    Real dp_;
    Vec2d center_;

  public:
    BeamBoundaryGeometry(SPHBody &body, const std::string &body_part_name, Real length, Real dp, const Vec2d &center)
        : BodyPartByParticle(body, body_part_name),
          half_length_(0.5 * length), dp_(dp), center_(center)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BeamBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (pos_[index_i].x() < center_.x() - half_length_ + dp_ ||
            pos_[index_i].x() > center_.x() + half_length_ - dp_)
            body_part_particles_.push_back(index_i);
    };
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

Real get_E_from_K_and_G(Real K, Real G) { return 9.0 * K * G / (3.0 * K + G); }
Real get_nu_from_K_and_G(Real K, Real G) { return (3 * K - 2 * G) / (2.0 * (3 * K + G)); }
Real get_E_from_K_and_nu(Real K, Real nu) { return 3.0 * K * (1.0 - 2.0 * nu); }
//------------------------------------------------------------------------------
void ball_beam_contact(int dp_factor, bool use_pseudo_stiffness)
{
    const Real end_time = 1;
    const Vec2d velocity = 3 * Vec2d::UnitY();

    // geometry
    const Real beam_length = 40;
    const Real beam_height = 2;
    const Real sphere_radius = 3;

    // resolution
    const Real dp = beam_height / (4.0 * dp_factor);

    // material properties
    const Real rho_beam = 1e3 * std::pow(unit_mm, 2); // 1e3 kg/m^3
    const Real bulk_modulus_beam = 3.33;              // 3.33 MPa
    const Real shear_modulus_beam = 0.34;             // 0.34 MPa
    const Real youngs_modulus_beam = get_E_from_K_and_G(bulk_modulus_beam, shear_modulus_beam);
    const Real poisson_ratio_beam = get_nu_from_K_and_G(bulk_modulus_beam, shear_modulus_beam);

    const Real rho_sphere = 1e3 * std::pow(unit_mm, 2); // 1e3 kg/m^3
    const Real poisson_ratio_sphere = 0.4;              // 0.4
    const Real youngs_modulus_sphere = use_pseudo_stiffness
                                           ? get_E_from_K_and_nu(bulk_modulus_beam, poisson_ratio_sphere)
                                           : 900; // 900 MPa

    // Material models
    auto material_beam = makeShared<NeoHookeanSolid>(rho_beam, youngs_modulus_beam, poisson_ratio_beam);
    auto material_sphere = makeShared<SaintVenantKirchhoffSolid>(rho_sphere, youngs_modulus_sphere, poisson_ratio_sphere);

    std::cout << "Youngs modulus of the beam: " << youngs_modulus_beam << std::endl;
    std::cout << "Sound speed of the beam: " << material_beam->ReferenceSoundSpeed() << std::endl;
    std::cout << "Youngs modulus of the sphere: " << youngs_modulus_sphere << std::endl;
    std::cout << "Sound speed of the sphere: " << material_sphere->ReferenceSoundSpeed() << std::endl;
    std::cout << "Impact speed: " << velocity << std::endl;

    const Real physical_viscosity_beam = get_physical_viscosity_general(rho_beam, youngs_modulus_beam, beam_height);
    const Real physical_viscosity_sphere = get_physical_viscosity_general(rho_sphere, youngs_modulus_sphere, sphere_radius);

    // Import meshes
    const Vec2d beam_translation = Vec2d::Zero();
    auto mesh_beam = std::make_shared<Beam>("Beam", beam_length, beam_height, beam_translation);

    const Vec2d sphere_translation = -(sphere_radius + 0.5 * beam_height) * Vec2d::UnitY();
    auto mesh_sphere = std::make_shared<Sphere>("Sphere", sphere_radius, sphere_translation);

    // System bounding box
    BoundingBox bb_system = union_bounding_box(mesh_beam->getBounds(), mesh_sphere->getBounds());

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody beam_body(system, mesh_beam);
    beam_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    beam_body.assignMaterial(material_beam.get());
    beam_body.generateParticles<BaseParticles, Lattice>();

    SolidBody sphere_body(system, mesh_sphere);
    sphere_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    sphere_body.assignMaterial(material_sphere.get());
    sphere_body.generateParticles<BaseParticles, Lattice>();

    // Inner relation
    InnerRelation inner_beam(beam_body);
    InnerRelation inner_sphere(sphere_body);

    // relax solid
    relax_solid(beam_body, inner_beam);
    relax_solid(sphere_body, inner_sphere);

    // Methods
    solid_algs beam_algs(inner_beam, physical_viscosity_beam, 0.2);
    solid_algs sphere_algs(inner_sphere, physical_viscosity_sphere, 0.2);

    // Contact
    contact_algs contact_beam_sphere(beam_body, {&sphere_body});
    contact_algs contact_sphere_beam(sphere_body, {&beam_body});

    // Boundary conditions
    BeamBoundaryGeometry beam_fix_bc_part(beam_body, "BeamClampingBC", beam_length, dp, beam_translation);
    SimpleDynamics<FixBodyPartConstraint> beam_fix_bc(beam_fix_bc_part);

    SphereBoundaryGeometry sphere_vel_bc_part(sphere_body, "SphereClampingBC", sphere_radius, dp, sphere_translation);
    SimpleDynamics<VelocityBoundaryCondition> sphere_vel_bc(sphere_vel_bc_part, velocity);

    // Check
    auto check_nan = [&](BaseParticles &particles)
    {
        const auto &pos_ = particles.ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("RepulsionForce");
    sphere_body.getBaseParticles().addVariableToWrite<Vec2d>("RepulsionForce");
    beam_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    sphere_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    beam_algs.corrected_config();
    sphere_algs.corrected_config();

    contact_beam_sphere.config_update();
    contact_sphere_beam.config_update();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min(beam_algs.time_step_size(), sphere_algs.time_step_size());
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

                contact_beam_sphere.density_update();
                contact_sphere_beam.density_update();
                contact_beam_sphere.force_update();
                contact_sphere_beam.force_update();

                dt = std::min(beam_algs.time_step_size(), sphere_algs.time_step_size());
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                beam_algs.stress_relaxation_first(dt);
                sphere_algs.stress_relaxation_first(dt);

                sphere_vel_bc.exec();
                beam_fix_bc.exec();

                beam_algs.damping(dt);
                sphere_algs.damping(dt);

                sphere_vel_bc.exec();
                beam_fix_bc.exec();

                beam_algs.stress_relaxation_second(dt);
                sphere_algs.stress_relaxation_second(dt);

                beam_body.updateCellLinkedList();
                sphere_body.updateCellLinkedList();

                contact_beam_sphere.config_update();
                contact_sphere_beam.config_update();

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                { // checking if any position has become nan
                    check_nan(beam_body.getBaseParticles());
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