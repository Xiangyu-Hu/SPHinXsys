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

void soft_stiff_contact(int dp_factor, bool use_pseudo_stiffness, double nu_min, int stiff_E_factor = 1);

int main(int ac, char *av[])
{
    soft_stiff_contact(1, true, 0.1, 1);
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

class Arc : public MultiPolygonShape
{
  public:
    explicit Arc(const std::string &shape_name, const Vec2d &center, double inner_radius, double thickness) : MultiPolygonShape(shape_name)
    {
        const Real outer_radius = inner_radius + thickness;
        multi_polygon_.addACircle(center, outer_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(center, inner_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addABox(Transform(center - outer_radius * Vec2d::UnitX()), outer_radius * Vec2d::Ones(), ShapeBooleanOps::sub);
        multi_polygon_.addABox(Transform(center - outer_radius * Vec2d::UnitY()), outer_radius * Vec2d::Ones(), ShapeBooleanOps::sub);
    }
};

class Tube : public MultiPolygonShape
{
  public:
    explicit Tube(const std::string &shape_name, const Vec2d &center, double width, double length) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(center), 0.5 * Vec2d(width, length), ShapeBooleanOps::add);
        multi_polygon_.addACircle(center + 0.5 * length * Vec2d::UnitY(), 0.5 * width, 100, ShapeBooleanOps::add);
    }
};

class VelocityBoundaryCondition : public BaseLocalDynamics<SPHBody>, public DataDelegateSimple
{
  private:
    StdLargeVec<Vec2d> *vel_;
    StdLargeVec<Vec2d> *pos_;
    Vec2d velocity_;
    std::function<bool(SPH::Vec2d)> contains_func_;

  public:
    VelocityBoundaryCondition(SPHBody &body, Vec2d velocity, std::function<bool(SPH::Vec2d)> contains)
        : BaseLocalDynamics<SPHBody>(body), DataDelegateSimple(body),
          vel_(this->particles_->template registerSharedVariable<Vec2d>("Velocity")),
          pos_(this->particles_->template registerSharedVariable<Vec2d>("Position")),
          velocity_(std::move(velocity)),
          contains_func_(std::move(contains)){};
    inline void update(size_t index_i, [[maybe_unused]] Real dt = 0.0)
    {
        if (contains_func_((*pos_)[index_i]))
            (*vel_)[index_i] = velocity_;
    }
};

class ArcBoundaryGeometry : public BodyPartByParticle
{
  private:
    Real dp_;
    Vec2d center_;

  public:
    ArcBoundaryGeometry(SPHBody &body, const std::string &body_part_name, Real dp, const Vec2d &center)
        : BodyPartByParticle(body, body_part_name),
          dp_(dp), center_(center)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ArcBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (pos_[index_i].y() < center_.y() + dp_)
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
    // beta: shape constant (0.4 for arc)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

Real get_E_from_K_and_G(Real K, Real G) { return 9.0 * K * G / (3.0 * K + G); }
Real get_nu_from_K_and_G(Real K, Real G) { return (3 * K - 2 * G) / (2.0 * (3 * K + G)); }
Real get_E_from_K_and_nu(Real K, Real nu) { return 3.0 * K * (1.0 - 2.0 * nu); }
Real get_nu_from_K_and_E(Real K, Real E) { return (3 * K - E) / (6.0 * K); }
//------------------------------------------------------------------------------
void soft_stiff_contact(int dp_factor, bool use_pseudo_stiffness, double nu_min, int stiff_E_factor)
{
    // geometry
    const Real arc_inner_radius = 60;
    const Real arc_thickness = 2;
    const Real tube_width = 3;
    const Real tube_length = 120;

    // resolution
    const Real dp = arc_thickness / (4.0 * dp_factor);

    // time
    const Real max_disp = 100;
    const Vec2d velocity = 5 * Vec2d::UnitY();
    const Real end_time = max_disp / velocity.norm();

    // material properties
    const Real rho_arc = 1e3 * std::pow(unit_mm, 2); // 1e3 kg/m^3
    const Real bulk_modulus_arc = 3.33;              // 3.33 MPa
    const Real shear_modulus_arc = 0.34;             // 0.34 MPa
    const Real youngs_modulus_arc = get_E_from_K_and_G(bulk_modulus_arc, shear_modulus_arc);
    const Real poisson_ratio_arc = get_nu_from_K_and_G(bulk_modulus_arc, shear_modulus_arc);

    const Real rho_tube = 1e3 * std::pow(unit_mm, 2); // 1e3 kg/m^3
    const Real poisson_ratio_tube_real = 0.3;         // 0.4
    const Real youngs_modulus_tube_real = 1322;       // 1322 MPa
    const Real youngs_modulus_tube = use_pseudo_stiffness ? youngs_modulus_tube_real / stiff_E_factor : youngs_modulus_tube_real;
    const Real poisson_ratio_tube = use_pseudo_stiffness
                                        ? SMAX(nu_min, get_nu_from_K_and_E(bulk_modulus_arc, youngs_modulus_tube))
                                        : poisson_ratio_tube_real;

    // Material models
    auto material_arc = makeShared<NeoHookeanSolid>(rho_arc, youngs_modulus_arc, poisson_ratio_arc);
    auto material_tube = makeShared<SaintVenantKirchhoffSolid>(rho_tube, youngs_modulus_tube, poisson_ratio_tube);

    std::cout << "Youngs modulus of the arc: " << youngs_modulus_arc << std::endl;
    std::cout << "Bulk modulus of the arc: " << bulk_modulus_arc << std::endl;
    std::cout << "Poisson ratio of the arc: " << poisson_ratio_arc << std::endl;
    std::cout << "Sound speed of the arc: " << material_arc->ReferenceSoundSpeed() << std::endl;
    std::cout << "Youngs modulus of the tube: " << youngs_modulus_tube << std::endl;
    std::cout << "Bulk modulus of the tube: " << material_tube->BulkModulus() << std::endl;
    std::cout << "Poisson ratio of the tube: " << poisson_ratio_tube << std::endl;
    std::cout << "Sound speed of the tube: " << material_tube->ReferenceSoundSpeed() << std::endl;
    std::cout << "Impact speed: " << velocity.norm() << std::endl;

    const Real physical_viscosity_arc = get_physical_viscosity_general(rho_arc, youngs_modulus_arc, arc_thickness);
    const Real physical_viscosity_tube = get_physical_viscosity_general(rho_tube, youngs_modulus_tube, tube_width);

    // Import meshes
    const Vec2d arc_translation = Vec2d::Zero();
    auto mesh_arc = std::make_shared<Arc>("Arc", arc_translation, arc_inner_radius, arc_thickness);

    const Vec2d tube_translation = (arc_inner_radius - tube_width) * Vec2d::UnitX() - 0.5 * tube_length * Vec2d::UnitY();
    auto mesh_tube = std::make_shared<Tube>("tube", tube_translation, tube_width, tube_length);

    // System bounding box
    BoundingBox bb_system = union_bounding_box(mesh_arc->getBounds(), mesh_tube->getBounds());

    // System
    SPHSystem system(bb_system, dp);
    IOEnvironment io_environment(system);

    // Create objects
    SolidBody arc_body(system, mesh_arc);
    arc_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    arc_body.assignMaterial(material_arc.get());
    arc_body.generateParticles<BaseParticles, Lattice>();

    SolidBody tube_body(system, mesh_tube);
    tube_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    tube_body.assignMaterial(material_tube.get());
    tube_body.generateParticles<BaseParticles, Lattice>();

    // Inner relation
    InnerRelation inner_arc(arc_body);
    InnerRelation inner_tube(tube_body);

    // relax solid
    relax_solid(arc_body, inner_arc);
    relax_solid(tube_body, inner_tube);

    // Methods
    solid_algs arc_algs(inner_arc, physical_viscosity_arc, 0.5);
    solid_algs tube_algs(inner_tube, physical_viscosity_tube, 0.5);

    // Contact
    contact_algs contact_arc_tube(arc_body, {&tube_body});
    contact_algs contact_tube_arc(tube_body, {&arc_body});

    // Boundary conditions
    ArcBoundaryGeometry arc_fix_bc_part(arc_body, "arcClampingBC", dp, arc_translation);
    SimpleDynamics<FixBodyPartConstraint> arc_fix_bc(arc_fix_bc_part);

    SimpleDynamics<VelocityBoundaryCondition> tube_vel_bc(tube_body, velocity,
                                                          [&](const Vec2d &pos)
                                                          { return pos.y() < 0; });

    // Check
    auto check_nan = [&](BaseParticles &particles)
    {
        const auto &pos_ = particles.ParticlePositions();
        for (const auto &pos : pos_)
            if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                throw std::runtime_error("position has become nan");
    };

    // Output
    arc_body.getBaseParticles().addVariableToWrite<Vec2d>("RepulsionForce");
    tube_body.getBaseParticles().addVariableToWrite<Vec2d>("RepulsionForce");
    arc_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    tube_body.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    arc_algs.corrected_config();
    tube_algs.corrected_config();

    contact_arc_tube.config_update();
    contact_tube_arc.config_update();

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 50.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = std::min(arc_algs.time_step_size(), tube_algs.time_step_size());
    auto run_simulation = [&]()
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                contact_arc_tube.density_update();
                contact_tube_arc.density_update();
                contact_arc_tube.force_update();
                contact_tube_arc.force_update();

                dt = std::min(arc_algs.time_step_size(), tube_algs.time_step_size());
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                arc_algs.stress_relaxation_first(dt);
                tube_algs.stress_relaxation_first(dt);

                tube_vel_bc.exec();
                arc_fix_bc.exec();

                arc_algs.damping(dt);
                tube_algs.damping(dt);

                tube_vel_bc.exec();
                arc_fix_bc.exec();

                arc_algs.stress_relaxation_second(dt);
                tube_algs.stress_relaxation_second(dt);

                arc_body.updateCellLinkedList();
                tube_body.updateCellLinkedList();

                contact_arc_tube.config_update();
                contact_tube_arc.config_update();

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                { // checking if any position has become nan
                    check_nan(arc_body.getBaseParticles());
                    check_nan(tube_body.getBaseParticles());
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