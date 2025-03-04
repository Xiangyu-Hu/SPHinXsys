/**
 * @file 	test_3d_shell_to_solid_coupling.cpp
 * @brief 	This is the case file for the test of solid-shell coupling.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "base_data_type.h"
#include "large_data_containers.h"
#include "sphinxsys.h"

using namespace SPH;

//-------Relaxation for the solid body-------------------------------------------------
inline void relax_solid(BaseInnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(inner.getSPHBody());
    RelaxationStepInner relaxation_step_inner(inner);
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
//-------Algorithm for solid body-------------------------------------------------
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
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> velocity_damping;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> rotation_damping;
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
//-------Boundary condition for solid body-------------------------------------------------
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
//-------Shell particle generator-------------------------------------------------
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
//-------Algorithm for shell body-------------------------------------------------
struct solid_coupling_algs
{
    ContactRelation contact_relation;
    solid_dynamics::CouplingPart part;
    SimpleDynamics<solid_dynamics::NearestNeighborSolidVelocityConstraint> vel_bc_nearest_neighbor;
    SimpleDynamics<solid_dynamics::InterpolationSolidVelocityConstraint> vel_bc_interpolation;

    solid_coupling_algs(SolidBody &body, RealBodyVector contact_bodies)
        : contact_relation(body, std::move(contact_bodies)),
          part(contact_relation, "CouplingPart"),
          vel_bc_nearest_neighbor(part, contact_relation),
          vel_bc_interpolation(part, contact_relation) {}

    void init(bool use_interpolation) { use_interpolation ? vel_bc_interpolation.init() : vel_bc_nearest_neighbor.init(); }
    void exec(bool use_interpolation) { use_interpolation ? vel_bc_interpolation.exec() : vel_bc_nearest_neighbor.exec(); }
};

struct shell_coupling_algs
{
    ContactRelation contact_relation;
    solid_dynamics::CouplingPart part;
    InteractionWithUpdate<solid_dynamics::NearestNeighborShellForceConstraint> force_bc_nearest_neighbor;
    InteractionWithUpdate<solid_dynamics::InterpolationShellForceConstraint> force_bc_interpolation;

    shell_coupling_algs(SolidBody &body, RealBodyVector contact_bodies)
        : contact_relation(body, std::move(contact_bodies)),
          part(contact_relation, "CouplingPart"),
          force_bc_nearest_neighbor(part, contact_relation),
          force_bc_interpolation(part, contact_relation) {}

    void init(bool use_interpolation) { use_interpolation ? force_bc_interpolation.init() : force_bc_nearest_neighbor.init(); }
    void exec(bool use_interpolation) { use_interpolation ? force_bc_interpolation.exec() : force_bc_nearest_neighbor.exec(); }
};
//-------Auxiliary functions-------------------------------------------------
inline Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

template <typename... Ts>
inline SPH::BoundingBox union_bounding_box(const SPH::BoundingBox &a, const SPH::BoundingBox &b, Ts &&...args)
{
    SPH::BoundingBox out = a;
    out.first_[0] = std::min(a.first_[0], b.first_[0]);
    out.first_[1] = std::min(a.first_[1], b.first_[1]);
    out.first_[2] = std::min(a.first_[2], b.first_[2]);
    out.second_[0] = std::max(a.second_[0], b.second_[0]);
    out.second_[1] = std::max(a.second_[1], b.second_[1]);
    out.second_[2] = std::max(a.second_[2], b.second_[2]);
    if constexpr (sizeof...(args) > 0)
        return union_bounding_box(out, args...);
    else
        return out;
}
//-------Boundary condition for the plate-------------------------------------------------
// 0: elongation, 1: shear
struct boundary_condition
{
    Real speed_;
    Real elongation_time_;

    Vec3d direction_;

    SolidBodyPart part_left_;
    SolidBodyPart part_right_;
    SimpleDynamics<VelocityBoundaryCondition> left_bc_;
    SimpleDynamics<VelocityBoundaryCondition> right_bc_;

    Vec3d get_velocity(Real time, const Vec3d &direction) const
    {
        Real vel = time < elongation_time_ ? speed_ : 0;
        return vel * direction;
    }

    boundary_condition(SolidBody &body, Real max_elongation, Real speed, int problem_type)
        : speed_(speed),
          elongation_time_(0.5 * max_elongation / speed),
          direction_(problem_type == 0 ? Vec3d::UnitX() : Vec3d::UnitZ()),
          part_left_(body, "LeftPart", [x0 = body.getInitialShape().getBounds().first_.x(), dp = body.getSPHBodyResolutionRef()](Vec3d &pos)
                     { return pos.x() < x0 + dp; }),
          part_right_(body, "RightPart", [x1 = body.getInitialShape().getBounds().second_.x(), dp = body.getSPHBodyResolutionRef()](Vec3d &pos)
                      { return pos.x() > x1 - dp; }),
          left_bc_(part_left_, [this](Real time)
                   { return get_velocity(time, -direction_); }),
          right_bc_(part_right_,
                    [this](Real time)
                    { return get_velocity(time, direction_); }) {};

    void exec()
    {
        left_bc_.exec();
        right_bc_.exec();
    }
};
//-------simulation wrappers-------------------------------------------------
struct solid_input_variables
{
    std::string name_;
    Real dp_;
    std::shared_ptr<Shape> mesh_;
    std::shared_ptr<ElasticSolid> material_;
    Real physical_viscosity_;
};

struct shell_input_variables
{
    std::string name_;
    Real dp_;
    Real thickness_;
    std::shared_ptr<Shape> mesh_;
    std::shared_ptr<ElasticSolid> material_;
    StdLargeVec<Vecd> pos_;
    StdLargeVec<Vecd> n_;
    Real physical_viscosity_;
};

struct solid_object
{
    SolidBody body_;
    std::unique_ptr<solid_algs> algs_;
    std::unique_ptr<solid_coupling_algs> coupling_algs_;
    std::unique_ptr<boundary_condition> bcs_;

    explicit solid_object(SPHSystem &system, const solid_input_variables &input)
        : body_(system, input.mesh_, input.name_) {}

    template <typename MaterialType>
    void particle_generation(const solid_input_variables &input, bool relax = false)
    {
        body_.defineAdaptationRatios(1.15, body_.getSPHSystem().ReferenceResolution() / input.dp_);
        body_.defineBodyLevelSetShape()->cleanLevelSet(0);
        body_.defineMaterial<MaterialType>(*dynamic_cast<MaterialType *>(input.material_.get()));
        body_.generateParticles<BaseParticles, Lattice>();
        algs_ = std::make_unique<solid_algs>(body_, input.physical_viscosity_);

        if (relax)
            relax_solid(algs_->inner_relation);
    }

    void add_coupling_algorithm(RealBodyVector contact_bodies)
    {
        coupling_algs_ = std::make_unique<solid_coupling_algs>(body_, contact_bodies);
    }

    void add_bcs(Real max_elongation, Real speed, int problem_type)
    {
        bcs_ = std::make_unique<boundary_condition>(body_, max_elongation, speed, problem_type);
    }
};

struct shell_object
{
    SolidBody body_;
    std::unique_ptr<shell_algs> algs_;
    std::unique_ptr<shell_coupling_algs> coupling_algs_;
    std::unique_ptr<boundary_condition> bcs_;

    explicit shell_object(SPHSystem &system, const shell_input_variables &input)
        : body_(system, input.mesh_, input.name_) {};

    template <typename MaterialType>
    void particle_generation(const shell_input_variables &input)
    {
        body_.defineAdaptationRatios(1.15, body_.getSPHSystem().ReferenceResolution() / input.dp_);
        body_.defineMaterial<MaterialType>(*dynamic_cast<MaterialType *>(input.material_.get()));
        body_.generateParticles<SurfaceParticles, ShellDirectGenerator>(input.pos_, input.n_, input.dp_, input.thickness_);
        algs_ = std::make_unique<shell_algs>(body_, input.physical_viscosity_);
    }

    void add_coupling_algorithm(RealBodyVector contact_bodies)
    {
        coupling_algs_ = std::make_unique<shell_coupling_algs>(body_, contact_bodies);
    }

    void add_bcs(Real max_elongation, Real speed, int problem_type)
    {
        bcs_ = std::make_unique<boundary_condition>(body_, max_elongation, speed, problem_type);
    }
};

struct simulation_system
{
    SPHSystem system;
    IOEnvironment io_environment;
    std::vector<std::unique_ptr<solid_object>> solid_objects_;
    std::vector<std::unique_ptr<shell_object>> shell_objects_;

    explicit simulation_system(Real dp = 1)
        : system(BoundingBox{}, dp),
          io_environment(system) {}

    void system_initialize()
    {
        system.initializeSystemCellLinkedLists();
        system.initializeSystemConfigurations();
    }

    template <typename MaterialType>
    void add_solid_object(const solid_input_variables &input, bool relax = false)
    {
        system.system_domain_bounds_ = union_bounding_box(system.system_domain_bounds_, input.mesh_->getBounds());
        solid_objects_.emplace_back(std::make_unique<solid_object>(system, input));
        solid_objects_.back()->particle_generation<MaterialType>(input, relax);
    };

    template <typename MaterialType>
    void add_shell_object(shell_input_variables &input)
    {
        shell_objects_.emplace_back(std::make_unique<shell_object>(system, input));
        shell_objects_.back()->particle_generation<MaterialType>(input);
    };

    void add_bcs(Real max_elongation, Real speed, int problem_type)
    {
        for (auto &solid : solid_objects_)
            solid->add_bcs(max_elongation, speed, problem_type);
        for (auto &shell : shell_objects_)
            shell->add_bcs(max_elongation, speed, problem_type);
    }

    void add_coupling_algs()
    {
        RealBodyVector solid_bodies;
        for (auto &solid : solid_objects_)
            solid_bodies.emplace_back(&solid->body_);
        RealBodyVector shell_bodies;
        for (auto &shell : shell_objects_)
            shell_bodies.emplace_back(&shell->body_);

        for (auto &solid : solid_objects_)
            solid->add_coupling_algorithm(shell_bodies);
        for (auto &shell : shell_objects_)
            shell->add_coupling_algorithm(solid_bodies);
    }

    void initialise_coupling_algs(bool use_interpolation)
    {
        for (auto &solid : solid_objects_)
            solid->coupling_algs_->init(use_interpolation);
        for (auto &shell : shell_objects_)
            shell->coupling_algs_->init(use_interpolation);
    }

    void init_config()
    {
        for (auto &solid : solid_objects_)
            solid->algs_->corrected_config();
        for (auto &shell : shell_objects_)
            shell->algs_->corrected_config();
    }

    template <typename... Type>
    void add_variable_to_write(std::string name)
    {
        for (auto &solid : solid_objects_)
            solid->body_.getBaseParticles().addVariableToWrite<Type...>(name);
        for (auto &shell : shell_objects_)
            shell->body_.getBaseParticles().addVariableToWrite<Type...>(name);
    }

    Real get_time_step_size()
    {
        Real dt = std::numeric_limits<Real>::max();
        for (auto &solid : solid_objects_)
            dt = std::min(dt, solid->algs_->time_step_size());
        for (auto &shell : shell_objects_)
            dt = std::min(dt, shell->algs_->time_step_size());
        return dt;
    }
};
//-------Problem parameters-------------------------------------------------
struct plate_parameters
{
    const Real unit_mm = 1e-3;
    const Real scale = 1.0;

    // Geometry
    const Real length = 40 * scale;
    const Real width = length;
    const Real height = 20 * scale;
    const Real thickness_shell = 0.5 * scale;

    // Material properties
    const Real rho = 1000 * pow(unit_mm, 2);
    const Real youngs_modulus_shell = 1e9 * pow(unit_mm, 2); // 10 GPa
    const Real poisson_ratio = 0.3;

    // Loading
    const Real maximum_elongation = 0.2 * length; // 20% strain
    const Real speed = 10;

    inline auto get_solid_extended_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, height + 2 * dp);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    }

    inline auto get_whole_mesh() const
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, height + 2 * thickness_shell);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    }

    inline auto get_upper_shell_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, dp);
        const Vec3d translation = 0.5 * (height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "upper_shell");
    }

    inline auto get_lower_shell_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(length, width, dp);
        const Vec3d translation = -0.5 * (height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "lower_shell");
    }
};