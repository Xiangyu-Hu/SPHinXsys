/**
 * @file 	test_3d_shell_to_solid_coupling.cpp
 * @brief 	This is the case file for the test of solid-shell coupling.
 * @author  Weiyi Kong (Virtonomy), Xiangyu Hu
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

struct solid_mr_algs
{
    AdaptiveInnerRelation inner_relation;
    AdaptiveSplittingInnerRelation damping_inner;
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration;
    Dynamics1Level<solid_dynamics::Integration1stHalfCauchy> stress_relaxation_first_half;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half;
    SimpleDynamics<NormalDirectionFromBodyShape> initial_normal_direction;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> normal_direction;
    DampingWithRandomChoice<InteractionAdaptiveSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size;

    explicit solid_mr_algs(SolidBody &body, Real physical_viscosity)
        : inner_relation(body),
          damping_inner(body),
          corrected_configuration(inner_relation),
          stress_relaxation_first_half(inner_relation),
          stress_relaxation_second_half(inner_relation),
          initial_normal_direction(body),
          normal_direction(body),
          damping(0.5, damping_inner, "Velocity", physical_viscosity),
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
//-------Material property initialization-------------------------------------------------
class SolidMaterialInitialization : public MaterialIdInitialization
{
  private:
    std::function<bool(Vec3d &)> contain_;

  public:
    SolidMaterialInitialization(SolidBody &solid_body, std::function<bool(Vec3d &)> contain)
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
/**@class CouplingPart
 * @brief Find the overlapping particles and set is_coupled to 1.
 * The overlapping criteria is set to distance < factor * average_spacing.
 */
class CouplingPart : public BodyPartByParticle
{
  private:
    int *is_coupled_;
    std::function<bool(Vec3d &)> contain_;

  public:
    CouplingPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec3d &)> contain)
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
    SimpleDynamics<solid_dynamics::TotalWeightComputation> total_weight;
    SimpleDynamics<solid_dynamics::InterpolationSolidVelocityConstraint> vel_bc;

    solid_coupling_algs(SolidBody &body,
                        RealBodyVector contact_bodies,
                        std::function<bool(Vec3d &)> contain,
                        const std::vector<Real> &factors = {})
        : contact_relation(body, std::move(contact_bodies), factors),
          part(body, "CouplingPart", std::move(contain)),
          total_weight(part, contact_relation),
          vel_bc(part, contact_relation) {}

    void initialize_total_weight() { total_weight.exec(); }
    void exec() { vel_bc.exec(); }
};

struct shell_coupling_algs
{
    NearestNeighborContactRelation contact_relation;
    CouplingPart part;
    InteractionWithUpdate<solid_dynamics::InterpolationShellForceConstraint> force_bc;

    shell_coupling_algs(SolidBody &body,
                        RealBodyVector contact_bodies,
                        std::function<bool(Vec3d &)> contain,
                        const std::vector<Real> &factors = {})
        : contact_relation(body, std::move(contact_bodies), factors),
          part(body, "CouplingPart", std::move(contain)),
          force_bc(part, contact_relation) {}

    void exec() { force_bc.exec(); }
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
    std::function<bool(Vec3d &)> contain_;
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
    std::function<bool(Vec3d &)> contain_;
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

    void add_coupling_algorithm(RealBodyVector contact_bodies, std::function<bool(Vec3d &)> contain, const std::vector<Real> &factors)
    {
        coupling_algs_ = std::make_unique<solid_coupling_algs>(body_, contact_bodies, std::move(contain), factors);
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

    void add_coupling_algorithm(RealBodyVector contact_bodies, std::function<bool(Vec3d &)> contain, const std::vector<Real> &factors)
    {
        coupling_algs_ = std::make_unique<shell_coupling_algs>(body_, contact_bodies, std::move(contain), factors);
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
    std::vector<std::function<bool(Vec3d &)>> solid_coupling_parts_;
    std::vector<std::function<bool(Vec3d &)>> shell_coupling_parts_;

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
        solid_coupling_parts_.emplace_back(std::move(input.contain_));
    };

    template <typename MaterialType>
    void add_shell_object(shell_input_variables &input)
    {
        shell_objects_.emplace_back(std::make_unique<shell_object>(system, input));
        shell_objects_.back()->particle_generation<MaterialType>(input);
        shell_coupling_parts_.emplace_back(std::move(input.contain_));
    };

    void add_bcs(Real max_elongation, Real speed, int problem_type)
    {
        for (auto &solid : solid_objects_)
            solid->add_bcs(max_elongation, speed, problem_type);
        for (auto &shell : shell_objects_)
            shell->add_bcs(max_elongation, speed, problem_type);
    }

    void add_coupling_algs(bool use_interpolation)
    {
        RealBodyVector solid_bodies;
        for (auto &solid : solid_objects_)
            solid_bodies.emplace_back(&solid->body_);
        RealBodyVector shell_bodies;
        for (auto &shell : shell_objects_)
            shell_bodies.emplace_back(&shell->body_);

        for (size_t i = 0; i < solid_objects_.size(); i++)
            use_interpolation ? solid_objects_[i]->add_coupling_algorithm(shell_bodies, solid_coupling_parts_[i], std::vector<Real>(shell_bodies.size(), 2.3))
                              : solid_objects_[i]->add_coupling_algorithm(shell_bodies, solid_coupling_parts_[i], {});

        for (size_t i = 0; i < shell_objects_.size(); i++)
            use_interpolation ? shell_objects_[i]->add_coupling_algorithm(solid_bodies, shell_coupling_parts_[i], std::vector<Real>(solid_bodies.size(), 2.3))
                              : shell_objects_[i]->add_coupling_algorithm(solid_bodies, shell_coupling_parts_[i], {});
    }

    void total_weight_initialization()
    {
        for (auto &solid : solid_objects_)
            solid->coupling_algs_->initialize_total_weight();
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

    Real get_time_step_size() const
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
};

template <typename material_type>
struct one_sided_problem
{
    plate_parameters params;

    inline auto get_solid_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, params.height + 2 * dp);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    }

    inline auto get_upper_shell_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, dp);
        const Vec3d translation = 0.5 * (params.height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "upper_shell");
    }

    inline auto get_lower_shell_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, dp);
        const Vec3d translation = -0.5 * (params.height + dp) * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "lower_shell");
    }

    std::pair<std::vector<solid_input_variables>, std::vector<shell_input_variables>> operator()(Real dp_solid, Real dp_shell, Real stiffness_ratio)
    {
        std::vector<solid_input_variables> solid_inputs;
        std::vector<shell_input_variables> shell_inputs;
        Real youngs_modulus_solid = params.youngs_modulus_shell * stiffness_ratio;
        Real eta = get_physical_viscosity_general(params.rho, params.youngs_modulus_shell, params.height + 2 * params.thickness_shell);
        solid_input_variables solid_input{
            .name_ = "solid",
            .dp_ = dp_solid,
            .mesh_ = get_solid_mesh(dp_solid),
            .material_ = std::make_shared<material_type>(params.rho, youngs_modulus_solid, params.poisson_ratio),
            .physical_viscosity_ = eta,
        };
        solid_input.contain_ = [z0 = solid_input.mesh_->getBounds().first_.z(),
                                z1 = solid_input.mesh_->getBounds().second_.z(),
                                dp_solid](Vec3d &pos)
        { return pos.z() < z0 + 0.7 * dp_solid || pos.z() > z1 - 0.7 * dp_solid; };
        solid_inputs.emplace_back(solid_input);

        shell_input_variables upper_shell_input{.name_ = "upper_shell",
                                                .dp_ = dp_shell,
                                                .thickness_ = params.thickness_shell,
                                                .mesh_ = get_upper_shell_mesh(dp_solid),
                                                .material_ = std::make_shared<material_type>(params.rho, params.youngs_modulus_shell, params.poisson_ratio),
                                                .physical_viscosity_ = eta,
                                                .contain_ = [](Vec3d &)
                                                { return true; }};
        shell_input_variables lower_shell_input{.name_ = "lower_shell",
                                                .dp_ = dp_shell,
                                                .thickness_ = params.thickness_shell,
                                                .mesh_ = get_lower_shell_mesh(dp_solid),
                                                .material_ = std::make_shared<material_type>(params.rho, params.youngs_modulus_shell, params.poisson_ratio),
                                                .physical_viscosity_ = eta,
                                                .contain_ = [](Vec3d &)
                                                { return true; }};
        {
            Real z_upper = 0.5 * params.height + 0.5 * dp_solid;
            Real z_lower = -0.5 * params.height - 0.5 * dp_solid;
            Real x = -0.5 * params.length + 0.5 * dp_shell;
            while (x < 0.5 * params.length)
            {
                Real y = -0.5 * params.width + 0.5 * dp_shell;
                while (y < 0.5 * params.width)
                {
                    upper_shell_input.pos_.emplace_back(x, y, z_upper);
                    upper_shell_input.n_.emplace_back(Vec3d::UnitZ());
                    lower_shell_input.pos_.emplace_back(x, y, z_lower);
                    lower_shell_input.n_.emplace_back(Vec3d::UnitZ());
                    y += dp_shell;
                }
                x += dp_shell;
            }
        }
        shell_inputs.emplace_back(upper_shell_input);
        shell_inputs.emplace_back(lower_shell_input);
        return std::make_pair(solid_inputs, shell_inputs);
    }
};

template <typename material_type>
struct two_sided_problem
{
    plate_parameters params;

    inline auto get_upper_solid_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, 0.5 * params.height + dp);
        const Vec3d translation = 0.25 * params.height * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "upper_solid");
    }

    inline auto get_lower_solid_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, 0.5 * params.height + dp);
        const Vec3d translation = -0.25 * params.height * Vec3d::UnitZ();
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize, "lower_solid");
    }

    inline auto get_shell_mesh(Real dp) const
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, dp);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "shell");
    }

    std::pair<std::vector<solid_input_variables>, std::vector<shell_input_variables>> operator()(Real dp_solid, Real dp_shell, Real stiffness_ratio)
    {
        std::vector<solid_input_variables> solid_inputs;
        std::vector<shell_input_variables> shell_inputs;
        Real youngs_modulus_solid = params.youngs_modulus_shell * stiffness_ratio;
        Real eta = get_physical_viscosity_general(params.rho, params.youngs_modulus_shell, params.height + params.thickness_shell);
        solid_inputs.emplace_back(solid_input_variables{.name_ = "upper_solid",
                                                        .dp_ = dp_solid,
                                                        .mesh_ = get_upper_solid_mesh(dp_solid),
                                                        .material_ = std::make_shared<material_type>(params.rho, youngs_modulus_solid, params.poisson_ratio),
                                                        .physical_viscosity_ = eta,
                                                        .contain_ = [dp_solid](Vec3d &pos)
                                                        {
                                                            return abs(pos.z()) < 0.5 * dp_solid;
                                                        }});
        solid_inputs.emplace_back(solid_input_variables{.name_ = "lower_solid",
                                                        .dp_ = dp_solid,
                                                        .mesh_ = get_lower_solid_mesh(dp_solid),
                                                        .material_ = std::make_shared<material_type>(params.rho, youngs_modulus_solid, params.poisson_ratio),
                                                        .physical_viscosity_ = eta,
                                                        .contain_ = [dp_solid](Vec3d &pos)
                                                        {
                                                            return abs(pos.z()) < 0.5 * dp_solid;
                                                        }});
        shell_input_variables shell_input{.name_ = "shell",
                                          .dp_ = dp_shell,
                                          .thickness_ = params.thickness_shell,
                                          .mesh_ = get_shell_mesh(dp_solid),
                                          .material_ = std::make_shared<material_type>(params.rho, params.youngs_modulus_shell, params.poisson_ratio),
                                          .physical_viscosity_ = eta,
                                          .contain_ = [](Vec3d &)
                                          { return true; }};

        {
            Real z = 0;
            Real x = -0.5 * params.length + 0.5 * dp_shell;
            while (x < 0.5 * params.length)
            {
                Real y = -0.5 * params.width + 0.5 * dp_shell;
                while (y < 0.5 * params.width)
                {
                    shell_input.pos_.emplace_back(x, y, z);
                    shell_input.n_.emplace_back(Vec3d::UnitZ());
                    y += dp_shell;
                }
                x += dp_shell;
            }
        }
        shell_inputs.emplace_back(shell_input);
        return std::make_pair(solid_inputs, shell_inputs);
    }
};
//-------Observers-------------------------------------------------
class InterfaceTotalForce : public LocalDynamicsReduce<ReduceSum<Vecd>>
{
  private:
    Vecd *force_;
    int *is_coupled_;

  public:
    explicit InterfaceTotalForce(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body),
          force_(particles_->getVariableDataByName<Vecd>("Force")),
          is_coupled_(particles_->getVariableDataByName<int>("IsCoupled"))
    {
        quantity_name_ = "TotalForce";
    }
    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        return is_coupled_[index_i] ? force_[index_i] : Vecd::Zero();
    }
};

class InterfaceTotalForcePrior : public LocalDynamicsReduce<ReduceSum<Vecd>>
{
  private:
    Vecd *force_prior_;
    int *is_coupled_;

  public:
    explicit InterfaceTotalForcePrior(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          is_coupled_(particles_->getVariableDataByName<int>("IsCoupled"))
    {
        quantity_name_ = "TotalForcePrior";
    }
    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        return is_coupled_[index_i] ? force_prior_[index_i] : Vecd::Zero();
    }
};

class InterfaceTotalEnergy : public LocalDynamicsReduce<ReduceSum<Real>>
{
  private:
    Vecd *pos_;
    Vecd *pos0_;
    Vecd *force_;
    int *is_coupled_;

  public:
    explicit InterfaceTotalEnergy(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          pos0_(particles_->getVariableDataByName<Vecd>("InitialPosition")),
          force_(particles_->getVariableDataByName<Vecd>("Force")),
          is_coupled_(particles_->getVariableDataByName<int>("IsCoupled"))
    {
        quantity_name_ = "TotalEnergy";
    }
    Real reduce(size_t index_i, Real dt = 0.0)
    {
        Vec3d displacement = pos_[index_i] - pos0_[index_i];
        return is_coupled_[index_i] ? force_[index_i].dot(displacement) : 0.0;
    }
};