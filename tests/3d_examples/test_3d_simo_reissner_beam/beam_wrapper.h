#include "simo_reissner_beam_dynamics.h"

struct BarGeometryParameters
{
    StdVec<Vec3d> pos_;
    StdVec<Vec3d> n0_;
    StdVec<Vec3d> b_n0_;
    StdVec<Real> width_;
    StdVec<Real> thickness_;
    Real dp_;
};

struct BarMaterialProperties
{
    StdVec<Mat3d> C_N_;
    StdVec<Mat3d> C_M_;
    StdVec<Real> rho_A_;
    StdVec<Mat3d> I_rho_l_;
};

namespace SPH
{
/** Define application dependent particle generator for thin structure. */
class BarParticleGenerator;
template <>
class ParticleGenerator<LinearParticles, BarParticleGenerator> : public ParticleGenerator<LinearParticles>
{
  private:
    BarGeometryParameters params_;

  public:
    ParticleGenerator(SPHBody &sph_body, LinearParticles &linear_particles, const BarGeometryParameters &params)
        : ParticleGenerator<LinearParticles>(sph_body, linear_particles), params_(params)
    {
        sph_body.getSPHAdaptation().getKernel()->reduceOnce();
    };
    void prepareGeometricData() override
    {
        // the beam and boundary
        for (size_t i = 0; i < params_.pos_.size(); i++)
        {
            addPositionAndVolumetricMeasure(params_.pos_[i], params_.dp_);
            addLineProperties(params_.n0_[i], params_.b_n0_[i], params_.thickness_[i], params_.width_[i]);
        }
    }
};
} // namespace SPH

struct bar_algorithm
{
    InteractionDynamics<slender_structure_dynamics::BarCorrectConfiguration> corrected_configuration;
    Dynamics1Level<slender_structure_dynamics::SimoReissnerRelaxationFirstHalf_v2> stress_relaxation_first_half;
    Dynamics1Level<slender_structure_dynamics::SimoReissnerRelaxationSecondHalf_v2> stress_relaxation_second_half;
    SimpleDynamics<slender_structure_dynamics::VelocityUpdate> velocity_update;
    ReduceDynamics<slender_structure_dynamics::BarAcousticTimeStepSize> computing_time_step_size;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> bar_position_damping;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> bar_rotation_damping;
    SimpleDynamics<slender_structure_dynamics::BeamInitialGeometry> initial_geometry;

    explicit bar_algorithm(InnerRelation &inner, Real physical_viscosity, bool hourglass_control = false)
        : corrected_configuration(inner),
          stress_relaxation_first_half(inner),
          stress_relaxation_second_half(inner),
          velocity_update(inner.getSPHBody()),
          computing_time_step_size(inner.getSPHBody()),
          bar_position_damping(0.5, inner, "Velocity", physical_viscosity),
          bar_rotation_damping(0.5, inner, "AngularVelocity", physical_viscosity),
          initial_geometry(inner) {};
};

struct BarParameters
{
    std::string body_name;
    std::shared_ptr<SaintVenantKirchhoffSolid> material;
    BarGeometryParameters geometry_params;
    BarMaterialProperties material_params;
    Real physical_viscosity;
    bool hourglass_control = false;
};

struct bar_object
{
    std::unique_ptr<SolidBody> bar_body;
    std::unique_ptr<InnerRelation> inner_relation;
    std::unique_ptr<bar_algorithm> alg;

    // Default move constructor and move assignment operator as noexcept
    bar_object(bar_object &&) noexcept = default;
    bar_object &operator=(bar_object &&) noexcept = default;
    bar_object(SPHSystem &system,
               BarParameters &bar_params)
        : bar_body(std::make_unique<SolidBody>(system, makeShared<DefaultShape>(bar_params.body_name)))
    {
        bar_body->defineMaterial<SaintVenantKirchhoffSolid>(*bar_params.material);
        bar_body->generateParticles<LinearParticles, BarParticleGenerator>(bar_params.geometry_params);

        register_variables(bar_params.material_params);

        inner_relation = std::make_unique<InnerRelation>(*bar_body);
        alg = std::make_unique<bar_algorithm>(
            *inner_relation,
            bar_params.physical_viscosity,
            bar_params.hourglass_control);
    }
    void update_relation() { inner_relation->updateConfiguration(); }
    void corrected_matrix() { alg->corrected_configuration.exec(); }
    Real get_time_step_size() { return alg->computing_time_step_size.exec(); }
    void stress_relaxation_first_half(Real dt) { alg->stress_relaxation_first_half.exec(dt); }
    void stress_relaxation_second_half(Real dt) { alg->stress_relaxation_second_half.exec(dt); }
    void update_velocity(Real dt) { alg->velocity_update.exec(dt); }
    void damping(Real dt)
    {
        alg->bar_position_damping.exec(dt);
        alg->bar_rotation_damping.exec(dt);
    }
    void compute_initial_geometry() { alg->initial_geometry.exec(); }

    // help functions
    void register_variables(BarMaterialProperties &material_params)
    {
        bar_body->getBaseParticles().registerStateVariableData<Mat3d>(
            "TensileConstitutiveMatrix", [&](size_t index_i) -> Matd
            { return material_params.C_N_[index_i]; });
        bar_body->getBaseParticles().registerStateVariableData<Mat3d>(
            "BendingConstitutiveMatrix", [&](size_t index_i) -> Matd
            { return material_params.C_M_[index_i]; });
        bar_body->getBaseParticles().registerStateVariableData<Mat3d>(
            "InertiaMatrix", [&](size_t index_i) -> Matd
            { return material_params.I_rho_l_[index_i]; });
        bar_body->getBaseParticles().registerStateVariableData<Real>(
            "RhoCrossSectionArea", [&](size_t index_i) -> Real
            { return material_params.rho_A_[index_i]; });
        bar_body->getBaseParticles().registerStateVariableData<Mat3d>(
            "TransformationFromInitialToCurrent", IdentityMatrix<Matd>::value);
    }

    // getters
    Mat3d *get_I_rho_l()
    {
        return bar_body->getBaseParticles().getVariableDataByName<Mat3d>("InertiaMatrix");
    }
    Mat3d *get_C_N()
    {
        return bar_body->getBaseParticles().getVariableDataByName<Mat3d>("TensileConstitutiveMatrix");
    }
    Mat3d *get_C_M()
    {
        return bar_body->getBaseParticles().getVariableDataByName<Mat3d>("BendingConstitutiveMatrix");
    }
    Real *get_A_rho()
    {
        return bar_body->getBaseParticles().getVariableDataByName<Real>("RhoCrossSectionArea");
    }
};

struct bar_simulation
{
    SPHSystem system;
    IOEnvironment io_environment;
    std::vector<bar_object> objects;
    std::function<void(size_t)> output_function = [&](size_t) {};
    std::function<void(Real)> acceleration_bc = [&](Real dt = 0) {};
    std::function<void(Real)> velocity_bc = [&](Real dt = 0) {};
    size_t output_number = 50;

    explicit bar_simulation(Real dp) : system(BoundingBox{}, dp),
                                       io_environment(system) {}

    void add_bar_object(BarParameters &bar_params)
    {
        objects.emplace_back(system, bar_params);
    }

    void initialize_system()
    {
        system.initializeSystemCellLinkedLists();
        system.initializeSystemConfigurations();
        for (auto &obj : objects)
        {
            obj.corrected_matrix();
            obj.compute_initial_geometry();
        }
    }

    Real get_time_step_size()
    {
        Real dt = std::numeric_limits<Real>::max();
        for (auto &obj : objects)
            dt = std::min(dt, obj.get_time_step_size());
        return dt;
    }

    void run_until(Real end_time, bool use_damping = true)
    {
        Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
        size_t ite = 0;
        size_t output_ite = 0;
        Real output_period = end_time / output_number;
        Real dt = 0.0;
        Real dt_ref = get_time_step_size();
        std::cout << "The reference time step size is " << dt_ref << std::endl;
        TickCount t1 = TickCount::now();
        TimeInterval interval;
        while (physical_time < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: " << physical_time << "	dt: " << dt << std::endl;
                }

                dt = get_time_step_size();
                if (dt < dt_ref / 1.1)
                {
                    output_function(output_ite++);
                    throw std::runtime_error("Time step size too small, simulation stopped.");
                }
                // dt = 0.25 * dt;
                for (auto &obj : objects)
                    obj.update_velocity(dt);
                velocity_bc(dt);
                if (use_damping)
                    for (auto &obj : objects)
                        obj.damping(dt);
                velocity_bc(dt);
                for (auto &obj : objects)
                    obj.stress_relaxation_first_half(dt);
                acceleration_bc(dt);
                for (auto &obj : objects)
                    obj.stress_relaxation_second_half(dt);
                velocity_bc(dt);

                ite++;
                integral_time += dt;
                physical_time += dt;
            }
            TickCount t2 = TickCount::now();
            output_ite++;
            output_function(output_ite);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    }
};

// Helpers
inline Real get_area_circular_shape(Real radius)
{
    return Pi * radius * radius;
}

inline Real get_area_rectangular_shape(Real width, Real height)
{
    return width * height;
}

inline Real get_moment_of_inertia_circular_shape(Real radius)
{
    return Pi * pow(radius, 4) / 4.0;
}

inline Real get_moment_of_inertia_rectangular_shape(Real width, Real height)
{
    return width * pow(height, 3) / 12.0;
}

size_t get_closest_particle(SPH::BaseParticles &particles, const SPH::Vec3d &point)
{
    const auto *pos = particles.getVariableDataByName<Vec3d>("Position");
    // get the closest particle to a point
    size_t id = 0;
    Real min_dist = INFINITY;
    for (size_t i = 0; i < particles.TotalRealParticles(); i++)
    {
        Real dist = (pos[i] - point).norm();
        if (dist < min_dist)
        {
            min_dist = dist;
            id = i;
        }
    }
    return id;
}

inline Real
get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}