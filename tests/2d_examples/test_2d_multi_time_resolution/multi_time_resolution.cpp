/**
 * @file 	multi_time_resolution.cpp
 * @brief 	contact of two solid bodies with different stiffness
 * @details This is the a case for test multi time resolution using dual loop
 * 			of two bodies in contact with significantly different time step sizes.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.

void solid_contact_multi_cycle(int res_factor_1, int res_factor_2, bool dual_loop = false);

int main(int ac, char *av[])
{
    solid_contact_multi_cycle(1, 1, true);
}

Real to_rad(Real angle) { return angle * Pi / 180; }

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

template <typename... Ts>
SPH::BoundingBox union_bounding_box(const SPH::BoundingBox &a, const SPH::BoundingBox &b, Ts &&...args)
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

class RectangularShapeGenerator : public MultiPolygonShape
{
  public:
    RectangularShapeGenerator(const std::string &shape_name, const Vec2d &center, const Vec2d &halfsize, Real angle = 0.0)
        : MultiPolygonShape(shape_name)
    {
        Transform transform2d(Rotation2d(-to_rad(angle)), center);
        multi_polygon_.addABox(transform2d, halfsize, ShapeBooleanOps::add);
    }
};

struct SolidAlgorithm
{
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>> position_damping_;
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> time_step_size_;

    SolidAlgorithm(InnerRelation &inner_relation, Real physical_viscosity)
        : corrected_configuration_(inner_relation),
          stress_relaxation_first_half_(inner_relation),
          stress_relaxation_second_half_(inner_relation),
          position_damping_(0.2, inner_relation, "Velocity", physical_viscosity),
          time_step_size_(inner_relation.getSPHBody())
    {
    }

    inline Real get_time_step_size() { return time_step_size_.exec(); }
    inline void correct_config() { corrected_configuration_.exec(); }
    inline void stress_first_half(Real dt) { stress_relaxation_first_half_.exec(dt); }
    inline void stress_second_half(Real dt) { stress_relaxation_second_half_.exec(dt); }
    inline void damping() { position_damping_.exec(); }

    inline void exec_all(Real dt)
    {
        stress_first_half(dt);
        damping();
        stress_second_half(dt);
    }
};

struct SolidContact
{
    SurfaceContactRelation contact_relation_;
    std::unique_ptr<InteractionDynamics<solid_dynamics::ContactDensitySummation>> contact_density_;
    std::unique_ptr<InteractionWithUpdate<solid_dynamics::ContactForce>> contact_forces_;

    explicit SolidContact(SolidBody &contact_to, SolidBody &contact_from)
        : contact_relation_(contact_to, {&contact_from})
    {
        contact_to.getBaseParticles().registerSharedVariable<Real>("RepulsionDensity");
        contact_from.getBaseParticles().registerSharedVariable<Real>("RepulsionDensity");
        contact_density_ = std::make_unique<InteractionDynamics<solid_dynamics::ContactDensitySummation>>(contact_relation_);
        contact_forces_ = std::make_unique<InteractionWithUpdate<solid_dynamics::ContactForce>>(contact_relation_);
    }
    inline void update_config() { contact_relation_.updateConfiguration(); }
    inline void exec()
    {
        contact_density_->exec();
        contact_forces_->exec();
    }
};

void solid_contact_multi_cycle(int res_factor_1, int res_factor_2, bool dual_loop)
{
    constexpr Real unit_mm = 1e-3; // mm
    constexpr Real scale = 1.0;

    // geometry properties
    const Real aorta_radius = 12.5 * scale;
    const Real aorta_half_length = 25 * scale;
    const Real aorta_thickness = 4 * scale;
    const Real catheter_radius = 2.5 * scale;
    const Real catheter_half_length = 10 * scale;
    const Real catheter_angle = 30;

    // resolution
    const Real res_ref = aorta_thickness / 4.0;
    const Real res_aorta = res_ref / Real(res_factor_1);
    const Real res_catheter = res_ref / Real(res_factor_2);

    // velocity
    Transform transform2d(Rotation2d(-to_rad(catheter_angle)));
    const Vec2d gravity_vector = transform2d.shiftBaseStationToFrame(2.0 * Vec2d::UnitX());
    std::cout << "gravity: " << gravity_vector.transpose() << std::endl;
    const Real end_time = 5.0;

    // material properties
    const Real rho0_s = 1000.0 * std::pow(unit_mm, 2);
    const Real poisson = 0.3;
    const Real Youngs_modulus_aorta = 0.2e6 * std::pow(unit_mm, 2);
    const Real Youngs_modulus_catheter = 1000e6 * std::pow(unit_mm, 2);

    // shape
    const Vec2d aorta_translation(aorta_radius, aorta_half_length);
    const Vec2d aorta_halfsize(0.5 * aorta_thickness, aorta_half_length);
    auto aorta_shape = std::make_shared<RectangularShapeGenerator>("Aorta", aorta_translation, aorta_halfsize);

    const Vec2d catheter_translation(0, catheter_half_length);
    const Vec2d catheter_halfsize(catheter_radius, catheter_half_length);
    auto catheter_shape = std::make_shared<RectangularShapeGenerator>("Catheter", catheter_translation, catheter_halfsize, catheter_angle);

    // System
    const auto bb_box = union_bounding_box(aorta_shape->getBounds(), catheter_shape->getBounds());
    SPHSystem sph_system(bb_box, res_aorta);
    IOEnvironment io_environment(sph_system);

    // Body
    SolidBody aorta(sph_system, aorta_shape);
    aorta.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus_aorta, poisson);
    aorta.generateParticles<BaseParticles, Lattice>();

    SolidBody catheter(sph_system, catheter_shape);
    catheter.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus_catheter, poisson);
    catheter.generateParticles<BaseParticles, Lattice>();

    // Inner relation
    InnerRelation aorta_inner(aorta);
    InnerRelation catheter_inner(catheter);

    // Solid solvers
    SolidAlgorithm aorta_algs(aorta_inner, get_physical_viscosity_general(rho0_s, Youngs_modulus_aorta, aorta_thickness));
    SolidAlgorithm catheter_algs(catheter_inner, get_physical_viscosity_general(rho0_s, Youngs_modulus_catheter, res_catheter));

    // Contact
    SolidContact aorta_contact(aorta, catheter);
    SolidContact catheter_contact(catheter, aorta);

    // Boundary conditions
    Gravity gravity(gravity_vector);
    SimpleDynamics<GravityForce> constant_gravity(catheter, gravity);

    // Output
    catheter.addBodyStateForRecording<Vec2d>("Velocity");
    aorta.addBodyStateForRecording<Vec2d>("Velocity");
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    body_states_recording.writeToFile(0);

    // Initialization
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    aorta_algs.correct_config();
    catheter_algs.correct_config();
    constant_gravity.exec();

    // Solver
    size_t output_iteration = 0;
    const Real num_outputs = 100;
    const Real output_period = end_time / num_outputs;

    auto run_single_simulation = [&]()
    {
        output_iteration = 0;
        const Real dt_solid_acoustic_ref =
            SMIN(aorta_algs.get_time_step_size(), catheter_algs.get_time_step_size());
        Real dt = 0;
        size_t ite = 0;

        GlobalStaticVariables::physical_time_ = 0;
        TickCount t1 = TickCount::now();
        TimeInterval interval;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_period)
            {
                if (ite % 100 == 0)
                    std::cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt
                              << "\n";

                aorta_contact.exec();
                catheter_contact.exec();

                // time step computation
                dt = SMIN(aorta_algs.get_time_step_size(), catheter_algs.get_time_step_size());
                { // check for extremely decreased dt
                    if (dt < dt_solid_acoustic_ref / 1e2)
                    {
                        std::cout << "dt = " << dt << " <<< dt_solid_acoustic_ref = " << dt_solid_acoustic_ref << "\n";
                        throw std::runtime_error("The time step decreased too much, stopping simulation\n");
                    }
                }

                // stress and damping
                aorta_algs.exec_all(dt);
                catheter_algs.exec_all(dt);

                // update
                aorta.updateCellLinkedList();
                catheter.updateCellLinkedList();
                aorta_contact.update_config();
                catheter_contact.update_config();

                // timestepping
                ++ite;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            // output
            ++output_iteration;
            TickCount t2 = TickCount::now();
            body_states_recording.writeToFile(output_iteration);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();
        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    StdLargeVec<Vec2d> contact_force_ave(aorta.getBaseParticles().total_real_particles_);
    auto reset_contact_force_ave = [&]()
    {
        particle_for(
            execution::ParallelPolicy(),
            IndexRange(0, aorta.getBaseParticles().total_real_particles_),
            [&](size_t index_i)
            { contact_force_ave[index_i] = Vec2d::Zero(); });
    };
    auto contact_force_ave_accumulate = [&](double dt)
    {
        const auto &contact_force = *aorta.getBaseParticles().getVariableByName<Vec2d>("RepulsionForce");
        particle_for(
            execution::ParallelPolicy(),
            IndexRange(0, aorta.getBaseParticles().total_real_particles_),
            [&](size_t index_i)
            { contact_force_ave[index_i] += contact_force[index_i] * dt; });
    };
    auto contact_force_ave_compute = [&](double dt_sum)
    {
        particle_for(
            execution::ParallelPolicy(),
            IndexRange(0, aorta.getBaseParticles().total_real_particles_),
            [&](size_t index_i)
            { contact_force_ave[index_i] /= dt_sum; });
    };
    auto run_dual_simulation = [&]()
    {
        output_iteration = 0;
        const Real dt_aorta_ref = aorta_algs.get_time_step_size();
        const Real dt_catheter_ref = catheter_algs.get_time_step_size();
        size_t ite = 0;
        double dt_aorta = 0;
        double dt_catheter = 0;

        GlobalStaticVariables::physical_time_ = 0;
        TickCount t1 = TickCount::now();
        TimeInterval interval;

        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_period)
            {
                auto force_prior = *aorta.getBaseParticles().getVariableByName<Vec2d>("ForcePrior");
                force_prior = contact_force_ave;

                dt_aorta = aorta_algs.get_time_step_size();
                if (dt_aorta < dt_aorta_ref / 1e2)
                    throw std::runtime_error("The aorta time step decreased too much, stopping simulation\n");

                aorta_algs.exec_all(dt_aorta);

                reset_contact_force_ave();
                Real dt_catheter_sum = 0.0;
                while (dt_catheter_sum < dt_aorta)
                {
                    catheter_contact.exec();
                    aorta_contact.exec();
                    contact_force_ave_accumulate(dt_catheter);

                    const Real dt_catheter_temp = catheter_algs.get_time_step_size();
                    if (dt_catheter_temp < dt_catheter_ref / 1e2)
                        throw std::runtime_error("The catheter time step decreased too much, stopping simulation\n");
                    dt_catheter = std::min(dt_catheter_temp, dt_aorta - dt_catheter_sum);
                    catheter_algs.exec_all(dt_catheter);

                    catheter.updateCellLinkedList();
                    catheter_contact.update_config();
                    aorta_contact.update_config();

                    dt_catheter_sum += dt_catheter;
                }
                contact_force_ave_compute(dt_catheter_sum);

                // timestepping
                ++ite;
                integration_time += dt_aorta;
                GlobalStaticVariables::physical_time_ += dt_aorta;

                if (ite % 10 == 0)
                    std::cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_
                              << "	dt_aorta: " << dt_aorta
                              << "	dt_catheter: " << dt_catheter << "\n";

                // update contact configurations
                aorta.updateCellLinkedList();
            }
            // output
            ++output_iteration;
            TickCount t2 = TickCount::now();
            body_states_recording.writeToFile(output_iteration);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };
    try
    { // try and catch errors, especially decreased timestep
        if (dual_loop)
            run_dual_simulation();
        else
            run_single_simulation();
    }
    catch (const std::exception &e)
    {
        std::cout << e.what();

        // output final configuration for debugging
        ++output_iteration;
        aorta.setNewlyUpdated();
        catheter.setNewlyUpdated();
        body_states_recording.writeToFile(output_iteration + 1);
    }
}
