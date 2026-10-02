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

struct SolidAlgorithm
{
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> position_damping_;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> time_step_size_;

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
    inline void damping(Real dt) { position_damping_.exec(dt); }
};

struct SolidContact
{
    SurfaceContactRelation contact_relation_;
    std::unique_ptr<InteractionDynamics<solid_dynamics::ContactFactorSummation>> contact_factor_;
    std::unique_ptr<InteractionWithUpdate<solid_dynamics::ContactForce>> contact_forces_;

    explicit SolidContact(SolidBody &contact_to, SolidBody &contact_from)
        : contact_relation_(contact_to, {&contact_from})
    {
        contact_to.getBaseParticles().registerStateVariable<Real>("RepulsionFactor");
        contact_from.getBaseParticles().registerStateVariable<Real>("RepulsionFactor");
        contact_factor_ = std::make_unique<InteractionDynamics<solid_dynamics::ContactFactorSummation>>(contact_relation_);
        contact_forces_ = std::make_unique<InteractionWithUpdate<solid_dynamics::ContactForce>>(contact_relation_);
    }
    inline void update_config() { contact_relation_.updateConfiguration(); }
    inline void update_factor() { contact_factor_->exec(); }
    inline void update_force() { contact_forces_->exec(); }
};

void relax_body(SPHBody &body, InnerRelation &inner_relation)
{
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepInner relaxation_step(inner_relation);

    random_particles.exec(0.25);
    relaxation_step.SurfaceBounding().exec();

    int ite_p = 0;
    while (ite_p < 2000)
    {
        relaxation_step.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
    }
}

class CatheterVelBC : public ForcePrior
{
  private:
    Vec2d targer_vel_;
    std::function<bool(SPH::Vec2d)> contains_func_ = nullptr;
    Vec2d *pos_;
    Vec2d *vel_;
    Real *mass_;

  public:
    CatheterVelBC(SPHBody &body, Vec2d target_vel, std::function<bool(SPH::Vec2d)> contains)
        : ForcePrior(body, "CatheterVelForce"),
          targer_vel_(std::move(target_vel)),
          contains_func_(std::move(contains)),
          pos_(body.getBaseParticles().getVariableDataByName<Vec2d>("Position")),
          vel_(body.getBaseParticles().registerStateVariable<Vec2d>("Velocity")),
          mass_(body.getBaseParticles().registerStateVariable<Real>("Mass")) {}
    inline void update(size_t index_i, Real dt = 0.0)
    {
        if (dt < std::numeric_limits<double>::epsilon())
            return;
        if (contains_func_(pos_[index_i]))
            current_force_[index_i] = (targer_vel_ - vel_[index_i]) / dt * mass_[index_i];
        else
            current_force_[index_i] = Vec2d::Zero();
        ForcePrior::update(index_i, dt);
    }
};

class FixParticles : public BodyPartByParticle
{
  private:
    std::function<bool(SPH::Vec2d)> contains_func_ = nullptr;

  public:
    FixParticles(SPHBody &body, const std::string &body_part_name, std::function<bool(SPH::Vec2d)> contains)
        : BodyPartByParticle(body, body_part_name),
          contains_func_(std::move(contains))
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixParticles::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (contains_func_(pos_[index_i]))
            body_part_particles_.push_back(index_i);
    };
};

void solid_contact_multi_cycle(int res_factor_1, int res_factor_2, bool dual_loop)
{
    constexpr Real unit_mm = 1e-3; // mm
    constexpr Real scale = 1.0;

    // geometry properties
    const Real aorta_radius = 12.5 * scale;
    const Real aorta_curve_radius = 25 * scale;
    const Real aorta_thickness = 4 * scale;
    const Real catheter_radius = 2.5 * scale;
    const Real catheter_half_length = 30 * scale;
    const Real insertion_length = 55 * scale;

    // resolution
    const Real res_ref = aorta_thickness / 4.0;
    const Real res_aorta = res_ref / Real(res_factor_1);
    const Real res_catheter = res_ref / Real(res_factor_2);

    // material properties
    const Real rho0_s = 1000.0 * std::pow(unit_mm, 2);
    const Real poisson_aorta = 0.45;
    const Real poisson_catheter = 0.4;
    const Real Youngs_modulus_aorta = 0.98;     // MPa
    const Real Youngs_modulus_catheter = 71205; // MPa

    auto aorta_material = std::make_shared<NeoHookeanSolid>(rho0_s, Youngs_modulus_aorta, poisson_aorta);
    auto catheter_material = std::make_shared<NeoHookeanSolid>(rho0_s, Youngs_modulus_catheter, poisson_catheter);

    // velocity
    const Vec2d velocity_vec = std::min(aorta_material->ReferenceSoundSpeed(), catheter_material->ReferenceSoundSpeed()) / 20.0 * Vec2d::UnitY();
    const Real end_time = insertion_length / velocity_vec.norm();

    std::cout << "Aorta sound speed: " << aorta_material->ReferenceSoundSpeed() << "\n";
    std::cout << "Catheter sound speed: " << catheter_material->ReferenceSoundSpeed() << "\n";
    std::cout << "Insertion speed: " << velocity_vec.norm() << "\n";
    std::cout << "End time: " << end_time << "\n";

    // shape
    auto aorta_mesh = [&]()
    {
        MultiPolygon aorta_shape{};
        Real half_length = aorta_curve_radius + aorta_radius + aorta_thickness;
        aorta_shape.addACircle(Vec2d::Zero(), aorta_curve_radius + aorta_radius + aorta_thickness, 100, ShapeBooleanOps::add);
        aorta_shape.addACircle(Vec2d::Zero(), aorta_curve_radius + aorta_radius, 100, ShapeBooleanOps::sub);
        aorta_shape.addABox(Transform(-half_length * Vec2d::UnitX()), half_length * Vec2d::Ones(), ShapeBooleanOps::sub);
        aorta_shape.addABox(Transform(-half_length * Vec2d::UnitY()), half_length * Vec2d::Ones(), ShapeBooleanOps::sub);
        return std::make_shared<MultiPolygonShape>(aorta_shape, "Aorta");
    }();

    auto catheter_mesh = [&]()
    {
        MultiPolygon catheter_shape{};
        const Vec2d halfsize(catheter_radius, catheter_half_length);
        const Vec2d translation(aorta_curve_radius + aorta_radius - catheter_radius - 2 * res_catheter, -catheter_half_length);
        catheter_shape.addABox(Transform(translation), halfsize, ShapeBooleanOps::add);
        const Vec2d tip_center = translation.x() * Vec2d::UnitX();
        catheter_shape.addACircle(tip_center, catheter_radius, 100, ShapeBooleanOps::add);
        return std::make_shared<MultiPolygonShape>(catheter_shape, "Catheter");
    }();

    // System
    const auto bb_box = union_bounding_box(aorta_mesh->getBounds(), catheter_mesh->getBounds());
    SPHSystem sph_system(bb_box, res_aorta);
    IOEnvironment io_environment(sph_system);

    // Body
    SolidBody aorta(sph_system, *aorta_mesh, "Aorta");
    aorta.defineBodyLevelSetShape();
    aorta.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus_aorta, poisson_aorta);
    aorta.generateParticles<BaseParticles, Lattice>();

    SolidBody catheter(sph_system, *catheter_mesh, "Catheter");
    catheter.defineAdaptationRatios(1.15, res_aorta / res_catheter);
    catheter.defineBodyLevelSetShape();
    catheter.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus_catheter, poisson_catheter);
    catheter.generateParticles<BaseParticles, Lattice>();

    // Inner relation
    InnerRelation aorta_inner(aorta);
    InnerRelation catheter_inner(catheter);

    // relaxation
    relax_body(aorta, aorta_inner);
    relax_body(catheter, catheter_inner);

    // Solid solvers
    SolidAlgorithm aorta_algs(aorta_inner, get_physical_viscosity_general(rho0_s, Youngs_modulus_aorta, aorta_thickness));
    SolidAlgorithm catheter_algs(catheter_inner, get_physical_viscosity_general(rho0_s, Youngs_modulus_catheter, res_catheter));

    // Contact
    SolidContact aorta_contact(aorta, catheter);
    SolidContact catheter_contact(catheter, aorta);

    // Boundary condition
    FixParticles fix_aorta_ids(aorta, "Aorta", [&](Vec2d pos)
                               { return pos.y() < 4 * res_aorta; });
    SimpleDynamics<FixBodyPartConstraint> fix_bc_aorta(fix_aorta_ids);
    SimpleDynamics<CatheterVelBC> catheter_vel_bc(catheter, velocity_vec, [&](Vec2d pos)
                                                  { return pos.y() < 0; });

    // Output
    catheter.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    aorta.getBaseParticles().addVariableToWrite<Vec2d>("Velocity");
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.writeToFile(0);

    // Initialization
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    aorta_algs.correct_config();
    catheter_algs.correct_config();

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

        Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
        TickCount t1 = TickCount::now();
        TimeInterval interval;
        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_period)
            {
                if (ite % 100 == 0)
                    std::cout << "N=" << ite << " Time: " << physical_time << "	dt: " << dt
                              << "\n";

                catheter_vel_bc.exec(dt);

                aorta_contact.update_factor();
                catheter_contact.update_factor();
                catheter_contact.update_force();
                aorta_contact.update_force();

                // time step computation
                dt = SMIN(aorta_algs.get_time_step_size(), catheter_algs.get_time_step_size());
                // check for extremely decreased dt
                if (dt < dt_solid_acoustic_ref / 1e2)
                {
                    std::cout << "dt = " << dt << " <<< dt_solid_acoustic_ref = " << dt_solid_acoustic_ref << "\n";
                    throw std::runtime_error("The time step decreased too much, stopping simulation\n");
                }

                // stress and damping
                aorta_algs.stress_first_half(dt);
                fix_bc_aorta.exec();
                aorta_algs.damping(dt);
                fix_bc_aorta.exec();
                aorta_algs.stress_second_half(dt);

                catheter_algs.stress_first_half(dt);
                catheter_algs.damping(dt);
                catheter_algs.stress_second_half(dt);

                // update
                aorta.updateCellLinkedList();
                catheter.updateCellLinkedList();
                aorta_contact.update_config();
                catheter_contact.update_config();

                // timestepping
                ++ite;
                integration_time += dt;
                physical_time += dt;
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

    auto run_dual_simulation = [&]()
    {
        output_iteration = 0;
        const Real dt_aorta_ref = aorta_algs.get_time_step_size();
        const Real dt_catheter_ref = catheter_algs.get_time_step_size();
        size_t ite = 0;
        double dt_aorta = 0;
        double dt_catheter = 0;

        Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
        TickCount t1 = TickCount::now();
        TimeInterval interval;

        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_period)
            {
                aorta_contact.update_factor();
                catheter_contact.update_factor();
                catheter_contact.update_force();
                aorta_contact.update_force();

                dt_aorta = aorta_algs.get_time_step_size();
                if (dt_aorta < dt_aorta_ref / 1e2)
                    throw std::runtime_error("The aorta time step decreased too much, stopping simulation\n");

                aorta_algs.stress_first_half(dt_aorta);
                fix_bc_aorta.exec();
                aorta_algs.damping(dt_aorta);
                fix_bc_aorta.exec();
                aorta_algs.stress_second_half(dt_aorta);

                Real dt_catheter_sum = 0.0;
                while (dt_catheter_sum < dt_aorta)
                {
                    catheter_vel_bc.exec(dt_catheter);

                    const Real dt_catheter_temp = catheter_algs.get_time_step_size();
                    if (dt_catheter_temp < dt_catheter_ref / 1e2)
                        throw std::runtime_error("The catheter time step decreased too much, stopping simulation\n");
                    dt_catheter = std::min(dt_catheter_temp, dt_aorta - dt_catheter_sum);

                    catheter_algs.stress_first_half(dt_catheter);
                    catheter_algs.damping(dt_catheter);
                    catheter_algs.stress_second_half(dt_catheter);

                    dt_catheter_sum += dt_catheter;
                }

                // update
                aorta.updateCellLinkedList();
                catheter.updateCellLinkedList();
                aorta_contact.update_config();
                catheter_contact.update_config();

                // timestepping
                ++ite;
                integration_time += dt_aorta;
                physical_time += dt_aorta;

                if (ite % 10 == 0)
                    std::cout << "N=" << ite << " Time: " << physical_time
                              << "	dt_aorta: " << dt_aorta
                              << "	dt_catheter: " << dt_catheter << "\n";
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
