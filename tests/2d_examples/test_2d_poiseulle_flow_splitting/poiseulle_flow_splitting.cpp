#include "sphinxsys.h"
using namespace SPH;

void poiseulle_flow();
void poiseulle_flow_mr();

int main(int ac, char *av[])
{
    poiseulle_flow_mr();
}

struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(1.0), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Real u_max = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        Real vel = 1.5 * u_max * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        return vel * Vecd::UnitX();
    }
};

void poiseulle_flow()
{
    // geometry
    const Real diameter = 1;
    const Real length = 2 * diameter;

    // resolution
    const Real dp_ref = 0.05;
    const Real dp_min = 0.5 * dp_ref;
    const Real wall_thickness = 4 * dp_ref;
    const Real emitter_length = 4 * dp_ref;

    // material properties
    const Real rho0 = 1;
    const Real U_max = 1;
    const Real c0 = 10 * U_max;
    const Real Re = 100;
    const Real mu = rho0 * U_max * diameter / Re;

    // generate shape
    Vec2d translation = 0.5 * (length - emitter_length) * Vec2d::UnitX();
    Vec2d fluid_halfsize = 0.5 * Vec2d(length + emitter_length, diameter);
    auto fluid_shape = [&]()
    {
        MultiPolygon shape;
        shape.addABox(Transform(translation), fluid_halfsize, ShapeBooleanOps::add);
        return std::make_shared<MultiPolygonShape>(shape, "fluid_shape");
    }();
    auto wall_shape = [&]()
    {
        MultiPolygon shape;
        shape.addABox(Transform(translation), fluid_halfsize + Vec2d(0, wall_thickness), ShapeBooleanOps::add);
        shape.addABox(Transform(translation), fluid_halfsize, ShapeBooleanOps::sub);
        return std::make_shared<MultiPolygonShape>(shape, "wall_shape");
    }();
    auto bbox = wall_shape->getBounds();

    // system
    SPHSystem sph_system(bbox, dp_ref);
    sph_system.setIOEnvironment();

    // create bodies
    FluidBody water_block(sph_system, fluid_shape);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0, c0, mu);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineAdaptationRatios(1.15, dp_ref / dp_min);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    // define inner shape
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    // define dynamics
    // wall
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    // fluid
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWallCorrection> viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionCorrectedComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    // inlet and outlet
    // emitter
    Vec2d emitter_halfsize = Vec2d(0.5 * emitter_length, 0.5 * diameter);
    Vec2d emitter_translation = -0.5 * emitter_length * Vec2d::UnitX();
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);

    // inlet velocity
    Vec2d inlet_flow_buffer_halfsize = Vec2d(0.5 * emitter_length, 0.5 * diameter);
    Vec2d inlet_flow_buffer_translation = -0.5 * emitter_length * Vec2d::UnitX();
    BodyAlignedBoxByCell inlet_flow_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inlet_flow_buffer);

    // outlet
    Vec2d disposer_halfsize = Vec2d(0.5 * emitter_length, 0.6 * diameter);
    Vec2d disposer_translation = (length - 0.5 * emitter_length) * Vec2d::UnitX();
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer);

    // Particle sorting
    ParticleSorting particle_sorting(water_block);

    // output
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Vec2d>(water_block, "Velocity");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.writeToFile(0);

    // initialization
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();

    // set up
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 100.0;  /**< End time. */
    Real Output_Time = 1.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */

    // cpu time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    // main loop
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            kernel_correction_complex.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
}

// class used to reset Vol_max
class ParticleVolumeLimit : public LocalDynamics
{
  private:
    Real Vol_0_; // dp_ref^Dimensions
    Real *Vol_max_;

  public:
    explicit ParticleVolumeLimit(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          Vol_0_(1e2 * std::pow(sph_body.getSPHBodyResolutionRef(), Dimensions)),
          Vol_max_(particles_->registerStateVariable<Real>("MaximumVolume", [this](size_t i)
                                                           { return Vol_0_; })) {}
    void update(size_t index_i, Real)
    {
        Vol_max_[index_i] = Vol_0_;
    }
};

// class to assign particle spacing in refined regions
// TODO: use a shape and smoothing for particles near the shape
class ParticleSplittingRegion : public LocalDynamics
{
  private:
    Shape *shape_;      // a shape to define the refined region
    double vol_target_; // maximum volume of the refined region
    Vecd *pos_;
    Real *Vol_max_;

  public:
    ParticleSplittingRegion(SPHBody &sph_body, Shape &shape, double dp_target, double coefficient = 1.15)
        : LocalDynamics(sph_body),
          shape_(&shape),
          vol_target_(coefficient * std::pow(dp_target, Dimensions)),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          Vol_max_(particles_->getVariableDataByName<Real>("MaximumVolume")) {}

    void update(size_t index_i, Real)
    {
        if (shape_->checkContain(pos_[index_i]))
            Vol_max_[index_i] = vol_target_;
    }
};

// base class for particle splitting
// TODO: use extrapolation to get new pressure and velocity
class BaseParticleSplitting : public LocalDynamics
{
  private:
    Fluid *fluid_;
    ParticleBuffer<Base> *buffer_;

    virtual Real get_epsilon() = 0;         // distance to the center: dx * epsilon
    virtual Real get_alpha() = 0;           // new smoothing length: alpha * h
    virtual int get_splitting_number() = 0; // number of daugther particles
    virtual Real get_mass_ratio(int k) = 0; // mass_daughter / mass_parent
    virtual Vecd get_rel_pos(int k) = 0;    // distance to the center / radius

    Real *h_ratio_;
    Real *mass_;
    Real *Vol_;
    Real *rho_;
    Vecd *pos_;
    Real *Vol_max_;
    Real h_ref_;
    BulkParticles within_scope_;

  public:
    explicit BaseParticleSplitting(SPHBody &sph_body, ParticleBuffer<Base> &buffer)
        : LocalDynamics(sph_body),
          fluid_(dynamic_cast<Fluid *>(&particles_->getBaseMaterial())),
          buffer_(&buffer),
          h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          Vol_max_(particles_->getVariableDataByName<Real>("MaximumVolume")),
          h_ref_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
          within_scope_(particles_)
    {
        buffer_->checkParticlesReserved();
    }
    virtual ~BaseParticleSplitting() = default;

    void update(size_t index_i, Real)
    {
        if (Vol_[index_i] > Vol_max_[index_i] && within_scope_(index_i)) // avoid spitting for surface particles
        {
            static std::mutex mutex_switch_to_real_;
            std::lock_guard lock(mutex_switch_to_real_);
            buffer_->checkEnoughBuffer(*particles_);
            Real h_ratio_new = h_ratio_[index_i] / get_alpha();
            Real radius = std::pow(Vol_[index_i], 1 / Real(Dimensions)) * get_epsilon();
            // generate new particles
            for (int k = 1; k < get_splitting_number(); ++k)
            {
                // copy properties to the daughter particle
                size_t index_new = particles_->createRealParticleFrom(index_i);
                // update properties
                pos_[index_new] += radius * get_rel_pos(k);
                mass_[index_new] *= get_mass_ratio(k);
                Vol_[index_new] = mass_[index_new] / rho_[index_new];
                h_ratio_[index_new] = h_ratio_new;
            }
            // update the first particle wihtout generating new particles
            pos_[index_i] += radius * get_rel_pos(0);
            mass_[index_i] *= get_mass_ratio(0);
            h_ratio_[index_i] = h_ratio_new;
            Vol_[index_i] = mass_[index_i] / rho_[index_i];
            // No need to call unlock explicitly, lock_guard will handle it
        }
    }
};

// 2d cubic pattern
class ParticleSplittingCubic : public BaseParticleSplitting
{
  private:
    Real get_epsilon() override { return 0.35; }
    Real get_alpha() override { return 0.6; }
    int get_splitting_number() override { return 4; }
    Real get_mass_ratio(int k) override { return 1.0 / Real(get_splitting_number()); }
    Vec2d get_rel_pos(int k) override
    {
        double theta = (2 * k + 1) * Pi / 4;
        return Vec2d(cos(theta), sin(theta));
    }

  public:
    using BaseParticleSplitting::BaseParticleSplitting;
};

// 2d hexagonal+1 pattern
class ParticleSplittingHexagon : public BaseParticleSplitting
{
  private:
    Real get_epsilon() override { return 0.4; }
    Real get_alpha() override { return 0.9; }
    int get_splitting_number() override { return 7; }
    Real get_mass_ratio(int k) override { return 1.0 / Real(get_splitting_number()); }
    Vec2d get_rel_pos(int k) override
    {
        if (k == 0)
            return Vec2d::Zero();
        double theta = k * Pi / 3;
        return Vec2d(cos(theta), sin(theta));
    }

  public:
    using BaseParticleSplitting::BaseParticleSplitting;
};

class MassExchange : public LocalDynamics, public DataDelegateInner
{
  private:
    Fluid *fluid_;

    int *is_exchange_;
    int *closest_id_;
    Real *h_ratio_;
    Real *Vol_max_;
    Vecd *pos_;
    Vecd *vel_;
    Vecd *force_;
    Real *mass_;
    Real *Vol_;
    Real *p_;
    Real *rho_;
    Real *drho_dt_;
    Real h_ref_;

  public:
    explicit MassExchange(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          fluid_(dynamic_cast<Fluid *>(&particles_->getBaseMaterial())),
          is_exchange_(particles_->registerStateVariable<int>("IsExchange")),
          closest_id_(particles_->registerStateVariable<int>("ClosestID", [](size_t i)
                                                             { return -1; })),
          h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")),
          Vol_max_(particles_->getVariableDataByName<Real>("MaximumVolume")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          force_(particles_->getVariableDataByName<Vecd>("Force")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          p_(particles_->getVariableDataByName<Real>("Pressure")),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          drho_dt_(particles_->getVariableDataByName<Real>("DensityChangeRate")),
          h_ref_(sph_body_.sph_adaptation_->ReferenceSmoothingLength())
    {
    }

    void initialization(size_t index_i, Real)
    {
        is_exchange_[index_i] = 0;
        closest_id_[index_i] = -1;
    }

    void interaction(size_t index_i, Real)
    {
        Real mass_max_i = Vol_max_[index_i] * rho_[index_i];
        if (mass_[index_i] >= 0.95 * mass_max_i)
            return;
        Real mass_diff = mass_max_i - mass_[index_i];
        Real h_i = h_ref_ / h_ratio_[index_i];
        const auto &inner_neighborhood = inner_configuration_[index_i];
        int closest_id = -1;
        Real r_ij_min = std::numeric_limits<Real>::max();
        for (size_t n = 0; n < inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real h_j = h_ref_ / h_ratio_[index_j];
            Real h_ave = 0.5 * (h_i + h_j);
            Real r_ij2 = inner_neighborhood.r_ij_[n] * inner_neighborhood.r_ij_[n];
            if (r_ij2 > h_ave ||
                mass_[index_j] > 0.95 * mass_max_i ||
                mass_[index_j] <= mass_diff ||
                std::abs(Vol_max_[index_i] - Vol_max_[index_j]) > Eps)
                continue;
            if (r_ij2 < r_ij_min)
            {
                r_ij_min = r_ij2;
                closest_id = index_j;
            }
        }
        closest_id_[index_i] = closest_id;
        if (closest_id != -1)
            is_exchange_[index_i] = 1;
    }

    void update(size_t index_i, Real)
    {
        if (closest_id_[index_i] == -1)
            return;
        size_t index_j = closest_id_[index_i];
        if (closest_id_[index_j] == -1)
            return;
        if (!is_exchange_[index_i] || !is_exchange_[index_j] || size_t(closest_id_[index_j]) != index_i)
            return;
        if (mass_[index_i] < mass_[index_j] || (abs(mass_[index_i] - mass_[index_j]) < Eps && index_i < size_t(index_j)))
            return;
        Real mass_max_i = Vol_max_[index_i] * rho_[index_i];
        Real mass_diff = mass_max_i - mass_[index_i];
        if (abs(mass_[index_j] - mass_diff) < Eps)
            return;

        static std::mutex mutex_switch_to_real_;
        std::lock_guard lock(mutex_switch_to_real_);
        is_exchange_[index_i] = 0;
        is_exchange_[index_j] = 0;

        Real morig = mass_[index_i]; // might change to other ratio
        Real mi = morig + mass_diff;
        mass_[index_i] = mi;
        mass_[index_j] -= mass_diff;
        Vol_[index_j] = mass_[index_j] / rho_[index_j];

        // Only big one moves and changes props.
        Real mj = mass_diff / mass_[index_i];
        mi = morig / mass_[index_i];
        pos_[index_i] = mi * pos_[index_i] + mj * pos_[index_j];
        vel_[index_i] = mi * vel_[index_i] + mj * vel_[index_j];
        p_[index_i] = mi * p_[index_i] + mj * p_[index_j];
        rho_[index_i] = fluid_->DensityFromPressure(p_[index_i]);
        force_[index_i] = mi * force_[index_i] + mj * force_[index_j];
        drho_dt_[index_i] = mi * drho_dt_[index_i] + mj * drho_dt_[index_j];
        Vol_[index_i] = mass_[index_i] / rho_[index_i];
    }
};

class DensityRenormalization : public LocalDynamics,
                               public DataDelegateInner,
                               public DataDelegateContact
{
  private:
    Real *rho_;
    Real *rho_re_;
    Real *mass_;
    Real *Vol_;
    Real rho0_;
    StdVec<Real *> contact_Vol_;
    Kernel *kernel_;
    Real *h_ratio_;
    int *indicator_;
    bool isNearFreeSurface(size_t index_i)
    {
        bool is_near_surface = false;
        const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            if (indicator_[inner_neighborhood.j_[n]] == 1)
            {
                is_near_surface = true;
                break;
            }
        }
        return is_near_surface;
    }

  public:
    explicit DensityRenormalization(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          DataDelegateContact(contact_relation),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          rho_re_(particles_->registerStateVariable<Real>("DensityRenormalization")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          rho0_(sph_body_.base_material_->ReferenceDensity()),
          kernel_(sph_body_.sph_adaptation_->getKernel()),
          h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")),
          indicator_(particles_->getVariableDataByName<int>("Indicator"))
    {
        for (size_t k = 0; k != contact_configuration_.size(); ++k)
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }

    void interaction(size_t index_i, Real)
    {
        Real W0 = kernel_->W0(h_ratio_[index_i], ZeroVecd);
        Real mass_sum = mass_[index_i] * W0;
        Real Vol_sum = Vol_[index_i] * W0;
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            mass_sum += inner_neighborhood.W_ij_[n] * mass_[inner_neighborhood.j_[n]];
            Vol_sum += inner_neighborhood.W_ij_[n] * Vol_[inner_neighborhood.j_[n]];
        }
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real *contact_Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                Real W_ijV_j = contact_neighborhood.W_ij_[n] * contact_Vol_k[contact_neighborhood.j_[n]];
                mass_sum += rho0_ * W_ijV_j;
                Vol_sum += W_ijV_j;
            }
        }
        Real rho_re = mass_sum / Vol_sum;
        // if (isNearFreeSurface(index_i))
        // {
        //     if (rho_re < rho_[index_i])
        //         rho_re = rho_re + SMAX(Real(0), (rho_[index_i] - rho_re)) * rho0_ / rho_[index_i];
        // }
        rho_re_[index_i] = rho_re;
    }

    void update(size_t index_i, Real)
    {
        rho_[index_i] = rho_re_[index_i];
        Vol_[index_i] = mass_[index_i] / rho_[index_i];
    }
};

class OutputKernelSummation : public LocalDynamics,
                              public DataDelegateInner,
                              public DataDelegateContact
{
  private:
    Real *W_ijV_j_sum_;
    Vecd *dW_ijV_j_e_ij_sum_;
    Real *Vol_;
    StdVec<Real *> contact_Vol_;
    Kernel *kernel_;
    Real *h_ratio_;

  public:
    OutputKernelSummation(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          DataDelegateContact(contact_relation),
          W_ijV_j_sum_(particles_->registerStateVariable<Real>("KernelSummation")),
          dW_ijV_j_e_ij_sum_(particles_->registerStateVariable<Vecd>("KernelGradientSummation")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          kernel_(sph_body_.sph_adaptation_->getKernel()),
          h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio"))
    {
        for (size_t k = 0; k != contact_configuration_.size(); ++k)
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
    void update(size_t index_i, Real)
    {
        Real W0 = kernel_->W0(h_ratio_[index_i], ZeroVecd);
        Real W_ijV_j_sum = W0 * Vol_[index_i];
        Vecd dW_ijV_j_e_ij_sum(0);
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            W_ijV_j_sum += inner_neighborhood.W_ij_[n] * Vol_[index_j];
            dW_ijV_j_e_ij_sum += inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        }
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real *contact_Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                W_ijV_j_sum += contact_neighborhood.W_ij_[n] * contact_Vol_k[index_j];
                dW_ijV_j_e_ij_sum += contact_neighborhood.dW_ij_[n] * contact_Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            }
        }
        W_ijV_j_sum_[index_i] = W_ijV_j_sum;
        dW_ijV_j_e_ij_sum_[index_i] = dW_ijV_j_e_ij_sum;
    }
};

void poiseulle_flow_mr()
{
    // geometry
    const Real diameter = 1;
    const Real length = 2 * diameter;

    // resolution
    const Real dp_ref = 0.05;
    const Real dp_min = 0.5 * dp_ref;
    const int refinement_level = std::floor(std::log2(dp_ref / dp_min));
    const Real wall_thickness = 4 * dp_ref;
    const Real emitter_length = 4 * dp_ref;
    const Real buffer_length = 10 * dp_ref;

    // material properties
    const Real rho0 = 1;
    const Real U_max = 1;
    const Real c0 = 10 * U_max;
    const Real Re = 100;
    const Real mu = rho0 * U_max * diameter / Re;

    // generate shape
    Vec2d translation_all = 0.5 * (length - emitter_length) * Vec2d::UnitX();
    Vec2d fluid_halfsize = 0.5 * Vec2d(length + emitter_length, diameter);
    auto fluid_shape = [&]()
    {
        MultiPolygon shape;
        shape.addABox(Transform(translation_all), fluid_halfsize, ShapeBooleanOps::add);
        return std::make_shared<MultiPolygonShape>(shape, "fluid_shape");
    }();
    auto wall_shape = [&]()
    {
        MultiPolygon shape;
        shape.addABox(Transform(translation_all), fluid_halfsize + Vec2d(0, wall_thickness), ShapeBooleanOps::add);
        shape.addABox(Transform(translation_all), fluid_halfsize, ShapeBooleanOps::sub);
        return std::make_shared<MultiPolygonShape>(shape, "wall_shape");
    }();
    auto bbox = wall_shape->getBounds();

    // system
    SPHSystem sph_system(bbox, dp_ref);
    sph_system.setIOEnvironment();

    // create bodies
    FluidBody water_block(sph_system, fluid_shape);
    water_block.defineAdaptation<ParticleWithLocalRefinement>(1.3, 1.0, refinement_level);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0, c0, mu);
    ParticleBuffer<ReserveSizeFactor> particle_buffer(7);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineAdaptationRatios(1.15, dp_ref / dp_min);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    // define inner shape
    AdaptiveInnerRelation water_block_inner(water_block);
    AdaptiveContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    // define dynamics
    // wall
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    // fluid
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplexAdaptive> update_density_by_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<DensityRenormalization> density_renormalization(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWallCorrection> viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseTransportVelocityCorrectionComplex<AdaptiveResolution, NoLimiter, NoKernelCorrection, BulkParticles>>
        transport_velocity_correction(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    // inlet and outlet
    // emitter
    Vec2d emitter_halfsize = Vec2d(0.5 * emitter_length, 0.5 * diameter);
    Vec2d emitter_translation = -0.5 * emitter_length * Vec2d::UnitX();
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, particle_buffer);

    // inlet velocity
    Vec2d inlet_flow_buffer_halfsize = Vec2d(0.5 * emitter_length, 0.5 * diameter);
    Vec2d inlet_flow_buffer_translation = -0.5 * emitter_length * Vec2d::UnitX();
    BodyAlignedBoxByCell inlet_flow_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inlet_flow_buffer);

    // outlet
    Vec2d disposer_halfsize = Vec2d(0.5 * emitter_length, 0.6 * diameter);
    Vec2d disposer_translation = (length - 0.5 * emitter_length) * Vec2d::UnitX();
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer);

    // particle splitting
    auto refinement_shape = [&]()
    {
        Vec2d halfsize = 0.5 * Vec2d(0.5 * length, 0.7 * diameter);
        Vec2d translation = (0.5 * length + halfsize.x()) * Vec2d::UnitX();
        MultiPolygon shape;
        shape.addABox(Transform(translation), halfsize, ShapeBooleanOps::add);
        return std::make_shared<MultiPolygonShape>(shape, "refinement_shape");
    }();
    SimpleDynamics<ParticleVolumeLimit> vol_max_reset(water_block);
    SimpleDynamics<ParticleSplittingRegion> splitting_region(water_block, *refinement_shape, dp_min);
    SimpleDynamics<ParticleSplittingHexagon> particle_splitting(water_block, particle_buffer);
    Dynamics1Level<MassExchange> mass_exchange(water_block_inner);

    // Particle sorting
    ParticleSorting particle_sorting(water_block);

    // Check
    auto check_nan = [&](BaseParticles *particles)
    {
        Vecd *pos = particles->getVariableDataByName<Vecd>("Position");
        for (size_t index_i = 0; index_i < particles->TotalRealParticles(); ++index_i)
            if (std::isnan(pos[index_i][0]) || std::isnan(pos[index_i][1]))
                throw std::runtime_error("position has become nan");
    };

    // initialization
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();

    // output
    SimpleDynamics<OutputKernelSummation> kernel_summation(water_block_inner, water_wall_contact);
    kernel_summation.exec();
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "SmoothingLengthRatio");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Real>(water_block, "Mass");
    body_states_recording.addToWrite<Real>(water_block, "VolumetricMeasure");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(water_block, "MaximumVolume");
    body_states_recording.addToWrite<Vec2d>(water_block, "Velocity");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "KernelSummation");
    body_states_recording.addToWrite<Vecd>(water_block, "KernelGradientSummation");
    body_states_recording.writeToFile(0);

    // set up
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 2.0;    /**< End time. */
    Real Output_Time = 1.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */

    // cpu time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    // main loop
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            // update_density_by_summation.exec();
            // density_renormalization.exec();
            kernel_correction_complex.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);
                check_nan(&water_block.getBaseParticles());
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            // particle splitting
            vol_max_reset.exec();
            splitting_region.exec();
            particle_splitting.exec();
            // mass_exchange.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        kernel_summation.exec();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
}