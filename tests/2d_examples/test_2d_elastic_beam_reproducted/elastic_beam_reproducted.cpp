/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases for
 * 			understanding SPH method for fluid-shell-interaction simulation.
 * @author 	Weiyi Kong
 */

#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real DL = 500.0 * scale;                   /**< Tank length. */
Real DH = 200.1 * scale;                   /**< Tank height. */
Real Dam_L = 100.0 * scale;                /**< Water block width. */
Real Dam_H = 140.0 * scale;                /**< Water block height. */
Real Gate_width = 5.0 * scale;             /**< Width of the gate. */
Real Base_bottom_position = 79.0 * scale;  /**< Position of gate base. (In Y direction) */
Real resolution_gate = Gate_width / 4.0;   /**< Initial reference particle spacing of the gate. */
Real resolution_ref = 2 * resolution_gate; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0;            /**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / resolution_gate) * resolution_gate);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/** Define the corner points of the gate constrain. */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); /**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW);           /**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW);                        /**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position);              /**< Right bottom. */
// force
Real max_force = 0.5;
Real load_time = 0.028;
//----------------------------------------------------------------------
//	Material parameters of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 1100;  /**< Reference density of gate. */
Real poisson = 0.47; /**< Poisson ratio. */
Real Youngs_modulus = 7.8e6;
Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * Gate_width * Gate_width;
//----------------------------------------------------------------------
//	create a gate shape
//----------------------------------------------------------------------
class GateParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit GateParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real x = DL - Dam_L - 0.5 * resolution_gate;
        Real y = 0.5 * resolution_gate;
        while (y < Dam_H + BW)
        {
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_gate);
            initializeSurfaceProperties(Vec2d(-1, 0), Gate_width);
            y += resolution_gate;
        }
    }
};
//----------------------------------------------------------------------
// Create the gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape;
    gate_constraint_shape.push_back(ConstrainP_lb);
    gate_constraint_shape.push_back(ConstrainP_lt);
    gate_constraint_shape.push_back(ConstrainP_rt);
    gate_constraint_shape.push_back(ConstrainP_rb);
    gate_constraint_shape.push_back(ConstrainP_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(gate_constraint_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
// Apply force
//----------------------------------------------------------------------
class PressureOnGate : public solid_dynamics::LoadingForce, public SolidDataSimple
{
    StdLargeVec<Vecd> &n_;
    StdLargeVec<Vecd> &pos0_;
    Real max_force_;
    Real load_time_;
    Real fix_position_;

  public:
    PressureOnGate(SPHBody &sph_body)
        : solid_dynamics::LoadingForce(sph_body, "PressureForceOnShell"),
          SolidDataSimple(sph_body),
          n_(particles_->n_), pos0_(particles_->pos0_),
          max_force_(max_force), load_time_(load_time),
          fix_position_(Base_bottom_position) {}
    void update(size_t index_i, Real dt = 0.0)
    {
        Real force_target = pos0_[index_i].y() >= fix_position_ ? max_force_ * (Dam_H - pos0_[index_i].y()) / (Dam_H - fix_position_) : max_force_;
        double force_ave = GlobalStaticVariables::physical_time_ < load_time_ ? force_target / load_time_ * GlobalStaticVariables::physical_time_ : force_target;
        loading_force_[index_i] = force_ave * n_[index_i];
        ForcePrior::update(index_i, dt);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody gate(sph_system, makeShared<DefaultShape>("Gate"));
    gate.defineAdaptationRatios(1.15, resolution_ref / resolution_gate);
    gate.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<GateParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation gate_inner(gate);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Algorithms of Elastic dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> gate_stress_relaxation_first_half(gate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> gate_stress_relaxation_second_half(gate_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> gate_computing_time_step_size(gate);
    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> gate_constraint(gate_constraint_part);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> gate_update_normal(gate);
    SimpleDynamics<PressureOnGate> apply_pressure_to_particles(gate);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.2, gate_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.2, gate_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    gate.addBodyStateForRecording<Vecd>("PressureForceOnShell");
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.5;
    Real output_interval = end_time / 200.0;
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            gate_update_normal.exec();
            apply_pressure_to_particles.exec();

            dt_s = 0.5 * gate_computing_time_step_size.exec();
            gate_stress_relaxation_first_half.exec(dt_s);
            gate_constraint.exec();
            shell_position_damping.exec(dt_s);
            shell_rotation_damping.exec(dt_s);
            gate_constraint.exec();
            gate_stress_relaxation_second_half.exec(dt_s);

            integration_time += dt_s;
            GlobalStaticVariables::physical_time_ += dt_s;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;
        }
        TickCount t2 = TickCount::now();
        write_real_body_states_to_vtp.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
