/**
 * @file 	oilcooling_2d.cpp
 * @brief 	2D example to show the filling process of oil cooling system in a running electric motor.
 * @details This is the feasibility verification for the 3D model based on the experimental study of the team
 * 			of Tanguy Davin.Coding is based on the programm filling_tank.cpp of Dr. Xiangyu Hu at TUM.
 * @author 	Lirong Zhuang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DMO = 0.2;  /**< Outer diameter of motor hausing. */
Real DM = 0.196; /**< Diameter of motor hausing. */
Real DRO = 0.12; /**< Diameter of outer rotor. */
Real DRI = 0.04; /**< Diameter of inner rotor. */
Real WL = 0.03;  /**< Winding length. */
Real WH = 0.025; /**< Winding height. */
Real AG = 0.001; /**< Air-gap. */
Real LL = 0.004; /**< Inflow region length. */
Real LH = 0.012; /**< Inflows region height. */
Real DO = 0.02;  /**< Outflow diameter. */

Real RM = 0.5 * DM;  /**< Radius of motor hausing. */
Real RR = 0.5 * DRI; /**< Radius of rotor. */
Real Wnum = 12;      /**< Winding Number. */
Real angle_increment = 2 * Pi / Wnum;
Real WD = 0.5 * DRO + AG + 0.5 * WH;     /**< Distance from center point to the center of a winding. */
Real ZW = 0.5 * DRO + AG;                /**< Distance from center point to the site of a winding. */
int resolution_circle = 60;              /**<Approximate the circle as the number of sides of the polygon. */
Real resolution_ref = 0.00075;           /**< Initial reference particle spacing. */
Real BW = 0.5 * (DMO - DM);              /**< Extending width for wall boundary. */
Real OH = LH;                            /**< Outflows region height. */
Real Lnum = 5;                           /**< Inflows number. */
Real Lstart = (Wnum - 2 * Lnum + 2) / 4; /**< Start location of inlets. */
Real Lend = (Wnum + 2 * Lnum - 2) / 4;   /**< End location of inlets. */
Real inlet_height = RM + BW - LH;        /**< Inflow location height */
Vec2d inlet_halfsize = Vec2d(0.5 * LH, 0.5 * LL);
Real R_inlet = RM + BW - 0.5 * LH;
Vec2d inlet_translation = Vec2d(0, R_inlet);
Vec2d outlet_halfsize = Vec2d(BW + 0.01, 0.5 * DO);
Vec2d outlet_translation = Vec2d(0, -RM + 0.01);
Transform outlet_transform(Rotation2d(-(Pi / 2)), outlet_translation);
Transform inlet_transform(Rotation2d(-(Pi / 2)), inlet_translation);
BoundingBoxd system_domain_bounds(Vec2d(-BW - RM, -BW - RM), Vec2d(RM + BW, RM + BW));
Vecd center(0.0, 0.0);
// observer location
StdVec<Vecd> slot_locations = {
    Vecd(WD, 0),
    Vecd(WD *cos(angle_increment), WD *sin(angle_increment)),
    Vecd(WD *cos(2 * angle_increment), WD *sin(2 * angle_increment)),
    Vecd(WD *cos(3 * angle_increment), WD *sin(3 * angle_increment)),
    Vecd(WD *cos(4 * angle_increment), WD *sin(4 * angle_increment)),
    Vecd(WD *cos(5 * angle_increment), WD *sin(5 * angle_increment)),
    Vecd(WD *cos(6 * angle_increment), WD *sin(6 * angle_increment)),
    Vecd(WD *cos(7 * angle_increment), WD *sin(7 * angle_increment)),
    Vecd(WD *cos(8 * angle_increment), WD *sin(8 * angle_increment)),
    Vecd(WD *cos(9 * angle_increment), WD *sin(9 * angle_increment)),
    Vecd(WD *cos(10 * angle_increment), WD *sin(10 * angle_increment)),
    Vecd(WD *cos(11 * angle_increment), WD *sin(11 * angle_increment))};
StdVec<Vecd> tooth_locations = {
    Vecd(WD, 0.5 * WL),
    Vecd(WD *cos(angle_increment) - 0.5 * WL * sin(angle_increment),
         WD *sin(angle_increment) + 0.5 * WL * cos(angle_increment)),
    Vecd(WD *cos(2 * angle_increment) - 0.5 * WL * sin(2 * angle_increment),
         WD *sin(2 * angle_increment) + 0.5 * WL * cos(2 * angle_increment)),
    Vecd(WD *cos(3 * angle_increment) - 0.5 * WL * sin(3 * angle_increment),
         WD *sin(3 * angle_increment) + 0.5 * WL * cos(3 * angle_increment)),
    Vecd(WD *cos(4 * angle_increment) - 0.5 * WL * sin(4 * angle_increment),
         WD *sin(4 * angle_increment) + 0.5 * WL * cos(4 * angle_increment)),
    Vecd(WD *cos(5 * angle_increment) - 0.5 * WL * sin(5 * angle_increment),
         WD *sin(5 * angle_increment) + 0.5 * WL * cos(5 * angle_increment)),
    Vecd(WD *cos(6 * angle_increment) - 0.5 * WL * sin(6 * angle_increment),
         WD *sin(6 * angle_increment) + 0.5 * WL * cos(6 * angle_increment)),
    Vecd(WD *cos(7 * angle_increment) - 0.5 * WL * sin(7 * angle_increment),
         WD *sin(7 * angle_increment) + 0.5 * WL * cos(7 * angle_increment)),
    Vecd(WD *cos(8 * angle_increment) - 0.5 * WL * sin(8 * angle_increment),
         WD *sin(8 * angle_increment) + 0.5 * WL * cos(8 * angle_increment)),
    Vecd(WD *cos(9 * angle_increment) - 0.5 * WL * sin(9 * angle_increment),
         WD *sin(9 * angle_increment) + 0.5 * WL * cos(9 * angle_increment)),
    Vecd(WD *cos(10 * angle_increment) - 0.5 * WL * sin(10 * angle_increment),
         WD *sin(10 * angle_increment) + 0.5 * WL * cos(10 * angle_increment)),
    Vecd(WD *cos(11 * angle_increment) - 0.5 * WL * sin(11 * angle_increment),
         WD *sin(11 * angle_increment) + 0.5 * WL * cos(11 * angle_increment))};
StdVec<Vecd> yoke_locations = {
    Vecd(WD + 0.5 * WH, 0.25 * WL),
    Vecd(WD *cos(angle_increment) + 0.5 * WH * cos(angle_increment) - 0.25 * WL * sin(angle_increment),
         WD *sin(angle_increment) + 0.5 * WH * sin(angle_increment) + 0.25 * WL * cos(angle_increment)),
    Vecd(WD *cos(2 * angle_increment) + 0.5 * WH * cos(2 * angle_increment) - 0.25 * WL * sin(2 * angle_increment),
         WD *sin(2 * angle_increment) + 0.5 * WH * sin(2 * angle_increment) + 0.25 * WL * cos(2 * angle_increment)),
    Vecd(WD *cos(3 * angle_increment) + 0.5 * WH * cos(3 * angle_increment) - 0.25 * WL * sin(3 * angle_increment),
         WD *sin(3 * angle_increment) + 0.5 * WH * sin(3 * angle_increment) + 0.25 * WL * cos(3 * angle_increment)),
    Vecd(WD *cos(4 * angle_increment) + 0.5 * WH * cos(4 * angle_increment) - 0.25 * WL * sin(4 * angle_increment),
         WD *sin(4 * angle_increment) + 0.5 * WH * sin(4 * angle_increment) + 0.25 * WL * cos(4 * angle_increment)),
    Vecd(WD *cos(5 * angle_increment) + 0.5 * WH * cos(5 * angle_increment) - 0.25 * WL * sin(5 * angle_increment),
         WD *sin(5 * angle_increment) + 0.5 * WH * sin(5 * angle_increment) + 0.25 * WL * cos(5 * angle_increment)),
    Vecd(WD *cos(6 * angle_increment) + 0.5 * WH * cos(6 * angle_increment) - 0.25 * WL * sin(6 * angle_increment),
         WD *sin(6 * angle_increment) + 0.5 * WH * sin(6 * angle_increment) + 0.25 * WL * cos(6 * angle_increment)),
    Vecd(WD *cos(7 * angle_increment) + 0.5 * WH * cos(7 * angle_increment) - 0.25 * WL * sin(7 * angle_increment),
         WD *sin(7 * angle_increment) + 0.5 * WH * sin(7 * angle_increment) + 0.25 * WL * cos(7 * angle_increment)),
    Vecd(WD *cos(8 * angle_increment) + 0.5 * WH * cos(8 * angle_increment) - 0.25 * WL * sin(8 * angle_increment),
         WD *sin(8 * angle_increment) + 0.5 * WH * sin(8 * angle_increment) + 0.25 * WL * cos(8 * angle_increment)),
    Vecd(WD *cos(9 * angle_increment) + 0.5 * WH * cos(9 * angle_increment) - 0.25 * WL * sin(9 * angle_increment),
         WD *sin(9 * angle_increment) + 0.5 * WH * sin(9 * angle_increment) + 0.25 * WL * cos(9 * angle_increment)),
    Vecd(WD *cos(10 * angle_increment) + 0.5 * WH * cos(10 * angle_increment) - 0.25 * WL * sin(10 * angle_increment),
         WD *sin(10 * angle_increment) + 0.5 * WH * sin(10 * angle_increment) + 0.25 * WL * cos(10 * angle_increment)),
    Vecd(WD *cos(11 * angle_increment) + 0.5 * WH * cos(11 * angle_increment) - 0.25 * WL * sin(11 * angle_increment),
         WD *sin(11 * angle_increment) + 0.5 * WH * sin(11 * angle_increment) + 0.25 * WL * cos(11 * angle_increment))};

Real rho0_f = 945;     /**< Reference density of fluid. */
Real gravity_g = 9.81; /**< Gravity force of fluid. */
std::string temperature_species_name = "Phi";
// dynamics informations of oil
Real flow_rate = 30;                                        /**< Oil flow rate [L / h] */
Real v_inlet = (flow_rate * 0.001 / 3600) / (Pi * LL * LL); /**< Inflow vilocity [m / s]. */
// dynamics informations of rotor
Real rotor_rotation_velocity = 150;                    /**<Angular velocity rpm. */
Real Omega = -(rotor_rotation_velocity * 2 * Pi / 60); /**<Angle of rotor. */
Real U_f = 2 * Pi * RR * rotor_rotation_velocity / 60; /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                 /**< Reference sound speed. */
// thermal parameters
Real mu_f = 0.00874;                                            /**< Dynamics viscosity [Pa * s]. */
Real phi_winding = 110.0;                                       /**< Temperature of winding at begin. */
Real phi_fluid_initial = 75.0;                                  /**< Temperature of oil at begin. */
Real k_oil = 1.0e-8;                                            /**< Diffusion coefficient of oil. */
Real k_winding = 1.0e-5;                                        /**< Diffusion coefficient of winding. */
Real k_contact = (2 * k_oil * k_winding) / (k_oil + k_winding); /**< Thermal conductivity between winding and oil. */
Real dq = 2.0;                                                  /**< Heating efficient of internal heat source [°C/s]. */
//----------------------------------------------------------------------
//	Geometrie of the othor 4 inlets.
//----------------------------------------------------------------------
Vec2d inlet2_translation = Vec2d(R_inlet * cos(5 * angle_increment), R_inlet *sin(5 * angle_increment));
Transform inlet2_transform(Rotation2d(-(angle_increment)), inlet2_translation);
Vec2d inlet3_translation = Vec2d(R_inlet * cos(4 * angle_increment), R_inlet *sin(4 * angle_increment));
Transform inlet3_transform(Rotation2d(-(2 * angle_increment)), inlet3_translation);
Vec2d inlet4_translation = Vec2d(R_inlet * cos(2 * angle_increment), R_inlet *sin(2 * angle_increment));
Transform inlet4_transform(Rotation2d(-(4 * angle_increment)), inlet4_translation);
Vec2d inlet5_translation = Vec2d(R_inlet * cos(angle_increment), R_inlet *sin(angle_increment));
Transform inlet5_transform(Rotation2d(-(5 * angle_increment)), inlet5_translation);
//----------------------------------------------------------------------
//	Case-dependent wall boundary
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, (RM + BW), resolution_circle, ShapeBooleanOps::add); /**< Outer wall of motor hausing. */
        multi_polygon_.addACircle(center, RM, resolution_circle, ShapeBooleanOps::sub);        /**< Inner wall of motor hausing. */
        multi_polygon_.addABox(inlet_transform, inlet_halfsize, ShapeBooleanOps::sub);         /**< Top Inlets. */
        multi_polygon_.addABox(outlet_transform, outlet_halfsize, ShapeBooleanOps::sub);       /**< Outlets. */
        multi_polygon_.addABox(inlet2_transform, inlet_halfsize, ShapeBooleanOps::sub);
        multi_polygon_.addABox(inlet3_transform, inlet_halfsize, ShapeBooleanOps::sub);
        multi_polygon_.addABox(inlet4_transform, inlet_halfsize, ShapeBooleanOps::sub);
        multi_polygon_.addABox(inlet5_transform, inlet_halfsize, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Case-dependent Rotor boundary
//----------------------------------------------------------------------
class RotorBoundary : public MultiPolygonShape
{
  public:
    explicit RotorBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, RR, resolution_circle, ShapeBooleanOps::add); /**< Rotor */
    }
};
//----------------------------------------------------------------------
//	Case-dependent Winding boundary
//----------------------------------------------------------------------
class WindingBoundary : public MultiPolygonShape
{
  public:
    explicit WindingBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Add the windings. */
        for (int i = 0; i < Wnum; ++i)
        {
            Real theta = i * angle_increment;
            Real center_x = WD * cos(theta);
            Real center_y = WD * sin(theta);
            Vec2d winding_translation(center_x, center_y);
            Vec2d winding_halfsize(WL / 2, WH / 2);
            Transform winding_transform(Rotation2d(theta - (Pi / 2)), winding_translation);
            multi_polygon_.addABox(winding_transform, winding_halfsize, ShapeBooleanOps::add);
        }
    }
};
//----------------------------------------------------------------------
//	Case-dependent Inletswater boundary
//----------------------------------------------------------------------
class FluidBoundary : public MultiPolygonShape
{
  public:
    explicit FluidBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(inlet_transform, inlet_halfsize, ShapeBooleanOps::add); /**< Top Inlets. */
        multi_polygon_.addABox(inlet2_transform, inlet_halfsize, ShapeBooleanOps::add);
        multi_polygon_.addABox(inlet3_transform, inlet_halfsize, ShapeBooleanOps::add);
        multi_polygon_.addABox(inlet4_transform, inlet_halfsize, ShapeBooleanOps::add);
        multi_polygon_.addABox(inlet5_transform, inlet_halfsize, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
  public:
    InletInflowCondition(AlignedBoxByParticle &aligned_box_part)
        : EmitterInflowCondition(aligned_box_part) {}

  protected:
    virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
    {
        return Vec2d(v_inlet, 0.0);
    }
};
//----------------------------------------------------------------------
//	Application dependent initial condition of winding.
//----------------------------------------------------------------------
class ThermoWindingInitialCondition : public LocalDynamics
{
  public:
    explicit ThermoWindingInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(temperature_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = phi_winding;
    };

  protected:
    Real *phi_;
};
//----------------------------------------------------------------------
//	The windings heat up due to the internal heat sources.
//----------------------------------------------------------------------
class ThermoWindingHeatSource : public LocalDynamics
{
  public:
    explicit ThermoWindingHeatSource(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(temperature_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] += dq * dt;
    }

  protected:
    Real *phi_;
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermofluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(temperature_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = phi_fluid_initial;
    };

  protected:
    Real *phi_;
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
using ThermalRelaxationComplex = DiffusionBodyRelaxationComplex<
    IsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Build up a SPHSystem */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody oil_body(sph_system, makeShared<FluidBoundary>("OilBody"));
    oil_body.getSPHAdaptation().resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    oil_body.defineClosure<WeaklyCompressibleFluid, Viscosity, IsotropicDiffusion>(
        ConstructArgs(rho0_f, c_f), mu_f, ConstructArgs(temperature_species_name, k_oil));
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(3500.0);
    oil_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall"));
    wall.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall.generateParticles<BaseParticles, Reload>(wall.getName())
        : wall.generateParticles<BaseParticles, Lattice>();

    SolidBody rotor(sph_system, makeShared<RotorBoundary>("Rotor"));
    rotor.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    rotor.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? rotor.generateParticles<BaseParticles, Reload>(rotor.getName())
        : rotor.generateParticles<BaseParticles, Lattice>();

    SolidBody winding(sph_system, makeShared<WindingBoundary>("Winding"));
    winding.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    winding.defineClosure<Solid, IsotropicDiffusion>(
        Solid(), ConstructArgs(temperature_species_name, k_winding));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? winding.generateParticles<BaseParticles, Reload>(winding.getName())
        : winding.generateParticles<BaseParticles, Lattice>();

    ObserverBody slot_observer(sph_system, "ObserverSlot");
    slot_observer.generateParticles<ObserverParticles>(slot_locations);
    ObserverBody tooth_observer(sph_system, "ObserverTooth");
    tooth_observer.generateParticles<ObserverParticles>(tooth_locations);
    ObserverBody yoke_observer(sph_system, "ObserverYoke");
    yoke_observer.generateParticles<ObserverParticles>(yoke_locations);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation wall_inner(wall);
        InnerRelation rotor_inner(rotor);
        InnerRelation winding_inner(winding);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_wall_particles(wall);
        SimpleDynamics<RandomizeParticlePosition> random_rotor_particles(rotor);
        SimpleDynamics<RandomizeParticlePosition> random_winding_particles(winding);
        RelaxationStepInner relaxation_step_inner_wall(wall_inner);
        RelaxationStepInner relaxation_step_inner_rotor(rotor_inner);
        RelaxationStepInner relaxation_step_inner_winding(winding_inner);
        BodyStatesRecordingToVtp write_wall_to_vtp(wall);
        BodyStatesRecordingToVtp write_rotor_to_vtp(rotor);
        BodyStatesRecordingToVtp write_winding_to_vtp(winding);
        ReloadParticleIO write_wall_particle_reload_files(wall);
        ReloadParticleIO write_rotor_particle_reload_files(rotor);
        ReloadParticleIO write_winding_particle_reload_files(winding);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_wall_particles.exec(0.25);
        random_rotor_particles.exec(0.25);
        random_winding_particles.exec(0.25);
        relaxation_step_inner_wall.SurfaceBounding().exec();
        relaxation_step_inner_rotor.SurfaceBounding().exec();
        relaxation_step_inner_winding.SurfaceBounding().exec();
        write_wall_to_vtp.writeToFile(0);
        write_rotor_to_vtp.writeToFile(0);
        write_winding_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner_wall.exec();
            relaxation_step_inner_rotor.exec();
            relaxation_step_inner_winding.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_wall_to_vtp.writeToFile(ite_p);
                write_rotor_to_vtp.writeToFile(ite_p);
                write_winding_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_wall_particle_reload_files.writeToFile(0);
        write_rotor_particle_reload_files.writeToFile(0);
        write_winding_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation oil_body_inner(oil_body);
    InnerRelation winding_inner(winding);
    ContactRelation oil_body_contact(oil_body, {&wall, &rotor, &winding});
    ContactRelation oil_contact(oil_body, {&winding});
    ContactRelation winding_contact(winding, {&oil_body});
    ContactRelation winding_slot_contact(slot_observer, {&winding});
    ContactRelation winding_tooth_contact(tooth_observer, {&winding});
    ContactRelation winding_yoke_contact(yoke_observer, {&winding});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation oil_body_complex(oil_body_inner, oil_body_contact);
    ComplexRelation oil_thermo_complex(oil_body_inner, oil_contact);
    ComplexRelation winding_thermo_complex(winding_inner, winding_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    SimpleDynamics<NormalDirectionFromBodyShape> rotor_normal_direction(rotor);
    SimpleDynamics<NormalDirectionFromBodyShape> winding_normal_direction(winding);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(oil_body_inner, oil_body_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(oil_body_inner, oil_body_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(oil_body_inner, oil_body_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(oil_body_inner, oil_body_contact);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(oil_body, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(oil_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(oil_body);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(oil_body_inner, oil_body_contact);

    AlignedBoxByParticle emitter(oil_body, AlignedBox(xAxis, inlet_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, inlet_buffer);
    AlignedBoxByParticle emitter2(oil_body, AlignedBox(xAxis, inlet2_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition2(emitter2);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection2(emitter2, inlet_buffer);
    AlignedBoxByParticle emitter3(oil_body, AlignedBox(xAxis, inlet3_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition3(emitter3);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection3(emitter3, inlet_buffer);
    AlignedBoxByParticle emitter4(oil_body, AlignedBox(xAxis, inlet4_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition4(emitter4);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection4(emitter4, inlet_buffer);
    AlignedBoxByParticle emitter5(oil_body, AlignedBox(xAxis, inlet5_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition5(emitter5);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection5(emitter5, inlet_buffer);
    AlignedBoxByCell outlet_disposer(oil_body, AlignedBox(xAxis, outlet_transform, outlet_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> outlet_disposer_outflow_deletion(outlet_disposer);

    GetDiffusionTimeStepSize get_thermal_time_step(oil_body);
    ThermalRelaxationComplex thermal_relaxation_complex_oil(oil_body_inner, oil_contact);
    ThermalRelaxationComplex thermal_relaxation_complex_winding(winding_inner, winding_contact);
    SimpleDynamics<ThermoWindingInitialCondition> thermowinding_condition(winding);
    SimpleDynamics<ThermofluidBodyInitialCondition> thermofluid_initial_condition(oil_body);
    SimpleDynamics<ThermoWindingHeatSource> heat_source(winding);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------

    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(oil_body, "Indicator");
    body_states_recording.addToWrite<Real>(oil_body, "Phi");
    body_states_recording.addToWrite<Real>(winding, "Phi");
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_slot_phi("Phi", winding_slot_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_tooth_phi("Phi", winding_tooth_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_yoke_phi("Phi", winding_yoke_contact);
    //----------------------------------------------------------------------
    //	Building of multibody for rotor rotation.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the rigid shell box. */
    SolidBodyPartForSimbody rotor_multibody(rotor, makeShared<RotorBoundary>("Rotor"));
    SimTK::Body::Rigid rigid_info(*rotor_multibody.body_part_mass_properties_);
    SimTK::MobilizedBody::Pin
        Rotor_Pin(matter.Ground(), SimTK::Transform(SimTKVec3(0)), rigid_info, SimTK::Transform(SimTKVec3(0)));
    /** Initial angle of rotation. */
    Rotor_Pin.setDefaultAngle(0.0);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    Rotor_Pin.setRate(state, Omega);
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.initialize(state);
    /** Coupling between SimBody and SPH.*/
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_rotor(rotor_multibody, MBsystem, Rotor_Pin, integ);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    rotor_normal_direction.exec();
    winding_normal_direction.exec();
    thermowinding_condition.exec();
    thermofluid_initial_condition.exec();
    indicate_free_surface.exec();
    constant_gravity.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    //----------------------------------------------------------------------
    //	Time stepping control parameters.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 200;
    Real end_time = 2.0;
    Real output_interval = 0.05;
    Real dt = 0.0; /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_slot_phi.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(SMIN(dt_thermal, get_fluid_time_step_size.exec()), Dt);

                integ.stepBy(dt);
                constraint_rotor.exec();

                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                inflow_condition2.exec();
                inflow_condition3.exec();
                inflow_condition4.exec();
                inflow_condition5.exec();
                density_relaxation.exec(dt);
                inflow_condition.exec();
                inflow_condition2.exec();
                inflow_condition3.exec();
                inflow_condition4.exec();
                inflow_condition5.exec();

                heat_source.exec(dt); /** Implement of heat recources. */
                thermal_relaxation_complex_winding.exec(dt);
                thermal_relaxation_complex_oil.exec(dt);

                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                TickCount current_time = TickCount::now();
                TimeInterval elapsed_time = current_time - t1 - interval;
                Real remaining_physical_time = end_time - physical_time;
                Real estimated_remaining_real_time = elapsed_time.seconds() *
                                                     (remaining_physical_time / physical_time);
                int hours = static_cast<int>(estimated_remaining_real_time) / 3600;
                int minutes = (static_cast<int>(estimated_remaining_real_time) % 3600) / 60;
                int seconds = static_cast<int>(estimated_remaining_real_time) % 60;
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << " Time = " << physical_time
                          << " Dt = " << Dt << " dt = " << dt
                          << " Remaining: " << hours << " h "
                          << minutes << " min " << seconds << " s\n";
            }
            number_of_iterations++;

            /** inflow emitter injection*/
            emitter_injection.exec();
            emitter_injection2.exec();
            emitter_injection3.exec();
            emitter_injection4.exec();
            emitter_injection5.exec();
            /** outflow delete*/
            outlet_disposer_outflow_deletion.exec();
            /** Update cell linked list and configuration. */
            oil_body.updateCellLinkedList();
            oil_body_complex.updateConfiguration();
            oil_thermo_complex.updateConfiguration();
            winding_thermo_complex.updateConfiguration();
            winding_slot_contact.updateConfiguration();
            winding_tooth_contact.updateConfiguration();
            winding_yoke_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        indicate_free_surface.exec();
        body_states_recording.writeToFile();
        write_slot_phi.writeToFile(number_of_iterations);
        write_tooth_phi.writeToFile(number_of_iterations);
        write_yoke_phi.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    // write_slot_phi.generateDataBase(0.005, 0.005);
    write_slot_phi.testResult();

    return 0;
}
