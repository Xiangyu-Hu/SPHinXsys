/**
 * @file 	oilcooling_2d.cpp
 * @brief 	2D example to show the filling process of oil cooling system in a running electric motor.
 * @details This is the sycl version.
 * @author 	sphinxsim@einsimo.com
 */
#include "all_simbody.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DMO = 0.2;      /**< Outer diameter of motor hausing. */
Real DM = 0.196;     /**< Diameter of motor hausing. */
Real DRO = 0.12;     /**< Diameter of outer rotor. */
Real DRI = 0.04;     /**< Diameter of inner rotor. */
Real WL = 0.03;      /**< Winding length. */
Real WH = 0.025;     /**< Winding height. */
Real AG = 0.001;     /**< Air-gap. */
Real LL = 0.004;     /**< Inflow region length. */
Real LH = 0.012;     /**< Inflows region height. */
Real DO = 0.02;      /**< Outflow diameter. */
Real RM = 0.5 * DM;  /**< Radius of motor hausing. */
Real RR = 0.5 * DRI; /**< Radius of rotor. */
Real Wnum = 12;      /**< Winding Number. */
Real angle_increment = 2 * Pi / Wnum;
Real WD = 0.5 * DRO + AG + 0.5 * WH; /**< Distance from center point to the center of a winding. */
Real ZW = 0.5 * DRO + AG;            /**< Distance from center point to the site of a winding. */
int resolution_circle = 60;          /**<Approximate the circle as the number of sides of the polygon. */
Real global_resolution = 0.00075;    /**< Initial reference particle spacing. */
Real BW = 0.5 * (DMO - DM);          /**< Extending width for wall boundary. */
BoundingBoxd system_domain_bounds(Vec2d(-BW - RM, -BW - RM), Vec2d(RM + BW, RM + BW));
//----------------------------------------------------------------------
//	Geometrie of the 5 inlets.
//----------------------------------------------------------------------
Real Lnum = 5;                           /**< Inflows number. */
Real Lstart = (Wnum - 2 * Lnum + 2) / 4; /**< Start location of inlets. */
Real Lend = (Wnum + 2 * Lnum - 2) / 4;   /**< End location of inlets. */
Real inlet_height = RM + BW - LH;        /**< Inflow location height */
Vec2d inlet_halfsize = Vec2d(0.5 * LH, 0.5 * LL);
Real R_inlet = RM + BW - 0.5 * LH;
Vec2d inlet1_translation = Vec2d(0, R_inlet);
Transform inlet1_transform(Rotation2d(-(Pi / 2)), inlet1_translation);
Vec2d inlet2_translation = Vec2d(R_inlet * cos(5 * angle_increment), R_inlet *sin(5 * angle_increment));
Transform inlet2_transform(Rotation2d(-(angle_increment)), inlet2_translation);
Vec2d inlet3_translation = Vec2d(R_inlet * cos(4 * angle_increment), R_inlet *sin(4 * angle_increment));
Transform inlet3_transform(Rotation2d(-(2 * angle_increment)), inlet3_translation);
Vec2d inlet4_translation = Vec2d(R_inlet * cos(2 * angle_increment), R_inlet *sin(2 * angle_increment));
Transform inlet4_transform(Rotation2d(-(4 * angle_increment)), inlet4_translation);
Vec2d inlet5_translation = Vec2d(R_inlet * cos(angle_increment), R_inlet *sin(angle_increment));
Transform inlet5_transform(Rotation2d(-(5 * angle_increment)), inlet5_translation);
//----------------------------------------------------------------------
//	Geometrie of the outlet.
//----------------------------------------------------------------------
Vec2d outlet_halfsize = Vec2d(BW + 0.01, 0.5 * DO);
Vec2d outlet_translation = Vec2d(0, -RM + 0.01);
Transform outlet_transform(Rotation2d(-(Pi / 2)), outlet_translation);
// observer location
StdVec<Vecd> slot_locations = {
    Vecd(WD, 0),
    Vecd(WD * cos(angle_increment), WD *sin(angle_increment)),
    Vecd(WD * cos(2 * angle_increment), WD *sin(2 * angle_increment)),
    Vecd(WD * cos(3 * angle_increment), WD *sin(3 * angle_increment)),
    Vecd(WD * cos(4 * angle_increment), WD *sin(4 * angle_increment)),
    Vecd(WD * cos(5 * angle_increment), WD *sin(5 * angle_increment)),
    Vecd(WD * cos(6 * angle_increment), WD *sin(6 * angle_increment)),
    Vecd(WD * cos(7 * angle_increment), WD *sin(7 * angle_increment)),
    Vecd(WD * cos(8 * angle_increment), WD *sin(8 * angle_increment)),
    Vecd(WD * cos(9 * angle_increment), WD *sin(9 * angle_increment)),
    Vecd(WD * cos(10 * angle_increment), WD *sin(10 * angle_increment)),
    Vecd(WD * cos(11 * angle_increment), WD *sin(11 * angle_increment))};
StdVec<Vecd> tooth_locations = {
    Vecd(WD, 0.5 * WL),
    Vecd(WD * cos(angle_increment) - 0.5 * WL * sin(angle_increment),
         WD *sin(angle_increment) + 0.5 * WL * cos(angle_increment)),
    Vecd(WD * cos(2 * angle_increment) - 0.5 * WL * sin(2 * angle_increment),
         WD *sin(2 * angle_increment) + 0.5 * WL * cos(2 * angle_increment)),
    Vecd(WD * cos(3 * angle_increment) - 0.5 * WL * sin(3 * angle_increment),
         WD *sin(3 * angle_increment) + 0.5 * WL * cos(3 * angle_increment)),
    Vecd(WD * cos(4 * angle_increment) - 0.5 * WL * sin(4 * angle_increment),
         WD *sin(4 * angle_increment) + 0.5 * WL * cos(4 * angle_increment)),
    Vecd(WD * cos(5 * angle_increment) - 0.5 * WL * sin(5 * angle_increment),
         WD *sin(5 * angle_increment) + 0.5 * WL * cos(5 * angle_increment)),
    Vecd(WD * cos(6 * angle_increment) - 0.5 * WL * sin(6 * angle_increment),
         WD *sin(6 * angle_increment) + 0.5 * WL * cos(6 * angle_increment)),
    Vecd(WD * cos(7 * angle_increment) - 0.5 * WL * sin(7 * angle_increment),
         WD *sin(7 * angle_increment) + 0.5 * WL * cos(7 * angle_increment)),
    Vecd(WD * cos(8 * angle_increment) - 0.5 * WL * sin(8 * angle_increment),
         WD *sin(8 * angle_increment) + 0.5 * WL * cos(8 * angle_increment)),
    Vecd(WD * cos(9 * angle_increment) - 0.5 * WL * sin(9 * angle_increment),
         WD *sin(9 * angle_increment) + 0.5 * WL * cos(9 * angle_increment)),
    Vecd(WD * cos(10 * angle_increment) - 0.5 * WL * sin(10 * angle_increment),
         WD *sin(10 * angle_increment) + 0.5 * WL * cos(10 * angle_increment)),
    Vecd(WD * cos(11 * angle_increment) - 0.5 * WL * sin(11 * angle_increment),
         WD *sin(11 * angle_increment) + 0.5 * WL * cos(11 * angle_increment))};
StdVec<Vecd> yoke_locations = {
    Vecd(WD + 0.5 * WH, 0.25 * WL),
    Vecd(WD * cos(angle_increment) + 0.5 * WH * cos(angle_increment) - 0.25 * WL * sin(angle_increment),
         WD *sin(angle_increment) + 0.5 * WH * sin(angle_increment) + 0.25 * WL * cos(angle_increment)),
    Vecd(WD * cos(2 * angle_increment) + 0.5 * WH * cos(2 * angle_increment) - 0.25 * WL * sin(2 * angle_increment),
         WD *sin(2 * angle_increment) + 0.5 * WH * sin(2 * angle_increment) + 0.25 * WL * cos(2 * angle_increment)),
    Vecd(WD * cos(3 * angle_increment) + 0.5 * WH * cos(3 * angle_increment) - 0.25 * WL * sin(3 * angle_increment),
         WD *sin(3 * angle_increment) + 0.5 * WH * sin(3 * angle_increment) + 0.25 * WL * cos(3 * angle_increment)),
    Vecd(WD * cos(4 * angle_increment) + 0.5 * WH * cos(4 * angle_increment) - 0.25 * WL * sin(4 * angle_increment),
         WD *sin(4 * angle_increment) + 0.5 * WH * sin(4 * angle_increment) + 0.25 * WL * cos(4 * angle_increment)),
    Vecd(WD * cos(5 * angle_increment) + 0.5 * WH * cos(5 * angle_increment) - 0.25 * WL * sin(5 * angle_increment),
         WD *sin(5 * angle_increment) + 0.5 * WH * sin(5 * angle_increment) + 0.25 * WL * cos(5 * angle_increment)),
    Vecd(WD * cos(6 * angle_increment) + 0.5 * WH * cos(6 * angle_increment) - 0.25 * WL * sin(6 * angle_increment),
         WD *sin(6 * angle_increment) + 0.5 * WH * sin(6 * angle_increment) + 0.25 * WL * cos(6 * angle_increment)),
    Vecd(WD * cos(7 * angle_increment) + 0.5 * WH * cos(7 * angle_increment) - 0.25 * WL * sin(7 * angle_increment),
         WD *sin(7 * angle_increment) + 0.5 * WH * sin(7 * angle_increment) + 0.25 * WL * cos(7 * angle_increment)),
    Vecd(WD * cos(8 * angle_increment) + 0.5 * WH * cos(8 * angle_increment) - 0.25 * WL * sin(8 * angle_increment),
         WD *sin(8 * angle_increment) + 0.5 * WH * sin(8 * angle_increment) + 0.25 * WL * cos(8 * angle_increment)),
    Vecd(WD * cos(9 * angle_increment) + 0.5 * WH * cos(9 * angle_increment) - 0.25 * WL * sin(9 * angle_increment),
         WD *sin(9 * angle_increment) + 0.5 * WH * sin(9 * angle_increment) + 0.25 * WL * cos(9 * angle_increment)),
    Vecd(WD * cos(10 * angle_increment) + 0.5 * WH * cos(10 * angle_increment) - 0.25 * WL * sin(10 * angle_increment),
         WD *sin(10 * angle_increment) + 0.5 * WH * sin(10 * angle_increment) + 0.25 * WL * cos(10 * angle_increment)),
    Vecd(WD * cos(11 * angle_increment) + 0.5 * WH * cos(11 * angle_increment) - 0.25 * WL * sin(11 * angle_increment),
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
//	Case-dependent wall boundary
//----------------------------------------------------------------------
Vecd center(0.0, 0.0);
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addCircle(center, (RM + BW), resolution_circle, GeometricOps::add); /**< Outer wall of motor hausing. */
        multi_polygon_.addCircle(center, RM, resolution_circle, GeometricOps::sub);        /**< Inner wall of motor hausing. */
        multi_polygon_.addBox(inlet1_transform, inlet_halfsize, GeometricOps::sub);        /**< Top Inlets. */
        multi_polygon_.addBox(outlet_transform, outlet_halfsize, GeometricOps::sub);       /**< Outlets. */
        multi_polygon_.addBox(inlet2_transform, inlet_halfsize, GeometricOps::sub);
        multi_polygon_.addBox(inlet3_transform, inlet_halfsize, GeometricOps::sub);
        multi_polygon_.addBox(inlet4_transform, inlet_halfsize, GeometricOps::sub);
        multi_polygon_.addBox(inlet5_transform, inlet_halfsize, GeometricOps::sub);
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
        multi_polygon_.addCircle(center, RR, resolution_circle, GeometricOps::add); /**< Rotor */
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
            multi_polygon_.addBox(winding_transform, winding_halfsize, GeometricOps::add);
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
        multi_polygon_.addBox(inlet1_transform, inlet_halfsize, GeometricOps::add); /**< Top Inlets. */
        multi_polygon_.addBox(inlet2_transform, inlet_halfsize, GeometricOps::add);
        multi_polygon_.addBox(inlet3_transform, inlet_halfsize, GeometricOps::add);
        multi_polygon_.addBox(inlet4_transform, inlet_halfsize, GeometricOps::add);
        multi_polygon_.addBox(inlet5_transform, inlet_halfsize, GeometricOps::add);
    }
};
//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
  public:
    InletInflowCondition(OrientedBoxByParticle &oriented_box_part)
        : EmitterInflowCondition(oriented_box_part) {}

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
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.setRunParticleRelaxation(true);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody oil_body(sph_system, makeShared<FluidBoundary>("OilBody"));
    oil_body.getSPHAdaptation().resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    oil_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    oil_body.addMaterialProperty<Viscosity>(mu_f);
    oil_body.addMaterialProperty<IsotropicDiffusion>(temperature_species_name, k_oil);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(3500.0);
    oil_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    OrientedBoxByParticle emitter1(oil_body, OrientedBox(xAxis, inlet1_transform, inlet_halfsize));
    OrientedBoxByParticle emitter2(oil_body, OrientedBox(xAxis, inlet2_transform, inlet_halfsize));
    OrientedBoxByParticle emitter3(oil_body, OrientedBox(xAxis, inlet3_transform, inlet_halfsize));
    OrientedBoxByParticle emitter4(oil_body, OrientedBox(xAxis, inlet4_transform, inlet_halfsize));
    OrientedBoxByParticle emitter5(oil_body, OrientedBox(xAxis, inlet5_transform, inlet_halfsize));
    StdVec<OrientedBoxByParticle *> emitters = {&emitter1, &emitter2, &emitter3, &emitter4, &emitter5};
    OrientedBoxByCell outlet_disposer(oil_body, OrientedBox(xAxis, outlet_transform, outlet_halfsize + 0.5 * Vec2d(BW, BW)));
    outlet_disposer.writeOrientedBoxToVtp();

    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall"));
    LevelSetShape &wall_level_set_shape = wall.defineBodyLevelSetShape().writeLevelSet();
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    RotorBoundary rotor_boundary("Rotor");
    SolidBody rotor(sph_system, rotor_boundary);
    LevelSetShape &rotor_level_set_shape = rotor.defineBodyLevelSetShape().writeLevelSet();
    rotor.defineMatterMaterial<Solid>();
    rotor.generateParticles<BaseParticles, Lattice>();

    SolidBody winding(sph_system, makeShared<WindingBoundary>("Winding"));
    LevelSetShape &winding_level_set_shape = winding.defineBodyLevelSetShape().writeLevelSet();
    winding.defineMatterMaterial<Solid>();
    winding.addMaterialProperty<IsotropicDiffusion>(temperature_species_name, k_winding);
    winding.generateParticles<BaseParticles, Lattice>();

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
        //  Setup a sub-system for particle relaxation and delete it after particle relaxation.
        RelaxationSystem relaxation_system(system_domain_bounds, global_resolution);
        //	Near surface shape for bounding particles.
        auto &near_wall_surface = wall.addBodyPart<NearShapeSurface>();
        auto &near_rotor_surface = rotor.addBodyPart<NearShapeSurface>();
        auto &near_winding_surface = winding.addBodyPart<NearShapeSurface>();
        StdVec<NearShapeSurface *> near_body_surfaces = {&near_wall_surface, &near_rotor_surface, &near_winding_surface};
        StdVec<RealBody *> real_bodies = {&wall, &rotor, &winding};
        //	Define body relation map.
        auto &wall_inner = relaxation_system.addInnerRelation(wall);
        auto &rotor_inner = relaxation_system.addInnerRelation(rotor);
        auto &winding_inner = relaxation_system.addInnerRelation(winding);
        //----------------------------------------------------------------------
        // Define SPH solver with particle methods and execution policies.
        // Generally, the host methods should be able to run immediately.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(relaxation_system);
        //----------------------------------------------------------------------
        // Define the numerical methods used in the simulation.
        //----------------------------------------------------------------------
        auto &host_methods = sph_solver.getHostMethodContainer();
        host_methods.addStateDynamics<RandomizeParticlePositionCK>(real_bodies).exec(); // host method able to run immediately
        //	Define simple file input and outputs functions.
        auto &main_methods = sph_solver.getMainMethodContainer();
        ParticleDynamicsGroup update_cell_linked_list;
        update_cell_linked_list.add(&main_methods.addCellLinkedListDynamics(wall));
        update_cell_linked_list.add(&main_methods.addCellLinkedListDynamics(rotor));
        update_cell_linked_list.add(&main_methods.addCellLinkedListDynamics(winding));
        ParticleDynamicsGroup update_relation;
        update_relation.add(&main_methods.addRelationDynamics(wall_inner));
        update_relation.add(&main_methods.addRelationDynamics(rotor_inner));
        update_relation.add(&main_methods.addRelationDynamics(winding_inner));
        ParticleDynamicsGroup update_configuration = update_cell_linked_list + update_relation;

        ParticleDynamicsGroup relaxation_residual;
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(wall_inner)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(wall, wall_level_set_shape));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(rotor_inner)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(rotor, rotor_level_set_shape));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(winding_inner)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(winding, winding_level_set_shape));

        ReduceDynamicsGroup relaxation_scaling = main_methods.addReduceDynamics<ReduceMin<Real>, RelaxationScalingCK>(real_bodies);
        ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(real_bodies);
        ParticleDynamicsGroup level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surfaces);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(relaxation_system);
        //----------------------------------------------------------------------
        //	Prepare for the time integration loop.
        //----------------------------------------------------------------------
        level_set_bounding.exec();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        body_state_recorder.writeToFile(0);
        while (ite_p < 2000)
        {
            update_configuration.exec();
            relaxation_residual.exec();
            Real relaxation_step = relaxation_scaling.exec();
            update_particle_position.exec(relaxation_step);
            level_set_bounding.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }

        host_methods.addStateDynamics<NormalFromBodyShapeCK>(real_bodies).exec(); // host method run after particle relaxation
        std::cout << "The physics relaxation process finish !" << std::endl;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &oil_body_inner = sph_system.addInnerRelation(oil_body);
    auto &winding_inner = sph_system.addInnerRelation(winding);

    auto &oil_wall_contact = sph_system.addContactRelation(oil_body, wall);
    auto &oil_rotor_contact = sph_system.addContactRelation(oil_body, rotor);
    auto &oil_winding_contact = sph_system.addContactRelation(oil_body, winding);

    auto &winding_oil_contact = sph_system.addContactRelation(winding, oil_body);

    auto &slot_observer_contact = sph_system.addContactRelation(slot_observer, winding);
    auto &tooth_observer_contact = sph_system.addContactRelation(tooth_observer, winding);
    auto &yoke_observer_contact = sph_system.addContactRelation(yoke_observer, winding);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    ParticleDynamicsGroup static_configuration;
    static_configuration.add(&main_methods.addCellLinkedListDynamics(wall));
    static_configuration.add(&main_methods.addCellLinkedListDynamics(winding));
    static_configuration.add(&main_methods.addRelationDynamics(winding_inner));

    ParticleDynamicsGroup update_configuration;
    update_configuration.add(&main_methods.addCellLinkedListDynamics(oil_body));
    update_configuration.add(&main_methods.addCellLinkedListDynamics(rotor));
    update_configuration.add(&main_methods.addRelationDynamics(
        oil_body_inner, oil_wall_contact, oil_rotor_contact, oil_winding_contact));
    update_configuration.add(&main_methods.addRelationDynamics(winding_oil_contact));
    auto &static_observer_relation = main_methods.addRelationDynamics(
        slot_observer_contact, tooth_observer_contact, yoke_observer_contact);
    auto &particle_sort = main_methods.addSortDynamics(oil_body);

    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(oil_body, gravity);
    auto &advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(oil_body);
    auto &update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(oil_body);

    auto &linear_correction_matrix =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(oil_body_inner, 0.5)
            .addPostContactInteraction(oil_wall_contact)
            .addPostContactInteraction(oil_rotor_contact)
            .addPostContactInteraction(oil_winding_contact);

    auto &acoustic_step_1st_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_body_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_wall_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_rotor_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_winding_contact);
    auto &acoustic_step_2nd_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_body_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_wall_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_rotor_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(oil_winding_contact);
    auto &density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(oil_body_inner)
            .addPostContactInteraction(oil_wall_contact)
            .addPostContactInteraction(oil_rotor_contact)
            .addPostContactInteraction(oil_winding_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeSurface>(oil_body);

    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(oil_body, U_f);
    auto &acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(oil_body);

    auto &inflow_condition = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowConditionCK, ConstantInflowSpeed>(emitters, v_inlet);
    auto &emitter_injection = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK>(emitters);
    auto &within_deposer_indication = main_methods.addStateDynamics<fluid_dynamics::WithinDisposerIndication>(outlet_disposer);
    auto &outflow_particle_deletion = main_methods.addStateDynamics<fluid_dynamics::OutflowParticleDeletion>(oil_body);

    IsotropicDiffusion oil_isotropic_diffusion(temperature_species_name, k_oil);
    IsotropicDiffusion winding_isotropic_diffusion(temperature_species_name, k_winding);
    auto &thermal_time_step = main_methods.addReduceDynamicsGroup<ReduceMin<Real>>();
    thermal_time_step.add(&main_methods.addReturnDynamics<GetDiffusionTimeStepSize>(oil_body, &oil_isotropic_diffusion));
    thermal_time_step.add(&main_methods.addReturnDynamics<GetDiffusionTimeStepSize>(winding, &winding_isotropic_diffusion));

    auto &all_runge_kutta_1st_stage = main_methods.addParticleDynamicsGroup();
    auto &all_runge_kutta_2nd_stage = main_methods.addParticleDynamicsGroup();

    auto &oil_runge_kutta_1st_stage =
        main_methods.addInteractionDynamicsOneLevel<
            DiffusionRelaxationCK, RungeKutta1stStage, IsotropicDiffusion, LinearCorrectionCK>(
            oil_body_inner, &oil_isotropic_diffusion);

    auto &winding_runge_kutta_1st_stage =
        main_methods.addInteractionDynamicsOneLevel<
            DiffusionRelaxationCK, RungeKutta1stStage, IsotropicDiffusion, LinearCorrectionCK>(
            winding_inner, &winding_isotropic_diffusion);

    //----------------------------------------------------------------------
    //	Define the multi-body system
    //----------------------------------------------------------------------
    SimbodySystem simbody_system;
    SolidBodyPartForSimbodyCK rotor_multibody(rotor, rotor_boundary);
    std::string rotor_name = simbody_system.createRigidBody(rotor_multibody);
    SimTK::MobilizedBody &rotor_mob = simbody_system.createMobilizedBody(rotor_name, "Ground", "Pin");
    simbody_system.setUForMobilizedPin(rotor_name, Omega);
    simbody_system.initializeStateForIntegrator();
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    auto &constraint_on_rotor = main_methods.addStateDynamics<
        solid_dynamics::ConstraintBodyPartBySimBodyCK>(rotor_multibody, simbody_system);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    size_t advection_steps = 0;
    int screening_interval = 100;
    Real end_time = 2.0;
    auto &state_recording = time_stepper.addTriggerByInterval(0.05);
    //----------------------------------------------------------------------
    //	Construct static cell linked list and Lagrangian configurations.
    //----------------------------------------------------------------------
    static_configuration.exec();
    update_configuration.exec();
    static_observer_relation.exec();
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    constant_gravity.exec();
    constraint_on_rotor.exec();
    inflow_condition.exec();

    density_regularization.exec();
    advection_step_setup.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    /** statistics for computing time. */
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime(end_time))
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        inflow_condition.exec();
        acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(advection_time_step))
        {
            advection_steps++;
            update_particle_position.exec();
            // simbody dynamics and constraint
            Real t_target = time_stepper.getPhysicalTime();
            if (t_target > simbody_system.getSimbodySystemTime())
            {
                simbody_system.stepSimbodySystemTo(t_target);
            }
            constraint_on_rotor.exec();

            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Physical Time = " << time_stepper.getPhysicalTime()
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (state_recording())
            {
                body_states_recording.writeToFile();
            }

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            emitter_injection.exec();
            within_deposer_indication.exec();
            outflow_particle_deletion.exec();
            if (advection_steps % 100)
            {
                particle_sort.exec();
            }
            update_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            density_regularization.exec();
            advection_step_setup.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for real computations.
    //----------------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}
