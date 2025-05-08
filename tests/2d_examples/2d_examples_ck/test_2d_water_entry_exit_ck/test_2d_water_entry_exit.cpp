/**
 * @file	water entry and exit.cpp
 * @brief	2D water entry and exit example with surface wetting considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach,
 *          especially when coupled with the wetting.
 * @author  Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys_ck.h" //SPHinXsys Library.
using namespace SPH;      // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real cylinder_radius = 0.055;                             /**< Cylinder radius. */
Real DL = 8.0 * cylinder_radius;                          /**< Water tank length. */
Real DH = 7.0 * cylinder_radius;                          /**< Water tank height. */
Real LL = DL;                                             /**< Water column length. */
Real LH = 3.0 * cylinder_radius;                          /**< Water column height. */
Real particle_spacing_ref = 2.0 * cylinder_radius / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                       /**< Thickness of tank wall. */
Vec2d cylinder_center(0.5 * DL, LH + 0.15);               /**< Location of the cylinder center. */
Vecd tethering_point(0.5 * DL, DH);                       /**< The tethering point. */
StdVec<Vecd> observer_location = {cylinder_center};       /**< Displacement observation point. */
StdVec<Vecd> wetting_observer_location =
    {cylinder_center - Vecd(0.0, cylinder_radius)}; /**< wetting observation point. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Fluid density. */
Real rho0_s = 0.5;                       /**< Cylinder density. */
Real gravity_g = 9.81;                   /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
Real mu_f = 8.9e-7;                      /**< Water dynamics viscosity. */
//----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 100.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real cylinder_moisture = 0.0;                                /**< cylinder moisture. */
Real wall_moisture = 1.0;                                    /**< wall moisture. */
//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block;
    water_block.push_back(Vecd(0.0, 0.0));
    water_block.push_back(Vecd(0.0, LH));
    water_block.push_back(Vecd(LL, LH));
    water_block.push_back(Vecd(LL, 0.0));
    water_block.push_back(Vecd(0.0, 0.0));

    return water_block;
}
class WettingFluidBody : public MultiPolygonShape
{
  public:
    explicit WettingFluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

class WettingFluidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingFluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = fluid_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall;
    outer_wall.push_back(Vecd(-BW, -BW));
    outer_wall.push_back(Vecd(-BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, -BW));
    outer_wall.push_back(Vecd(-BW, -BW));

    return outer_wall;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall;
    inner_wall.push_back(Vecd(0.0, 0.0));
    inner_wall.push_back(Vecd(0.0, DH));
    inner_wall.push_back(Vecd(DL, DH));
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));

    return inner_wall;
}
class WettingWallBody : public MultiPolygonShape
{
  public:
    explicit WettingWallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
class WettingWallBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingWallBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = wall_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for cylinder body
//----------------------------------------------------------------------
class WettingCylinderBody : public MultiPolygonShape
{
  public:
    explicit WettingCylinderBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    }
};
class WettingCylinderBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingCylinderBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = cylinder_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//------------------------------------------------------------------------------
// Constrained part for Simbody
//------------------------------------------------------------------------------
MultiPolygon createSimbodyConstrainShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WettingFluidBody>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WettingWallBody>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder(sph_system, makeShared<WettingCylinderBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 1.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineMaterial<Solid>(rho0_s);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();

    ObserverBody cylinder_observer(sph_system, "CylinderObserver");
    cylinder_observer.generateParticles<ObserverParticles>(observer_location);

    ObserverBody wetting_observer(sph_system, "WettingObserver");
    wetting_observer.generateParticles<ObserverParticles>(wetting_observer_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(cylinder);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(cylinder);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_block_contact(water_block, {&wall_boundary, &cylinder});
    Contact<> cylinder_contact(cylinder, {&water_block});
    Contact<> cylinder_observer_contact(cylinder_observer, {&cylinder});
    Contact<> wetting_observer_contact(wetting_observer, {&cylinder});
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> cylinder_cell_linked_list(cylinder);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>>
        water_block_update_complex_relation(water_block_inner, water_block_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>>
        cylinder_update_contact_relation(cylinder_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>>
        cylinder_observer_update_contact_relation(cylinder_observer_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>>
        wetting_observer_update_contact_relation(wetting_observer_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);

    IsotropicDiffusion diffusion(diffusion_species_name, diffusion_species_name, diffusion_coeff);
    GetDiffusionTimeStepSize get_thermal_time_step(cylinder, &diffusion);
    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, UniformDistribution<Real>>>
        cylinder_initial_condition(cylinder, diffusion_species_name, wall_moisture);
    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, UniformDistribution<Real>>>
        water_block_initial_condition(water_block, diffusion_species_name, fluid_moisture);
    DynamicsSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<Contact<OneLevel, RungeKutta1stStage, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<Contact<OneLevel, RungeKutta2ndStage, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(cylinder_contact, &diffusion), DynamicsArgs(cylinder_contact, &diffusion));

    return 0;
};
