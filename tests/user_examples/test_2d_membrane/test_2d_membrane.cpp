/**
 * @file 	membrane.cpp
 * @brief 	fluid diffusion coupled with porous structure deformation  
 * @details This is the one of the test cases using fluid diffusion inside the porous elastic media model.                              *
             In this case, the multi-time step algorithm is used.   
 * @author 	Xiaojing Tang and Xiangyu Hu
 */                                           
#include "particle_momentum_dissipation.hpp"
#include "porous_media_dynamics.h"
#include "porous_media_solid.h"
#include "porous_solid_particles.h"
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 10.0;  // membrane length
Real PH = 0.125; // membrane thickness
Real BC = PL * 0.15; // half length of the saturation

int y_num = 8;
// reference particle spacing
Real resolution_ref = PH / y_num;

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-PL, -PL),
                                 Vec2d(2.0 * PL, PL));

//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho_0 = 2.0;  // reference solid density non-dimensional 
Real poisson = 0.26316;  /**< Poisson ratio. */
Real Youngs_modulus = 8.242e6;  /**< Youngs modulus. */
Real physical_viscosity = 5000.0;
Real saturation = 0.4;   /**< solid porosity. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real diffusivity_constant_ = 1.0e-4; // reference diffusion constant non-dimensional 
Real fluid_initial_density_ = 1.0; // reference fluid density non-dimensional
Real water_pressure_constant_ = 3.0e6;

// the criterion to represent the static state 
Real refer_density_energy = 0.5 * water_pressure_constant_;

//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a membrane base shape
std::vector<Vecd> membrane_base_shape{
    Vecd(-resolution_ref * 3.0, -PH / 2.0), Vecd(-resolution_ref * 3.0, PH / 2.0), Vecd(0.0, PH / 2.0),
    Vecd(0.0, -PH / 2.0), Vecd(-resolution_ref * 3.0, -PH / 2.0)};

// a membrane shape
std::vector<Vecd> membrane_shape{Vecd(0.0, -PH / 2.0), Vecd(0.0, PH / 2.0),
                             Vecd(PL, PH / 2.0), Vecd(PL, -PH / 2.0), Vecd(0.0, -PH / 2.0)};

// a membrane end shape
std::vector<Vecd> membrane_end_shape{
    Vecd(PL, -PH / 2.0), Vecd(PL, PH / 2.0),
    Vecd(PL + 4.0 * resolution_ref, PH / 2.0), Vecd(PL + 4.0 * resolution_ref, -PH / 2.0),
    Vecd(PL, -PH / 2.0)};

// a membrane saturation shape
std::vector<Vecd> membrane_saturation_shape{
    Vecd(PL / 2.0 - BC, 0.0), Vecd(PL / 2.0 - BC, PH / 2.0), Vecd(PL / 2.0 + BC, PH / 2.0),
    Vecd(PL / 2.0 + BC, 0.0), Vecd(PL / 2.0 - BC, 0.0)};

// membrane observer location
StdVec<Vecd> observation_location = {Vecd(PL / 2.0, 0.0)};

//----------------------------------------------------------------------
//	Define the membrane body
//----------------------------------------------------------------------
class Membrane : public MultiPolygonShape
{
  public:
    explicit Membrane(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(membrane_base_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(membrane_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(membrane_end_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	define the membrane base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createMembraneConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(membrane_base_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(membrane_shape, ShapeBooleanOps::sub);
    multi_polygon.addAPolygon(membrane_end_shape, ShapeBooleanOps::add);
    return multi_polygon;
};

//----------------------------------------------------------------------
//	define the membrane part which the saturation condition is applied.
//----------------------------------------------------------------------
MultiPolygon createSaturationConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(membrane_saturation_shape, ShapeBooleanOps::add);
    return multi_polygon;
};

//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class SaturationInitialCondition : public multi_species_continuum::PorousMediaSaturationDynamicsInitialCondition
{
  public:
    SaturationInitialCondition(BodyPartByParticle &body_part) : multi_species_continuum::PorousMediaSaturationDynamicsInitialCondition(body_part){};
    virtual ~SaturationInitialCondition(){};

  protected:
    void update(size_t index_i, Real dt = 0.0)
    {
        fluid_saturation_[index_i] = saturation;
        fluid_mass_[index_i] = saturation * fluid_initial_density_ * Vol_update_[index_i];
        total_mass_[index_i] = rho_n_[index_i] * Vol_update_[index_i] + fluid_mass_[index_i];
    };
};

//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
       //	Creating body, materials and particles.
       //----------------------------------------------------------------------
    SolidBody membrane(sph_system, makeShared<Membrane>("MembraneBody"));
    membrane.defineParticlesAndMaterial<multi_species_continuum::PorousMediaParticles, multi_species_continuum::PorousMediaSolid>(
        rho_0, Youngs_modulus, poisson, diffusivity_constant_, fluid_initial_density_, water_pressure_constant_);
    membrane.generateParticles<ParticleGeneratorLattice>();

    ObserverBody membrane_observer(sph_system, "MembraneObserver");
    membrane_observer.defineAdaptationRatios(1.15, 2.0);
    membrane_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation membrane_body_inner(membrane);
    ContactRelation membrane_observer_contact(membrane_observer, {&membrane});
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    // corrected strong configuration
    InteractionWithUpdate<CorrectedConfigurationInner> membrane_corrected_configuration(membrane_body_inner);
    // time step size calculation for solid stress relaxation
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(membrane);
    // time step size calculation for fluid diffusion relaxation
    ReduceDynamics<multi_species_continuum::GetSaturationTimeStepSize> saturation_time_step_size(membrane);

    // stress relaxation for the membrane
    Dynamics1Level<multi_species_continuum::PorousMediaStressRelaxationFirstHalf> stress_relaxation_first_half(membrane_body_inner);
    Dynamics1Level<multi_species_continuum::PorousMediaStressRelaxationSecondHalf> stress_relaxation_second_half(membrane_body_inner);
    // fluid diffusion relaxation inside the membrane
     Dynamics1Level<multi_species_continuum::SaturationRelaxationInPorousMedia> saturation_relaxation(membrane_body_inner);

    // clamping a solid body part using momentum constraint  
    BodyRegionByParticle membrane_base(membrane, makeShared<MultiPolygonShape>(createMembraneConstrainShape()));
    SimpleDynamics<multi_species_continuum::MomentumConstraint> clamp_constrain_membrane_base(membrane_base);

    // applying the saturation boundary condition. 
    BodyRegionByParticle membrane_saturation(membrane, makeShared<MultiPolygonShape>(createSaturationConstrainShape()));
    SimpleDynamics<SaturationInitialCondition> constrain_membrane_saturation(membrane_saturation);

    //total mechanical energy to check the static state is achieved or not
    ReduceDynamics<TotalMechanicalEnergy> get_kinetic_energy(membrane);

    /** Damping applied on momentum */
    DampingWithRandomChoice<InteractionSplit<multi_species_continuum::PorousMediaDampingPairwiseInner<Vec2d>>>
        membrane_damping(0.5, membrane_body_inner, "TotalMomentum", physical_viscosity);
    
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp write_membrane_states(io_environment, sph_system.real_bodies_);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_membrane_tip_position("Position", io_environment, membrane_observer_contact);
    
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    constrain_membrane_saturation.exec();
    membrane_corrected_configuration.exec();

	//----------------------------------------------------------------------
	//	Setup computing time-step controls.
	//----------------------------------------------------------------------
    int ite = 0;
    int total_ite = 0;
    GlobalStaticVariables::physical_time_ = 0.0;

    Real End_Time = 100.0;
    // time for applying saturation condition
    Real setup_saturation_time_ = End_Time * 0.1;
    // time step size for output file
    Real D_Time = End_Time / 100.0;
    Real dt = 0.0; // default acoustic time step sizes

    // statistics for computing time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    write_membrane_states.writeToFile(0);
    write_membrane_tip_position.writeToFile(0);
    // computation loop starts
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;
        // the outer loop for fluid diffusion, also integrate time (loop) until the next output time
        while (integration_time < D_Time)
        {
            Real Dt = saturation_time_step_size.exec();
            if (GlobalStaticVariables::physical_time_ < setup_saturation_time_)
            {
                constrain_membrane_saturation.exec();
            }
            saturation_relaxation.exec(Dt);

            int stress_ite = 0; // solid stress relaxation times
            Real relaxation_time = 0.0;
            Real total_kinetic_energy = 1000.0;

            while (relaxation_time < Dt)
            {
                //the inner loop for solid stress relaxation
                if (total_kinetic_energy > (5e-9 * refer_density_energy)) 
                {
                    stress_relaxation_first_half.exec(dt);
                    clamp_constrain_membrane_base.exec();
                    membrane_damping.exec(dt);
                    clamp_constrain_membrane_base.exec();
                    stress_relaxation_second_half.exec(dt);

                    //calculate the total kinetic energy to check the static state
                    total_kinetic_energy = get_kinetic_energy.exec();
                    ite++;
                    stress_ite++;
                    dt = SMIN(computing_time_step_size.exec(), Dt);

                    if (ite % 10000 == 0)
                    {
                        std::cout << "N=" << ite << " Time: "
                                  << GlobalStaticVariables::physical_time_ << "  Dt:" << Dt << "	dt: "
                                  << dt << "  Dt/ dt:" << Dt / dt << "\n";
                    }
                }

                total_ite++;
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            write_membrane_tip_position.writeToFile(ite);
            std::cout << "One Diffusion finishes   "
                       << "total_kinetic_energy =  " << total_kinetic_energy
                       << "     stress_ite = " << stress_ite << std::endl;   
        }

        TickCount t2 = TickCount::now();
        write_membrane_states.writeToFile(ite);      

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds."
              << "  Iterations:  " << ite << std::endl;
    std::cout << "Total iterations computation:  " << GlobalStaticVariables::physical_time_ / dt
              << "  Total iterations:  " << total_ite << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_membrane_tip_position.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
    }
    else
    {
        write_membrane_tip_position.testResult();
    }
    return 0;
}
