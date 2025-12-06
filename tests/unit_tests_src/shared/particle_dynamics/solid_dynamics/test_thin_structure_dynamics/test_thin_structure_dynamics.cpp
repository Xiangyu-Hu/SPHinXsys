/**
 * @file 	test_thin_structure_dynamics.cpp
 * @brief 	Shell model validation
 * @details A rigid rotation of a 3D plate is simulated and tested by GTest to validate the shell model.
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 10.0;                                   /** Length of the square plate. */
Real PH = 10.0;                                   /** Width of the square plate. */
Real PT = 1.0;                                    /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                 /** Pseudo-normal. */
int particle_number = 40;                         /** Particle number in the direction of the length */
Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;                                      /** Width of the boundary layer measured by number of particles. */
Real BW = resolution_ref * (Real)BWD;             /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
                                 Vec3d(PL + BW, PH + BW, 0.5 * resolution_ref));

/** For material properties of the solid. */
Real rho0_s = 1.0;                 /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;                /** Poisson ratio. */

StdVec<int> random_index;
Vecd *pseudo_normal = nullptr;
Vecd *normal = nullptr;
StdVec<Real> von_mises_strain;
TEST(Plate, RigidRotationTest)
{
    for (size_t i = 0; i < random_index.size(); i++)
    {
        EXPECT_NEAR(pseudo_normal[random_index[i]][0], normal[random_index[i]][0], 1.0e-3);
        EXPECT_NEAR(pseudo_normal[random_index[i]][1], normal[random_index[i]][1], 1.0e-3);
        EXPECT_NEAR(pseudo_normal[random_index[i]][2], normal[random_index[i]][2], 1.0e-3);
        EXPECT_NEAR(von_mises_strain[i], 0.0, 1.0e-3);

        std::cout << "pseudo_normal[" << random_index[i] << "]: " << pseudo_normal[random_index[i]] << std::endl;
        std::cout << "normal[" << random_index[i] << "]: " << normal[random_index[i]] << std::endl;
        std::cout << "von_mises_strain[" << random_index[i] << "]: " << von_mises_strain[i] << std::endl;
    }
}

namespace SPH
{
/** Define application dependent particle generator for thin structure. */
class Plate;
template <>
class ParticleGenerator<SurfaceParticles, Plate> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    virtual void prepareGeometricData() override
    {
        // the plate and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {
            for (int j = 0; j < (particle_number + 2 * BWD); j++)
            {
                Real x = resolution_ref * i - BW + resolution_ref * 0.5 - PL * 0.5;
                Real y = resolution_ref * j - BW + resolution_ref * 0.5 - PH * 0.5;
                addPositionAndVolumetricMeasure(Vecd(x, y, 0.0), resolution_ref * resolution_ref);
                addSurfaceProperties(n_0, PT);
            }
        }
    }
};
/** Define the controlled geometry. */
class ControlledGeometry : public BodyPartByParticle
{
  public:
    ControlledGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ControlledGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~ControlledGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        return true;
    };
};
/** Define the controlled rotation. */
class ControlledRotation : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    ControlledRotation(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          angular_vel_(particles_->getVariableDataByName<Vecd>("AngularVelocity")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {};
    virtual ~ControlledRotation() {};

  protected:
    Vecd *vel_, *angular_vel_, *pos_;
    Real *physical_time_;
    Real rotation_v = Pi;
    void update(size_t index_i, Real dt = 0.0)
    {
        Real current_time = *physical_time_;
        if (current_time <= 0.5)
        {
            vel_[index_i] = Vecd(0.0, -rotation_v * pos_[index_i][2], rotation_v * pos_[index_i][1]);
            angular_vel_[index_i] = Vecd(rotation_v, 0.0, 0.0);
        }
        else
        {
            vel_[index_i] = Vecd(-rotation_v * pos_[index_i][1], rotation_v * pos_[index_i][0], 0.0);
            angular_vel_[index_i] = Vecd(0.0, rotation_v, 0.0);
        }
    };
};
} // namespace SPH
/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem system(system_domain_bounds, resolution_ref);

    /** create a plate body. */
    SolidBody plate_body(system, makeShared<DefaultShape>("PlateBody"));
    plate_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    plate_body.generateParticles<SurfaceParticles, Plate>();
    auto shell_particles = dynamic_cast<SurfaceParticles *>(&plate_body.getBaseParticles());

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation plate_body_inner(plate_body);

    /**
     * This section define all numerical methods will be used in this case.
     */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(plate_body_inner);

    /** active-passive stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(plate_body_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(plate_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
    /** Constrain the Boundary. */
    ControlledGeometry controlled_geometry(plate_body);
    SimpleDynamics<ControlledRotation> controlled_rotation(controlled_geometry);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal(plate_body);
    /** File and screen outputs */
    BodyStatesRecordingToVtp write_states(system);
    write_states.addToWrite<Vecd>(plate_body, "PseudoNormal");
    write_states.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(plate_body);
    Real *all_von_mises_strain = shell_particles->getVariableDataByName<Real>("VonMisesStrain");

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();
    controlled_rotation.exec();

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    write_states.writeToFile(0);

    /** Setup physical parameters. */
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 1.1;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (physical_time < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            if (physical_time > 0.5 && physical_time < 0.51)
            {
                controlled_rotation.exec();
            }
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            physical_time += dt;
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    for (int i = 0; i < 10; i++)
    {
        random_index.push_back(rand_uniform(0.0, 1.0) * shell_particles->TotalRealParticles());
        von_mises_strain.push_back(all_von_mises_strain[random_index[i]]);
    }

    update_normal.exec();

    pseudo_normal = shell_particles->getVariableDataByName<Vecd>("PseudoNormal");
    normal = shell_particles->getVariableDataByName<Vecd>("NormalDirection");

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
