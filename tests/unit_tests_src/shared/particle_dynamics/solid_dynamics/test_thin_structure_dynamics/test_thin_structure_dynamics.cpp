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
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
                                 Vec3d(PL + BW, PH + BW, 0.5 * resolution_ref));

/** For material properties of the solid. */
Real rho0_s = 1.0;                 /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;                /** Poisson ratio. */

StdVec<int> rondom_index;
StdLargeVec<Vecd> pseudo_normal;
StdLargeVec<Vecd> normal;
StdVec<Real> von_mises_strain;
TEST(Plate, RigidRotationTest)
{
    for (size_t i = 0; i < rondom_index.size(); i++)
    {
        EXPECT_NEAR(pseudo_normal[rondom_index[i]][0], normal[rondom_index[i]][0], 1.0e-3);
        EXPECT_NEAR(pseudo_normal[rondom_index[i]][1], normal[rondom_index[i]][1], 1.0e-3);
        EXPECT_NEAR(pseudo_normal[rondom_index[i]][2], normal[rondom_index[i]][2], 1.0e-3);
        EXPECT_NEAR(von_mises_strain[i], 0.0, 1.0e-3);

        std::cout << "pseudo_normal[" << rondom_index[i] << "]: " << pseudo_normal[rondom_index[i]] << std::endl;
        std::cout << "normal[" << rondom_index[i] << "]: " << normal[rondom_index[i]] << std::endl;
        std::cout << "von_mises_strain[" << rondom_index[i] << "]: " << von_mises_strain[i] << std::endl;
    }
}

/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit PlateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        // the plate and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {
            for (int j = 0; j < (particle_number + 2 * BWD); j++)
            {
                Real x = resolution_ref * i - BW + resolution_ref * 0.5 - PL * 0.5;
                Real y = resolution_ref * j - BW + resolution_ref * 0.5 - PH * 0.5;
                initializePositionAndVolumetricMeasure(Vecd(x, y, 0.0), resolution_ref * resolution_ref);
                initializeSurfaceProperties(n_0, PT);
            }
        }
    }
};
/** Define the controled geometry. */
class ControledGeometry : public BodyPartByParticle
{
  public:
    ControledGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ControledGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~ControledGeometry(){};

  private:
    void tagManually(size_t index_i)
    {
        body_part_particles_.push_back(index_i);
    };
};
/** Define the controled rotation. */
class ControledRotation : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    ControledRotation(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part),
          vel_(particles_->vel_), angular_vel_(particles_->angular_vel_), pos_(particles_->pos_){};
    virtual ~ControledRotation(){};

  protected:
    StdLargeVec<Vecd> &vel_, &angular_vel_, &pos_;
    Real ratation_v = Pi;
    void update(size_t index_i, Real dt = 0.0)
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        if (current_time <= 0.5)
        {
            vel_[index_i] = Vecd(0.0, -ratation_v * pos_[index_i][2], ratation_v * pos_[index_i][1]);
            angular_vel_[index_i] = Vecd(ratation_v, 0.0, 0.0);
        }
        else
        {
            vel_[index_i] = Vecd(-ratation_v * pos_[index_i][1], ratation_v * pos_[index_i][0], 0.0);
            angular_vel_[index_i] = Vecd(0.0, ratation_v, 0.0);
        }
    };
};
/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem system(system_domain_bounds, resolution_ref);

    /** create a plate body. */
    SolidBody plate_body(system, makeShared<DefaultShape>("PlateBody"));
    plate_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    plate_body.generateParticles<PlateParticleGenerator>();
    auto shell_particles = dynamic_cast<ShellParticles *>(&plate_body.getBaseParticles());
    plate_body.addBodyStateForRecording<Vecd>("PseudoNormal");

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation plate_body_inner(plate_body);

    /**
     * This section define all numerical methods will be used in this case.
     */
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
        corrected_configuration(plate_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
    /** active-passive stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
        stress_relaxation_first_half(plate_body_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
        stress_relaxation_second_half(plate_body_inner);
    /** Constrain the Boundary. */
    ControledGeometry controled_geometry(plate_body, "ControledGeometry");
    SimpleDynamics<ControledRotation> controled_rotaton(controled_geometry);
    /** Output */
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();
    controled_rotaton.exec();

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    write_states.writeToFile(0);

    /** Setup physical parameters. */
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            if (GlobalStaticVariables::physical_time_ > 0.5 && GlobalStaticVariables::physical_time_ < 0.51)
            {
                controled_rotaton.exec();
            }
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
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
        rondom_index.push_back((double)rand() / (RAND_MAX)*shell_particles->total_real_particles_);
        von_mises_strain.push_back(shell_particles->getVonMisesStrain(rondom_index[i]));
    }
    pseudo_normal = shell_particles->pseudo_n_;
    normal = shell_particles->n_;

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
