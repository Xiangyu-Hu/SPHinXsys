/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case by testing the relaxation with evolution method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;
Real DH = 1.0;
Real insert_circle_radius = 1;
Vec2d insert_circle_center(0.0, 0.0);
Real resolution_ref = 1 / 10.0;
Real BW = resolution_ref * 4;
BoundingBox system_domain_bounds(Vec2d(-2 * BW - DL, -2 * BW - DH), Vec2d(DL + 2 * BW, DH + 2 * BW));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
std::vector<Vecd> createBlockDomain()
{
    std::vector<Vecd> blockDomain;
    blockDomain.push_back(Vecd(-BW, -BW));
    blockDomain.push_back(Vecd(-BW, DH + BW));
    blockDomain.push_back(Vecd(DL + BW, DH + BW));
    blockDomain.push_back(Vecd(DL + BW, -BW));
    blockDomain.push_back(Vecd(-BW, -BW));

    return blockDomain;
};

std::vector<Vecd> createAverageDomain()
{
    std::vector<Vecd> averageDomain;
    averageDomain.push_back(Vecd(0.0, 0.0));
    averageDomain.push_back(Vecd(0.0, DH));
    averageDomain.push_back(Vecd(DL, DH));
    averageDomain.push_back(Vecd(DL, 0.0));
    averageDomain.push_back(Vecd(0.0, 0.0));

    return averageDomain;
};

class Block : public MultiPolygonShape
{
public:
    explicit Block(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createBlockDomain(), ShapeBooleanOps::add);
    }
};

MultiPolygon averageDomain()
{
    MultiPolygon multi_polygon;
    //multi_polygon.addAPolygon(createAverageDomain(), ShapeBooleanOps::add);
    multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
    return multi_polygon;
}

class Circle : public ComplexShape
{
public:
    explicit Circle(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon multi_polygon;
        multi_polygon.addACircle(insert_circle_center, insert_circle_radius + BW, 100, ShapeBooleanOps::add);
        add<MultiPolygonShape>(multi_polygon);
    }
};

class TestingInitialCondition : public fluid_dynamics::FluidInitialCondition
{
public: 
    explicit TestingInitialCondition(SPHBody& sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_)
    {
        particles_->registerVariable(scalar_, "Scalar");
        particles_->registerVariable(vector_, "Vector");
        particles_->registerVariable(matrix_, "Matrix");
    };

    void update(size_t index_i, Real dt)
    {
        /* initial pressure distribution. */
        //scalar_[index_i] = sin(2 * Pi * pos_[index_i][0]);
        scalar_[index_i] = exp(-pos_[index_i][0] * pos_[index_i][0] / 0.1);
    }

protected:
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Real> scalar_;
    StdLargeVec<Vecd> vector_;
    StdLargeVec<Matd> matrix_;
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(true);
    sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody body(sph_system, makeShared<Circle>("Body"));
    body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1, 1, 1);
    body.addBodyStateForRecording<Vecd>("Position");
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
        : body.generateParticles<ParticleGeneratorLattice>();

    BodyRegionByParticle average_domain(body, makeShared<MultiPolygonShape>(averageDomain(), "average_domain"));
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map.
        //----------------------------------------------------------------------
        InnerRelation body_inner(body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_body_to_vtp(io_environment, {&body});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&body});
        /* Relaxation method including based on the 0th and 1st order consistency. */
        InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> correction_matrix(body_inner);
        body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
        relax_dynamics::RelaxationStepInner relaxation_inner(body_inner, true);
        relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner_implicit(body_inner, true);
        /* Update the relaxation residual. */
        SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
        InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_skgc_realization(body_inner, true);
        InteractionDynamics<relax_dynamics::CheckReverseConsistencyRealization> check_rkgc_realization(body_inner, true);
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_ac(average_domain, "ACTERMNORM");
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_as(average_domain, "ASTERMNORM");
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_c(average_domain, "CTERMNORM");
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_ab(average_domain, "ABTERMNORM");
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_ar(average_domain, "ARTERMNORM");
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
            calculate_body_average_b(average_domain, "BTERMNORM");
        /* Update the kinetic energy for stopping the relaxation. */
        SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(body_inner);
        ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>> 
            calculate_body_average_kinetic_energy(average_domain, "ParticleKineticEnergy");
        ReduceDynamics<QuantityMaximum<Real>> calculate_body_maximum_kinetic_energy(body, "ParticleKineticEnergy");
        std::string filefullpath_error_analysis = io_environment.output_folder_ + "/" + "error_analysis.dat";
        std::ofstream out_file_error_analysis(filefullpath_error_analysis.c_str(), std::ios::app);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.22);
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        testing_initial_condition.exec();
        write_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        TickCount t1 = TickCount::now();
        int ite = 0;

        Real body_average_kinetic_energy = 100.0;
        Real body_maximum_kinetic_energy = 100.0;
        Real ac_term = 0;
        Real as_term = 0;
        Real c_term = 0;
        Real ab_term = 0;
        Real ar_term = 0;
        Real b_term = 0;

        /* Initial error analysis */
        correction_matrix.exec();
        check_skgc_realization.exec();
        check_rkgc_realization.exec();
        ac_term = calculate_body_average_ac.exec();
        as_term = calculate_body_average_as.exec();
        c_term = calculate_body_average_c.exec();
        ab_term = calculate_body_average_ab.exec();
        ar_term = calculate_body_average_ar.exec();
        b_term = calculate_body_average_b.exec();
        out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " <<
            ac_term << " " << as_term << " " << c_term << " " <<
            ab_term << " " << ar_term << " " << b_term << "\n";

        GlobalStaticVariables::physical_time_ = ite;
        while (ite < 10000)
        {
            correction_matrix.exec();
            relaxation_inner_implicit.exec();
            ite++;

            if (ite % 50 == 0)
            {
                testing_initial_condition.exec();
                correction_matrix.exec();
                check_skgc_realization.exec();
                check_rkgc_realization.exec();
                update_body_kinetic_energy.exec();
                body_average_kinetic_energy = calculate_body_average_kinetic_energy.exec();
                body_maximum_kinetic_energy = calculate_body_maximum_kinetic_energy.exec();
                ac_term = calculate_body_average_ac.exec();
                as_term = calculate_body_average_as.exec();
                c_term = calculate_body_average_c.exec();
                ab_term = calculate_body_average_ab.exec();
                ar_term = calculate_body_average_ar.exec();
                b_term = calculate_body_average_b.exec();
                std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
                std::cout << "Body: "
                    << " Average: " << body_average_kinetic_energy
                    << " Maximum: " << body_maximum_kinetic_energy << std::endl;
                out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " << 
                    ac_term << " " << as_term << " " << c_term << " " <<
                    ab_term << " " << ar_term << " " << b_term << "\n";
                write_body_to_vtp.writeToFile(ite);
            }
        }
        ite++;

        write_body_to_vtp.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);
        TickCount t2 = TickCount::now();
        TickCount::interval_t tt;
        tt = t2 - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    /* Body relation. */
    InnerRelation body_inner(body);
    /** Write the body state to Vtp file. */
    BodyStatesRecordingToVtp write_body_to_vtp(io_environment, { &body });
    /** Random reset the insert body particle position. */
    SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
    /* Kernel correction matrix. */
    InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> correction_matrix(body_inner);
    body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    /* Initialize the field. */
    SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
    /* Update reduce error. */
    InteractionDynamics<relax_dynamics::CheckL2NormError> check_l2_error(body_inner);
    ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
        calculate_nkgc_error(average_domain, "NKGCError");
    ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
        calculate_skgc_error(average_domain, "SKGCError");
    ReduceDynamics<Average<QuantitySummationPartly<Real, BodyPartByParticle>>>
        calculate_ckgc_error(average_domain, "CKGCError");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correction_matrix.exec();
    testing_initial_condition.exec();
    //random_insert_body_particles.exec(0.1);
    //----------------------------------------------------------------------
    //	Setup for error initialization
    //----------------------------------------------------------------------
    Real nkgc_error = 0.0;
    Real skgc_error = 0.0;
    Real ckgc_error = 0.0;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    std::string filefullpath_error_information = io_environment.output_folder_ + "/" + "error_information.dat";
    std::ofstream out_file_nonopt_temperature(filefullpath_error_information.c_str(), std::ios::app);
    check_l2_error.exec();
    nkgc_error = calculate_nkgc_error.exec();
    skgc_error = calculate_skgc_error.exec();
    ckgc_error = calculate_ckgc_error.exec();

    std::cout << "nkgc error is " << nkgc_error << std::endl;
    std::cout << "skgc error is " << skgc_error << std::endl;
    std::cout << "ckgc error is " << ckgc_error << std::endl;
    write_body_to_vtp.writeToFile(0);

    out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "nkgc error is " << "   " << nkgc_error << "\n";
    out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "skgc error is " << "   " << skgc_error << "\n";
    out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "ckgc error is " << "   " << ckgc_error << "\n";
    out_file_nonopt_temperature.close();

    return 0;
};

