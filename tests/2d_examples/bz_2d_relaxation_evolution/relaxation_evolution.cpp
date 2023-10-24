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
Real LL = 1.0;
Real LH = 1.0;
Real resolution_ref = LH / 50.0;
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(LL, LH));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d water_block_translation = water_block_halfsize;

class Insert : public ComplexShape
{
public:
	explicit Insert(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform(water_block_halfsize), water_block_translation);
	}
};

class TestingInitialCondition
	: public fluid_dynamics::FluidInitialCondition
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
		scalar_[index_i] = pos_[index_i][0];
        vector_[index_i] = pos_[index_i];
        matrix_[index_i] = Vec2d(1,1) * pos_[index_i].transpose();
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
    FluidBody body(sph_system, makeShared<Insert>("WaterBody"));
    body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1, 1, 1);
    body.addBodyStateForRecording<Vecd>("Position");
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
        : body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map.
        //----------------------------------------------------------------------
        InnerRelation insert_body_inner(body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
        BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, {&body});
        ReloadParticleIO write_particle_reload_files(io_environment, {&body});

        /* Relaxation method: including based on the 0th and 1st order consistency. */
        InteractionWithUpdate<KernelCorrectionMatrixInner> calculate_correction_matrix(insert_body_inner, false);
        relax_dynamics::RelaxationStepInner relaxation_inner(insert_body_inner, false);
        relax_dynamics::RelaxationStepInnerImplicit relaxation_inner_implicit(insert_body_inner, false);
       
        /* Update the relaxation residual. */
        SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
        InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_consistency_realization(insert_body_inner, false);
        body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        testing_initial_condition.exec();
        write_insert_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        TickCount t1 = TickCount::now();
        int ite = 0; //iteration step for the total relaxation step.

        std::string filefullpath_all_information = io_environment.output_folder_ + "/" + "all_information.dat";
        std::ofstream out_file_all_information(filefullpath_all_information.c_str(), std::ios::app);

        GlobalStaticVariables::physical_time_ = ite;
        /* The procedure to obtain uniform particle distribution that satisfies the 0th order consistency. */
        while (ite < 500)
        {
            calculate_correction_matrix.exec();
            relaxation_inner.exec();
            ite++;

            if (ite % 100 == 0)
            {
                testing_initial_condition.exec();
                calculate_correction_matrix.exec();
                check_consistency_realization.exec();
                std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
                write_insert_body_to_vtp.writeToFile(ite);
            }
        }

        ite++;
        write_insert_body_to_vtp.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);
        TickCount t2 = TickCount::now();
        TickCount::interval_t tt;
        tt = t2 - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        return 0;
    }
};

