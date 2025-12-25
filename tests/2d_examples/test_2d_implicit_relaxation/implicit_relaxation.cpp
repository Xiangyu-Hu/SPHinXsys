/**
 * @file 	implicit_relaxation.cpp
 * @brief 	This is the first case by testing the relaxation with implicit method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string input_body = "./input/TurbineBlade.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 10.0;
Real DH = 10.0;
Real global_resolution = 1 / 25.0;
Real BW = global_resolution * 4;
BoundingBoxd system_domain_bounds(Vec2d(-BW - DL, -BW - DH), Vec2d(BW + DL, BW + DH));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
class InputBody : public ComplexShape
{
  public:
    explicit InputBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon turbine_blade;
        turbine_blade.addAPolygonFromFile(input_body, ShapeBooleanOps::add);
        add<MultiPolygonShape>(turbine_blade);
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody insert_body(sph_system, makeShared<InputBody>("Body"));
    insert_body.defineBodyLevelSetShape()->writeLevelSet();
    insert_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map used for particle relaxation.
    //----------------------------------------------------------------------
    InnerRelation insert_body_inner(insert_body);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(insert_body);
    RelaxationStepLevelSetCorrectionInnerImplicit relaxation_step_inner(insert_body_inner);
    BodyStatesRecordingToVtp write_insert_body_to_vtp(insert_body);
    ReloadParticleIO write_particle_reload_files(insert_body);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<Average<QuantitySummation<Real>>>>
        write_average_kinetic_energy(insert_body, "ParticleKineticEnergy");
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_insert_body_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    write_insert_body_to_vtp.writeToFile(0);
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
            write_insert_body_to_vtp.writeToFile(ite_p);
            write_average_kinetic_energy.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    /** Output results. */
    write_particle_reload_files.writeToFile(0);

    if (sph_system.GenerateRegressionData())
    {
        write_average_kinetic_energy.generateDataBase(0.05);
    }
    else
    {
        write_average_kinetic_energy.testResult();
    }
    return 0;
};
