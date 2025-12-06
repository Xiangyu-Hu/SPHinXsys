/**
 * @file 	test_3d_shell_particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by levelset for a complex thin structures geometry (3D).
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_geometry = "./input/curved_tube.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(12, 14, 446);
Vec3d domain_upper_bound(1315, 1317, 1302);
Real dp_0 = 25.0;
Real thickness = 50.0;
// level set resolution much higher than that of particles is required
Real level_set_refinement_ratio = dp_0 / (0.1 * thickness);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ImportedShellModel : public ComplexShape
{
  public:
    explicit ImportedShellModel(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_geometry, Vecd::Zero(), 1.0);
    }
};
//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody imported_model(sph_system, makeShared<ImportedShellModel>("ImportedShellModel"));
    imported_model.defineBodyLevelSetShape(level_set_refinement_ratio, UsageType::Surface)
        ->writeLevelSet(sph_system);
    imported_model.generateParticles<SurfaceParticles, Lattice>(thickness);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp({imported_model});
    write_imported_model_to_vtp.addToWrite<Vecd>(imported_model, "NormalDirection");
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation imported_model_inner(imported_model);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
    /** A  Physics relaxation step. */
    ShellRelaxationStep relaxation_step_inner(imported_model_inner);
    ShellNormalDirectionPrediction shell_normal_prediction(imported_model_inner, thickness);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_step_inner.MidSurfaceBounding().exec();
    write_imported_model_to_vtp.writeToFile(0.0);
    imported_model.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
            write_imported_model_to_vtp.writeToFile(ite_p);
        }
        relaxation_step_inner.exec();
        ite_p += 1;
    }
    shell_normal_prediction.exec();
    write_imported_model_to_vtp.writeToFile(ite_p);
    std::cout << "The physics relaxation process of imported model finish !" << std::endl;

    return 0;
}
