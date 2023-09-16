/**
 * @file 	particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/teapot.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-9.0, -6.0, 0.0);
Vec3d domain_upper_bound(9.0, 6.0, 9.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 12.5;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	define a body from the imported model.
//----------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, dp_0);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody imported_model(system, makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
    imported_model.defineAdaptation<ParticleRefinementNearSurface>(1.15, 1.0, 3);
    imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment);
    imported_model.defineParticlesAndMaterial();
    imported_model.generateParticles<ParticleGeneratorMultiResolution>();
    imported_model.addBodyStateForRecording<Real>("SmoothingLengthRatio");
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp(io_environment, {imported_model});
    MeshRecordingToPlt cell_linked_list_recording(io_environment, imported_model.getCellLinkedList());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    AdaptiveInnerRelation imported_model_inner(imported_model);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
    /** A  Physics relaxation step. */
    relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
    SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(imported_model);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    update_smoothing_length_ratio.exec();
    write_imported_model_to_vtp.writeToFile();
    imported_model.updateCellLinkedList();
    cell_linked_list_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        update_smoothing_length_ratio.exec();
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
            write_imported_model_to_vtp.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process of imported model finish !" << std::endl;

    return 0;
}
