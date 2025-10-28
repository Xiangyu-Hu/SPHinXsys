/**
 * @file 	load_image.cpp
 * @brief 	This is the test of using distance map to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yijin Mao, Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
#include "image_shape.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_image = "./input/sphere.mhd";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-25.0, -25.0, -25.0);
Vec3d domain_upper_bound(25.0, 25.0, 25.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 50.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ImageShapeFromFile>(full_path_to_image);
    }
};
//----------------------------------------------------------------------
//	Main program begins here
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody imported_model(sph_system, makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
    imported_model.defineAdaptation<ParticleRefinementNearSurface>(1.15, 1.0, 2);
    imported_model.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    imported_model.generateParticles<BaseParticles, Lattice, Adaptive>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp({imported_model});
    write_imported_model_to_vtp.addToWrite<Real>(imported_model, "SmoothingLengthRatio");
    MeshRecordingToPlt cell_linked_list_recording(sph_system, imported_model.getCellLinkedList());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    AdaptiveInnerRelation imported_model_inner(imported_model);
    // BaseInnerRelation* imported_model_inner(imported_model);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
    /** A  Physics relaxation step. */
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(imported_model_inner);
    SimpleDynamics<UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(imported_model);
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
        // update_smoothing_length_ratio.exec();
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
