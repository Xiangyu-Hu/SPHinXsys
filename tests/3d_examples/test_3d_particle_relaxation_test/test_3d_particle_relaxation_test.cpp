/**
 * @file 	particle_relaxation_test.cpp
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

std::string full_path_to_file = "./input/gear3.2.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-0.1, -0.1, -0.1);
Vec3d domain_upper_bound(0.1, 0.1, 0.1);
Vecd translation(0.0, 0.0, 0.0);
Real scaling = 1.0;
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
Real dp_0 = 0.002;
//--------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
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
    // level set shape is used for particle relaxation
    imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet()->writeLevelSet(sph_system);
    imported_model.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp({imported_model});
    MeshRecordingToPlt write_cell_linked_list(sph_system, imported_model.getCellLinkedList());
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
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(imported_model_inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    write_imported_model_to_vtp.writeToFile(0.0);
    imported_model.updateCellLinkedList();
    write_cell_linked_list.writeToFile(0.0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
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
