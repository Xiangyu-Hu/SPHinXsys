/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the first geometry.
//	To use this, please commenting the setting for the second geometry.
//----------------------------------------------------------------------
// std::string full_path_to_file = "./input/SPHinXsys.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
/*Vec3d domain_lower_bound(-2.3, -0.1, -0.3);
Vec3d domain_upper_bound(2.3, 4.5, 0.3);
Vecd translation(0.0, 0.0, 0.0);
Real scaling = 1.0; */

//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/triangle_prism.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.75; // Domain width.
Real DH = 1.3;  // Domain height.
Real DW = 1.2;  // Domain width.
Vec3d domain_lower_bound(-0.2, -0.2, -0.2);
Vec3d domain_upper_bound(DL, DH, DW);
Vecd translation(0.5 * DL, 0.5 * DH, 0.5 * DW);
Real scaling = 2.5;
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
//--------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(4.0 * dp_0, full_path_to_file, translation, scaling);
        subtract<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
// std::string full_path_to_file = "./input/fluid.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
// Vec3d domain_lower_bound(-0.036, -0.046, -0.011);
// Vec3d domain_upper_bound(0.036, 0.093, 0.072);
// Vecd translation(0.0, 0.0, 0.0);
// Real scaling = 1.0;

//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
// BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
// Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
//----------------------------------------------------------------------
//	define the imported model.
//----------------------------------------------------------------------
// class SolidBodyFromMesh : public ComplexShape
//{
//  public:
//    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
//    }
//};
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
    imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet();
    imported_model.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp({imported_model});
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
