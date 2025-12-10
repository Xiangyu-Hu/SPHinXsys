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
    RealBody input_body(sph_system, makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
    // level set shape is used for particle relaxation
    LevelSetShape *level_set_shape = input_body.defineBodyLevelSetShape(par_ck)
                                         ->correctLevelSetSign()
                                         ->writeLevelSet(sph_system);
    input_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    NearShapeSurface near_body_surface(input_body);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    Inner<> input_body_inner(input_body);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &input_body_cell_linked_list = main_methods.addCellLinkedListDynamics(input_body);
    auto &input_body_update_inner_relation = main_methods.addRelationDynamics(input_body_inner);
    auto &random_input_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(input_body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(input_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(input_body, *level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(input_body);
    auto &update_particle_position =
        main_methods.addStateDynamics<PositionRelaxationCK>(input_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(input_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    random_input_body_particles.exec();

    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile(0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        input_body_cell_linked_list.exec();
        input_body_update_inner_relation.exec();

        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();

        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            body_state_recorder.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;

    return 0;
}
