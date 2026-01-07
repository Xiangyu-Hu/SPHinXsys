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
Real global_resolution = (domain_upper_bound[0] - domain_lower_bound[0]) / 12.5;
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
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
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &imported_model = sph_system.addAdaptiveBody<RealBody>(
        AdaptiveNearSurface(global_resolution, 1.15, 1.0, 3), makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
    LevelSetShape *level_set_shape =
        imported_model.defineBodyLevelSetShape()->cleanLevelSet()->writeLevelSet();
    imported_model.generateParticles<BaseParticles, Lattice>();
    auto &near_body_surface = imported_model.addBodyPart<NearShapeSurface>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &imported_model_inner = sph_system.addInnerRelation(imported_model);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    host_methods.addStateDynamics<RandomizeParticlePositionCK>(imported_model).exec();
    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(imported_model);
    auto &update_inner_relation = main_methods.addRelationDynamics(imported_model_inner);

    auto &relaxation_residual =
        main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(imported_model_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(imported_model, *level_set_shape);

    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(imported_model);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    auto &update_smoothing_length_ratio = main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(imported_model, *level_set_shape);

    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(imported_model);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(imported_model, "SmoothingLengthRatio");
    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        update_cell_linked_list.exec();
        update_inner_relation.exec();

        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();
        update_smoothing_length_ratio.exec();

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
