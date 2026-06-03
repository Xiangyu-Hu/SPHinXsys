/**
 * @file 	sdf_shape_particle_relaxation.cpp
 * @brief 	This is the test of using signed distance function primitives to generate shape for particle relaxation.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//-----------------------------------------------------------------------------------------------------------
//	Geometric elements.
//-----------------------------------------------------------------------------------------------------------
BoundingBoxd system_domain_bounds(Vec3d::Constant(2.0));
Real spacing_ref = system_domain_bounds.MinimumDimension() / Real(10);
AdaptiveNearSurface adaptive_near_surface(spacing_ref, 1.15, 1.0, 3);
SDFBall sdf_ball(0.8);
SDFCappedCone sdf_capped_cone(1.0, 1.0, 0.5);
SDFTransform sdf_transform(Rotation3d(Pi / 4.0, Vec3d::UnitY()), Vec3d(-0.5, 0.0, 0.0));
SDFExtension sdf_extend(sdf_capped_cone, sdf_transform);
SDFSmoothAddition sdf_addition(adaptive_near_surface.MinimumSpacing());
SDFOperation sdf_operation(sdf_addition, sdf_ball, sdf_extend);
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &sdf_shape = sph_system.addShape<SDFShape>(adaptive_near_surface.MinimumSpacing(), "SDFShape");
    sdf_shape.insertSDFPrimitive("Ball", sdf_operation, GeometricOps::add);
    auto &my_body = sph_system.addAdaptiveBody<RealBody>(adaptive_near_surface, sdf_shape);
    LevelSetShape &level_set_shape = my_body.defineBodyLevelSetShape().writeLevelSet();
    my_body.generateParticles<BaseParticles, Lattice>();
    auto &near_body_surface = my_body.addBodyPart<NearShapeSurface>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &my_body_inner = sph_system.addInnerRelation(my_body);
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
    host_methods.addStateDynamics<RandomizeParticlePositionCK>(my_body).exec();
    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(my_body);
    auto &update_inner_relation = main_methods.addRelationDynamics(my_body_inner);

    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(my_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(my_body, level_set_shape);

    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(my_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    auto &update_smoothing_length_ratio = main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(my_body, level_set_shape);

    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(my_body);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(my_body, "SmoothingLengthRatio");
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
