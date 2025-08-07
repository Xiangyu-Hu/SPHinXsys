/**
 * @file 	particle_generator_single_resolution.cpp
 * @brief 	This is the test of using level set to generate particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by level set for a complex geometry (2D).
 *			Before particle generation, we clean the level set, then do re-initialization.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys_ck.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string input_body = "./input/SPHinXsys-2d.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Real DL = 2.5;                          /**< InputBody length right part. */
Real DL1 = 2.5;                         /**< InputBody length left part. */
Real DH = 5.0;                          /**< InputBody height. */
Real resolution_ref = (DL + DL1) / 120; /**< Reference resolution. */
BoundingBox system_domain_bounds(Vec2d(-DL1, -0.5), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    MultiPolygon original_logo;
    original_logo.addAPolygonFromFile(input_body, ShapeBooleanOps::add);
    ComplexShape input_shape("SPHInXsysLogo");
    input_shape.add<ExtrudeShape<MultiPolygonShape>>(4.0 * resolution_ref, original_logo);
    input_shape.subtract<MultiPolygonShape>(original_logo);
    RealBody input_body(sph_system, input_shape);
    LevelSetShape *level_set_shape =
        input_body.defineBodyLevelSetShape(2.0)
            ->addMeshVariableToWrite<Real>("KernelWeight")
            ->writeLevelSet(sph_system)
            ->addBackgroundMeshVariableToWrite<UnsignedInt>("CellPackageIndex")
            ->writeBackgroundMesh(sph_system);
    input_body.generateParticles<BaseParticles, Lattice>();

    MultiPolygonShape filler_shape(original_logo, "Filler");
    RealBody filler(sph_system, filler_shape);
    filler.generateParticles<BaseParticles, Lattice>();
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
    Inner<> filler_inner(filler);
    Contact<> filler_contact(filler, {&input_body});
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par);
    auto &host_methods = sph_solver.addParticleMethodContainer(par);
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
    auto &filler_cell_linked_list = host_methods.addCellLinkedListDynamics(filler);
    auto &input_body_update_inner_relation = main_methods.addRelationDynamics(input_body_inner);
    auto &filler_update_complex_relation = host_methods.addRelationDynamics(filler_inner, filler_contact);

    auto &random_input_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(input_body);
    auto &random_filler_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(filler);

    auto &relaxation_residue =
        main_methods.addInteractionDynamics<RelaxationResidueCK, NoKernelCorrectionCK>(input_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(input_body, *level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(input_body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(input_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);

    auto &filler_relaxation_residue =
        main_methods.addInteractionDynamics<RelaxationResidueCK, NoKernelCorrectionCK>(filler_inner)
            .addContactInteraction<Boundary, NoKernelCorrectionCK>(filler_contact);
    auto &filler_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(filler);
    auto &filler_update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(filler);

    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    random_input_body_particles.exec();
    random_filler_particles.exec();
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
        filler_cell_linked_list.exec();
        input_body_update_inner_relation.exec();
        filler_update_complex_relation.exec();

        relaxation_residue.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();

        filler_relaxation_residue.exec();
        Real filler_relaxation_step = filler_relaxation_scaling.exec();
        filler_update_particle_position.exec(filler_relaxation_step);

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
