/**
 * @file 	particle_generator_single_resolution.cpp
 * @brief 	This is the test of using level set to generate particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by level set for a complex geometry (2D).
 *			Before particle generation, we clean the level set, then do re-initialization.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"

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
BoundingBoxd system_domain_bounds(Vec2d(-DL1, -0.5), Vec2d(DL, DH));
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
    LevelSetShape *level_set_shape = input_body.defineBodyLevelSetShape(par_ck, 2.0)
                                         ->addPackageVariableToWrite<Real>("KernelWeight")
                                         ->addCellVariableToWrite<UnsignedInt>("CellPackageIndex")
                                         ->addCellVariableToWrite<int>("CellContainID")
                                         ->writeLevelSet(sph_system);
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
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    StdVec<RealBody *> real_bodies = {&input_body, &filler};

    host_methods.addStateDynamics<RandomizeParticlePositionCK>(real_bodies).exec(); // host method able to run immediately

    ParticleDynamicsGroup update_cell_linked_list = main_methods.addCellLinkedListDynamics(real_bodies);
    ParticleDynamicsGroup update_relation;
    update_relation.add(&main_methods.addRelationDynamics(input_body_inner));
    update_relation.add(&main_methods.addRelationDynamics(filler_inner, filler_contact));
    ParticleDynamicsGroup update_configuration = update_cell_linked_list + update_relation;

    ParticleDynamicsGroup relaxation_residual;
    relaxation_residual.add(&main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(input_body_inner)
                                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(input_body, *level_set_shape));
    relaxation_residual.add(&main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(filler_inner)
                                 .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(filler_contact));

    ReduceDynamicsGroup relaxation_scaling = main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(real_bodies);

    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(real_bodies);
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(near_body_surface));
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    BaseCellLinkedList &input_body_cell_linked_list = input_body.getCellLinkedList();
    input_body_cell_linked_list.addCellVariableToWrite<UnsignedInt>("CurrentListSize");
    MeshRecordingToPlt cell_linked_list_recording(sph_system, input_body_cell_linked_list);
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
        update_configuration.exec();
        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);

        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            body_state_recorder.writeToFile(ite_p);
            cell_linked_list_recording.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;

    return 0;
}
