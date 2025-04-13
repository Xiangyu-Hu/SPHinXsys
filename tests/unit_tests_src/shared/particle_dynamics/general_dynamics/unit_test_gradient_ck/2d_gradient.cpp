/**
 * @file 	2d_gradient.cpp
 * @brief 	test the linear reproducing approximation of interpolation.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys_ck.h"
#include <gtest/gtest.h>
using namespace SPH;
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real mu_f = 1.0e-1;      /**< Viscosity. */
Real U_max = 1.0;        // make sure the maximum anticipated speed
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real width = 1.0;
Real height = 0.5;
Real particle_spacing = 0.01;
Real boundary_width = particle_spacing * 4; // boundary width
//----------------------------------------------------------------------
//	Google test item.
//----------------------------------------------------------------------
Matd approximated_gradient(Matd::Zero());
Matd reference_gradient(Matd::Identity());
TEST(LinearGradient, Error)
{
    EXPECT_LT((reference_gradient - approximated_gradient).norm(), 1.0e-6);
    std::cout << "Reference Gradient: " << reference_gradient << " and "
              << "Predicted Gradient: " << approximated_gradient << std::endl;
};

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(container);
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin), container);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd container_outer(0.5 * width + boundary_width, 0.5 * height + boundary_width);
        Vecd container(0.5 * width + 2.0 * boundary_width, 0.5 * height);
        Transform translate_to_origin_outer(Vec2d(-boundary_width, -boundary_width) + container_outer);
        Transform translate_to_origin_inner(Vec2d(-boundary_width, 0.0) + container);

        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin_outer), container_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin_inner), container);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-boundary_width * 2, -boundary_width * 2),
                                     Vecd(width + boundary_width * 2, height + boundary_width * 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_airfoil_particles(water_block);
    random_airfoil_particles.exec(0.25);

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    Vecd random_coordinate(rand_uniform(0.0, width), rand_uniform(0.0, height));
    StdVec<Vecd> observation_location = {random_coordinate};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Relation<Inner<>> water_block_inner(water_block);
    Relation<Contact<>> water_wall_contact(water_block, {&wall});
    Relation<Contact<>> fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
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
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_water_block_inner(water_block_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_water_wall_contact(water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_fluid_observer_contact(fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        fluid_linear_correction_matrix(DynamicsArgs(water_block_inner, 0.5), water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>, Contact<Vecd>>>
        position_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Position")),
            DynamicsArgs(water_wall_contact, std::string("Position")));
    ObservedQuantityRecording<MainExecutionPolicy, Matd, RestoringCorrection>
        fluid_observer_position("PositionGradient", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    update_water_block_inner.exec();
    update_water_wall_contact.exec();
    update_fluid_observer_contact.exec();

    fluid_linear_correction_matrix.exec();
    position_linear_gradient.exec();

    approximated_gradient = *fluid_observer_position.getObservedQuantity();
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
