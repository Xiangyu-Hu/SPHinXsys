/**
 * @file 	2d_gradient.cpp
 * @brief 	test the linear and parabolic reproducing gradient and hessian
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h"
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
//	Geometric shapes used in the test
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(container);
        add<GeometricShapeBox>(Transform(translate_to_origin), container);
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

        add<GeometricShapeBox>(Transform(translate_to_origin_outer), container_outer);
        subtract<GeometricShapeBox>(Transform(translate_to_origin_inner), container);
    }
};
//----------------------------------------------------------------------
//	Google test items
//----------------------------------------------------------------------
Vec2d random_observation(rand_uniform(0.0, width), rand_uniform(0.0, height));

Matd approximated_gradient(Matd::Zero());
Matd reference_gradient(Matd::Identity());
TEST(LinearGradient, Error)
{
    EXPECT_LT((reference_gradient - approximated_gradient).norm(), 1.0e-5);
    std::cout << "Reference Gradient: " << reference_gradient << " and "
              << "Predicted Gradient: " << approximated_gradient << std::endl;
};

Vec2d first_coefficient(rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0));
VecMat2d second_coefficient{rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0)};
class ParabolicProfile : public ReturnFunction<Real>
{
    Vec2d first_coefficient_;
    VecMat2d second_coefficient_;

  public:
    ParabolicProfile() : first_coefficient_(first_coefficient),
                         second_coefficient_(second_coefficient) {};

    Real operator()(const Vec2d &position)
    {
        return first_coefficient_.dot(position) +
               second_coefficient_.dot(vectorizeTensorSquare(position));
    }
};

VecMat2d approximated_hessian = VecMat2d::Zero();
VecMat2d reference_hessian = 2.0 * second_coefficient;
TEST(Hessian, Error)
{
    EXPECT_LT((reference_hessian - approximated_hessian).norm(), 1.0e-4);
    std::cout << "Reference Hessian: " << reference_hessian << " and "
              << "Predicted Hessian: " << approximated_hessian << std::endl;
};

Vec2d approximated_2nd_order_gradient = Vec2d::Zero();
Vec2d reference_2nd_order_gradient =
    first_coefficient +
    Vec2d(2.0 * random_observation[0] * second_coefficient[0] + random_observation[1] * second_coefficient[2],
          2.0 * random_observation[1] * second_coefficient[1] + random_observation[0] * second_coefficient[2]);
TEST(SecondOrderGradient, Error)
{
    EXPECT_LT((reference_2nd_order_gradient - approximated_2nd_order_gradient).norm(), 1.0e-5);
    std::cout << "Reference Second Order Gradient: " << reference_2nd_order_gradient << " and "
              << "Predicted Second Order Gradient: " << approximated_2nd_order_gradient << std::endl;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-boundary_width * 2, -boundary_width * 2),
                                     Vecd(width + boundary_width * 2, height + boundary_width * 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_fluid_particles(water_block);
    random_fluid_particles.exec(0.25); // randomize particle to avoid the symmetry

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vec2d> observation_location = {random_observation};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall});
    Contact<> fluid_observer_contact(fluid_observer, {&water_block});
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
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_water_block_inner(water_block_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_water_wall_contact(water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_fluid_observer_contact(fluid_observer_contact);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        linear_correction_matrix(DynamicsArgs(water_block_inner, 0.0), water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>, Contact<Vecd>>>
        position_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Position")),
            DynamicsArgs(water_wall_contact, std::string("Position")));
    ObservedQuantityRecording<MainExecutionPolicy, Matd, RestoringCorrection>
        observed_position_gradient("PositionGradient", fluid_observer_contact);

    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>, Contact<>>>
        displacement_matrix_gradient(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        hessian_correction_matrix(DynamicsArgs(water_block_inner, 0.0), water_wall_contact);

    StateDynamics<MainExecutionPolicy, VariableAssignment<SpatialDistribution<ParabolicProfile>, SPHBody>>
        water_block_initial_condition(water_block, "Phi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SpatialDistribution<ParabolicProfile>, SPHBody>>
        wall_initial_condition(wall, "Phi");
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        variable_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>, Contact<Real>>>
        variable_hessian(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian("PhiHessian", fluid_observer_contact);

    InteractionDynamicsCK<MainExecutionPolicy, SecondOrderGradient<Inner<Real>, Contact<Real>>>
        variable_2nd_order_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    ObservedQuantityRecording<MainExecutionPolicy, Vec2d, RestoringCorrection>
        observed_2nd_order_gradient("PhiGradient", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    update_water_block_inner.exec();
    update_water_wall_contact.exec();
    update_fluid_observer_contact.exec();

    water_block_initial_condition.exec();
    wall_initial_condition.exec();

    linear_correction_matrix.exec();
    position_linear_gradient.exec();
    observed_position_gradient.writeToFile(0);
    approximated_gradient = *observed_position_gradient.getObservedQuantity();

    displacement_matrix_gradient.exec();
    hessian_correction_matrix.exec();
    variable_linear_gradient.exec();
    variable_hessian.exec();
    observed_hessian.writeToFile(0);
    approximated_hessian = *observed_hessian.getObservedQuantity();

    variable_2nd_order_gradient.exec();
    observed_2nd_order_gradient.writeToFile(0);
    approximated_2nd_order_gradient = *observed_2nd_order_gradient.getObservedQuantity();
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
