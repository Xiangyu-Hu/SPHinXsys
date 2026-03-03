/**
 * @file 	velocity_gradient.cpp
 * @brief 	test the approximation of velocity gradient.
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
//	Google test item.
//----------------------------------------------------------------------
Real min_computed(0.0);
Real max_computed(0.0);
Real reference = 2.0;
TEST(VelocityGradient, MaxErrorNorm)
{
    Real max_error = ABS(min_computed - reference) / reference;
    max_error = SMAX(max_error, ABS(max_computed - reference) / reference);
    EXPECT_LT(max_error, 0.05);
    std::cout << "Reference VelocityGradientNorm: " << reference << " and "
              << "MaxErrorNorm: " << max_error << std::endl;
}
//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class UpperBoundary : public ComplexShape
{
  public:
    explicit UpperBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width + boundary_width, 0.5 * boundary_width);
        Transform translate_to_origin(scaled_container);
        Vecd transform(-boundary_width, height);
        Transform translate_to_position(transform + scaled_container);
        add<GeometricShapeBox>(Transform(translate_to_position), scaled_container);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container_outer(0.5 * width + boundary_width, 0.5 * height + boundary_width);
        Vecd scaled_container(0.5 * width + 2.0 * boundary_width, 0.5 * height);
        Transform translate_to_origin_outer(Vec2d(-boundary_width, -boundary_width) + scaled_container_outer);
        Transform translate_to_origin_inner(Vec2d(-boundary_width, 0.0) + scaled_container);

        add<GeometricShapeBox>(Transform(translate_to_origin_outer), scaled_container_outer);
        subtract<GeometricShapeBox>(Transform(translate_to_origin_inner), scaled_container);
    }
};
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(scaled_container);
        add<GeometricShapeBox>(Transform(translate_to_origin), scaled_container);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class CouetteFlowInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit CouetteFlowInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        Vecd velocity = ZeroData<Vecd>::value;
        velocity[0] = pos_[index_i][1] / height;
        vel_[index_i] = velocity;
    }
};

class BoundaryVelocity : public BodyPartMotionConstraint
{
  public:
    explicit BoundaryVelocity(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd velocity = ZeroData<Vecd>::value;
        velocity[0] = 1.0;
        vel_[index_i] = velocity;
    };
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

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<CouetteFlowInitialCondition> initial_condition(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(water_wall_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> vel_grad_calculation(water_block_inner, water_wall_contact);
    BodyRegionByParticle upper_wall(wall_boundary, makeShared<UpperBoundary>("UpperWall"));
    SimpleDynamics<BoundaryVelocity> upper_wall_velocity(upper_wall);
    ReduceDynamics<VariableNorm<Matd, ReduceMax>> maximum_velocity_gradient_norm(water_block, "VelocityGradient");
    ReduceDynamics<VariableNorm<Matd, ReduceMin>> minimum_velocity_gradient_norm(water_block, "VelocityGradient");

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Matd>(water_block, "VelocityGradient");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    upper_wall_velocity.exec();
    wall_boundary_normal_direction.exec();
    distance_to_wall.exec();
    corrected_configuration_fluid.exec();
    vel_grad_calculation.exec();

    body_states_recording.writeToFile(0);
    max_computed = maximum_velocity_gradient_norm.exec(),
    min_computed = minimum_velocity_gradient_norm.exec(),
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
