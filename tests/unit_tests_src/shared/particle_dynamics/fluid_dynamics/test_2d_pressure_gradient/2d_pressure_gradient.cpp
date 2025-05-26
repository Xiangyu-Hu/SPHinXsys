#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real mu_f   = 1.0e-1;      /**< Viscosity. */
Real U_max  = 1.0;         // make sure the maximum anticipated speed
Real c_f    = 10.0 * U_max; /**< Reference sound speed. */
//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real width            = 1.0;
Real height           = 0.5;
Real particle_spacing = 0.01;
Real boundary_width   = particle_spacing * 4; // boundary width
//----------------------------------------------------------------------
//  Google test item.
//----------------------------------------------------------------------
Real min_computed(0.0);
Real max_computed(0.0);
Real reference = 2.0;
TEST(PressureGradient, MaxErrorNorm)
{
    Real max_error = ABS(min_computed - reference) / reference;
    max_error = SMAX(max_error, ABS(max_computed - reference) / reference);
    EXPECT_LT(max_error, 0.05);
    std::cout << "Reference PressureGradientNorm: " << reference
              << " and MaxErrorNorm: " << max_error << std::endl;
}
//----------------------------------------------------------------------
//  Complex shapes for wall boundary
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
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_position), scaled_container);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container_outer(0.5 * width + boundary_width,
                                    0.5 * height + boundary_width);
        Vecd scaled_container(0.5 * width + 2.0 * boundary_width,
                              0.5 * height);
        Transform translate_to_origin_outer(
            Vec2d(-boundary_width, -boundary_width) + scaled_container_outer);
        Transform translate_to_origin_inner(
            Vec2d(-boundary_width, 0.0) + scaled_container);

        add<TransformShape<GeometricShapeBox>>(
            Transform(translate_to_origin_outer), scaled_container_outer);
        subtract<TransformShape<GeometricShapeBox>>(
            Transform(translate_to_origin_inner), scaled_container);
    }
};

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(scaled_container);
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin),
                                               scaled_container);
    }
};

//----------------------------------------------------------------------
//  Application dependent initial condition
//----------------------------------------------------------------------
class LinearPressureProfile : public fluid_dynamics::FluidInitialCondition
{
  public:
    Real* p_;

    explicit LinearPressureProfile(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body)
    {
        p_ = particles_->template getVariableDataByName<Real>("Pressure");
    }

    void update(size_t index_i, Real dt)
    {
        Real y = pos_[index_i][1];  // Use y coordinate like velocity gradient
        p_[index_i] = y / height;   // Linear profile: p(y) = y/h (like v(y) = y/h in velocity)
    }
};

class BoundaryPressure : public BodyPartMotionConstraint
{
  public:
    Real* p_;

    explicit BoundaryPressure(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part)
    {
        p_ = particles_->template getVariableDataByName<Real>("Pressure");
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        p_[index_i] = 1.0;  // Fixed pressure at boundary, like fixed velocity in velocity test
    }
};

//----------------------------------------------------------------------
//  Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    BoundingBox system_domain_bounds(
        Vecd(-boundary_width * 2, -boundary_width * 2),
        Vecd(width + boundary_width * 2, height + boundary_width * 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.setIOEnvironment();

    // Creating bodies with materials and particles
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(
        ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    water_block.getBaseParticles().registerStateVariable<Real>("Pressure");
    water_block.getBaseParticles().registerStateVariable<Vecd>("PressureGradient");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    wall_boundary.getBaseParticles().registerStateVariable<Real>("Pressure");

    // Body relation map
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    // Numerical methods
    SimpleDynamics<LinearPressureProfile> initial_condition(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(water_wall_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::PressureGradientWithWall<LinearGradientCorrection>> pressure_grad_calculation(water_block_inner, water_wall_contact);
    
    BodyRegionByParticle upper_wall(wall_boundary, makeShared<UpperBoundary>("UpperWall"));
    SimpleDynamics<BoundaryPressure> upper_wall_pressure(upper_wall);
    
    ReduceDynamics<VariableNorm<Vecd, ReduceMax>> maximum_pressure_gradient_norm(water_block, "PressureGradient");
    ReduceDynamics<VariableNorm<Vecd, ReduceMin>> minimum_pressure_gradient_norm(water_block, "PressureGradient");

    // I/O and regression test setup
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(water_block, "PressureGradient");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");

    // Prepare and execute
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    upper_wall_pressure.exec();
    wall_boundary_normal_direction.exec();
    distance_to_wall.exec();
    corrected_configuration_fluid.exec();
    pressure_grad_calculation.exec();

    body_states_recording.writeToFile(0);
    max_computed = maximum_pressure_gradient_norm.exec();
    min_computed = minimum_pressure_gradient_norm.exec();

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
