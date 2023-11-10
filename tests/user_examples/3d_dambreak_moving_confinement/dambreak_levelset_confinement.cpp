/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
#include "body_part_by_cell_tracing.h"
#include "level_set_confinement.h"
#include "math.h"
#include "time_step_size_moving_velocity.h"
using namespace SPH;

// general parameters for geometry
Real resolution_ref = 0.5;   // particle spacing
Real BW = resolution_ref * 4; // boundary width
Real DL = 30.0;              // tank length
Real DH = 50.0;                // tank height
Real DW = 30.0;                // tank width
Real LL = 30.0;                // liquid length
Real LH = 30.0;                // liquid height
Real LW = 30.0;                // liquid width
Real CL = 5.0;                // liquid length
Real CH = 5.0;                // liquid height
Real CW = 5.0;                // liquid width

// for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

Vecd axis_point_1 (0.0, 0.0, -0.5 * LW);
Vecd axis_point_2 (0.0, 0.0, 0.5 * LW);

Vecd ball_center(0.0, 0.0, 0.0);
Real ball_radius = 5.0;
Real angular_velocity = 0.2 * Pi;
//	define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_water(halfsize_water);
        Vecd halfsize_cubic(0.5 * CL, 0.5 * CH, 0.5 * CW);
        Transform translation_start_point(Vecd(0.5 * LL, 8.0, 0.5 * LW));
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
        add<GeometricShapeBox>(halfsize_water);
        //subtract<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
        //subtract<TransformShape<GeometricShapeBox>>(Transform(translation_start_point), halfsize_cubic);
        subtract<GeometricShapeBox>(halfsize_cubic);
        //subtract<GeometricShapeBall>(ball_center, ball_radius);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_wall(halfsize_inner);
        Transform translation_wall_original_axis(Vecd(0.0, 10.0, 0.0));
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall_original_axis), halfsize_inner);
    }
};

//	define the cubic wall boundary shape
class Cubic : public ComplexShape
{
  public:
    explicit Cubic(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_cubic(0.5 * CL, 0.5 * CH, 0.5 * CW);
        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_start_point(Vecd(0.5 * LL, 8.0, 0.5 * LW));
        Transform translation_water(halfsize_water);
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_cubic);
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_start_point), halfsize_cubic);
        add<GeometricShapeBox>(halfsize_cubic);
        //add<GeometricShapeBall>(ball_center, ball_radius);
    }
};


//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        // add observation points
        positions_.push_back(Vecd(DL, 0.01, 0.5 * DW));
        positions_.push_back(Vecd(DL, 0.1, 0.5 * DW));
        positions_.push_back(Vecd(DL, 0.2, 0.5 * DW));
        positions_.push_back(Vecd(DL, 0.24, 0.5 * DW));
        positions_.push_back(Vecd(DL, 0.252, 0.5 * DW));
        positions_.push_back(Vecd(DL, 0.266, 0.5 * DW));
    }
};

class HorizontalMovement: public BaseTracingMethod
{
 public:
     HorizontalMovement() {};
     virtual ~HorizontalMovement(){};

     virtual Vecd tracingPosition (Vecd previous_position, Real current_time = 0.0) override
     {
         Real run_time = GlobalStaticVariables::physical_time_;
         Vecd current_position (0.0, 0.0, 0.0);
         current_position[1]= previous_position[1] - 0.3 * run_time;
         //current_position[1] = previous_position[1];
         /*if(run_time <= 10.0)
         {
             current_position[1] = previous_position[1] + 0.05 * run_time;
         }
         else{
             current_position[1] = previous_position[1] - 0.05 * (run_time-20.0);
         }
         */
         current_position[0] = previous_position[0];
         current_position[2] = previous_position[2];

         return current_position;
     }
};


class RotationMovement : public BaseTracingMethod
{
public:
    RotationMovement(Vecd axis_point_1, Vecd axis_point_2, Real rotation_velocity) : axis_point_1_(axis_point_1), 
        axis_point_2_(axis_point_2), rotation_v_(rotation_velocity)
    {
        rotation_axis_[0] = axis_point_1_[0] - axis_point_2_[0];
        rotation_axis_[1] = axis_point_1_[1] - axis_point_2_[1];
        rotation_axis_[2] = axis_point_1_[2] - axis_point_2_[2];
    };
    virtual ~RotationMovement() {};
   
    virtual Vecd tracingPosition(Vecd previous_position, Real current_time = 0.0) override
    {
        Real rho = (previous_position - findProjectionPoint(previous_position)).norm();
        Real difference_0 = previous_position[0] - findProjectionPoint(previous_position)[0];
        Real difference_1 = previous_position[1] - findProjectionPoint(previous_position)[1];
        Real difference_2 = previous_position[2] - findProjectionPoint(previous_position)[2];
        Real theta = atan2(difference_0, difference_1);
        Real run_time = GlobalStaticVariables::physical_time_;
        Vecd current_position(0.0, 0.0, 0.0);
        current_position[0] = rotation_axis_[0] + cos(theta + rotation_v_ * run_time) * rho;
        current_position[2] = previous_position[2];
        current_position[1] = rotation_axis_[1] + sin(theta + rotation_v_ * run_time) * rho;
       
        return current_position;
    }

    virtual Vecd updateNormalForVector(Vecd previous_vector) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real magnitude = (previous_vector - findProjectionPoint(previous_vector)).norm();
        Real theta = atan2(previous_vector[0], previous_vector[2]);
        Vecd current_vector(0.0, 0.0, 0.0);
        current_vector[0] = magnitude * cos(theta + run_time * rotation_v_);
        current_vector[1] = magnitude * sin(theta + run_time * rotation_v_);
        current_vector[2] = previous_vector[2];

        return current_vector;
    }

protected:
    Vecd findProjectionPoint(Vecd previous_position)
    {
        Vecd point_vector (0.0,0.0,0.0);
        Vecd projection_vector = Vecd::Zero();
        point_vector[0] = previous_position[0] - axis_point_1_[0];
        point_vector[1] = previous_position[1] - axis_point_1_[1];
        point_vector[2] = previous_position[2] - axis_point_1_[2];

        Real axisLengthSquared = projection_vector[0] * rotation_axis_[0] + rotation_axis_[1] * rotation_axis_[1] + rotation_axis_[2] * rotation_axis_[2];
        Real dotProduct = rotation_axis_[0] * point_vector[0] + rotation_axis_[1] * point_vector[1] + rotation_axis_[2] * point_vector[2];

        projection_vector[0] = rotation_axis_[0] * (dotProduct / axisLengthSquared);
        projection_vector[1] = rotation_axis_[1] * (dotProduct / axisLengthSquared);
        projection_vector[2] = rotation_axis_[2] * (dotProduct / axisLengthSquared);

        return projection_vector;
    }
    Vecd rotation_axis_; //two points axis
    Real rotation_v_;
    Vecd axis_point_1_;
    Vecd axis_point_2_;
};

class ThreeDRotation : public BaseTracingMethod
{
public:
    ThreeDRotation(Vecd axis_point_1, Vecd axis_point_2, Real rotation_velocity) :
        axis_point_A_(axis_point_1), angular_v_(rotation_velocity)
    {
        axis_ = (axis_point_1 - axis_point_2).normalized();
        
    }

    virtual Vecd tracingPosition(Vecd previous_position, Real current_time = 0.0) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        Eigen::Quaterniond rotation_position;
        rotation_position = Eigen::AngleAxisd(angular_v_ * run_time_, axis_);
        Eigen::Quaterniond point_quaternion_previos(0, previous_position.x() - axis_point_A_.x(), previous_position.y() - axis_point_A_.y(), previous_position.z() - axis_point_A_.z());
        Eigen::Quaterniond rotated_point_quaternion = rotation_position * point_quaternion_previos * rotation_position.conjugate();
        Vecd new_position(rotated_point_quaternion.x() + axis_point_A_.x(),
            rotated_point_quaternion.y() + axis_point_A_.y(),
            rotated_point_quaternion.z() + axis_point_A_.z());
        return new_position;
    }

    virtual Vecd updateNormalForVector(Vecd previous_vector) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        Eigen::Quaterniond rotation_vector;
        rotation_vector = Eigen::AngleAxisd(angular_v_ * run_time_, axis_);
        Eigen::Quaterniond vector_quaternion_previous(0, previous_vector.x(), previous_vector.y(), previous_vector.z());
        Eigen::Quaterniond rotated_quaternion = rotation_vector * vector_quaternion_previous * rotation_vector.conjugate();
        Vecd new_vector(rotated_quaternion.x(), rotated_quaternion.y(), rotated_quaternion.z());
        return new_vector;
    }
protected:
    Vecd axis_;
    Vecd axis_point_A_;
    Real angular_v_;
    
    
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-50.0, -50.0, -50.0), Vecd(50.0, 50.0, 50.0));
    SPHSystem system(system_domain_bounds, resolution_ref);
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    //water_block.defineBodyLevelSetShape();
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    /*SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody cubic_boundary(system, makeShared<Cubic>("Cubic"));
    cubic_boundary.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    cubic_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    cubic_boundary.generateParticles<ParticleGeneratorLattice>();
    cubic_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");*/

    ObserverBody fluid_observer(system, "FluidObserver");
    fluid_observer.generateParticles<WaterObserverParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    //ComplexRelation water_block_complex(water_block, RealBodyVector{&wall_boundary, &cubic_boundary});
    InnerRelation water_block_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> density_relaxation(water_block_inner);
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSizeMovingVelocity> get_fluid_advection_time_step_size(water_block, U_f, angular_velocity);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSizeMovingVelocity> get_fluid_time_step_size(water_block, angular_velocity, 0.2);
    //SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_block_inner);

    /** Define the confinement condition for wall. */
    NearShapeSurface near_surface_wall(water_block, makeShared<WallBoundary>("Wall"));
    near_surface_wall.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinementGeneral confinement_condition_wall(near_surface_wall);

    //RotationMovement circle_movement(axis_point_1, axis_point_2, 0.2*Pi);
    ThreeDRotation circle_movement(axis_point_1, axis_point_2, angular_velocity);
    //HorizontalMovement horizaontal_movement;
    NearShapeSurfaceTracing near_surface_cubic(water_block, makeShared<InverseShape<Cubic>>("Cubic"), circle_movement);
    near_surface_cubic.level_set_shape_.writeLevelSet(io_environment);
    //fluid_dynamics::StaticConfinement confinement_condition_cubic(near_surface_cubic);
    fluid_dynamics::MovingConfinementGeneral confinement_condition_cubic(near_surface_cubic);

    /** Define the confinement condition for structure. */
    update_density_by_summation.post_processes_.push_back(&confinement_condition_wall.density_summation_);
    update_density_by_summation.post_processes_.push_back(&confinement_condition_cubic.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_wall.pressure_relaxation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_cubic.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_cubic.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.surface_bounding_);
    density_relaxation.post_processes_.push_back(&confinement_condition_cubic.surface_bounding_);
    viscous_acceleration.post_processes_.push_back(&confinement_condition_wall.viscous_acceleration_);
    viscous_acceleration.post_processes_.push_back(&confinement_condition_cubic.viscous_acceleration_);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_water_block_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    //wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 60.0;
    Real output_interval = end_time / 20.0;
    Real dt = 0.0; // default acoustic time step sizes
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_water_block_states.writeToFile(0);
    write_water_mechanical_energy.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {

                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                //write_water_block_states.writeToFile();
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                write_water_block_states.writeToFile();
            }
            number_of_iterations++;

            water_block.updateCellLinkedListWithParticleSort(100);
            //water_block_complex.updateConfiguration();
            water_block_inner.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            near_surface_cubic.updateCellList();
            write_recorded_water_pressure.writeToFile(number_of_iterations);
            
        }

        write_water_mechanical_energy.writeToFile(number_of_iterations);

        TickCount t2 = TickCount::now();
        //write_water_block_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.generate_regression_data_)
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        write_recorded_water_pressure.generateDataBase(1.0e-3);
    }
    else
    {
        write_water_mechanical_energy.testResult();
        write_recorded_water_pressure.testResult();
    }

    return 0;
}
