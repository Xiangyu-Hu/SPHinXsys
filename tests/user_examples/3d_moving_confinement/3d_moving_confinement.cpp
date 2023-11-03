/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
#include "body_part_by_cell_tracing.h"
#include "level_set_confinement.h"
#include "math.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the first geometry.
//	To use this, please commenting the setting for the second geometry.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/propeller.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 400.0; // Domain width.
Real DH = 400.0;  // Domain height.
Real DW = 400.0;  // Domain width.
Real DH_T = 600.0; 
Vec3d domain_lower_bound(-300, -300.0, -300.0);
Vec3d domain_upper_bound(300, 500, 300);
Vecd translation(0.0, 0.0, 0.0);
Real scaling = 1.0;

//Vecd insert_circle_center (2.0, 1.0);
Real insert_circle_radius = 0.25;

Vecd rotation_axis_1(0.0, 0.0, 0.0);
Vecd rotation_axis_2(0.0, 100.0, 0.0);

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity force of fluid. */
Real U_max = 2.0 * sqrt(gravity_g * DH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */

//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
Real BW = dp_0 * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	define the imported model.
//----------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
    }
};
/** create a water block shape */
//	define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_water(Vecd::Zero());
        Vecd halfsize_cubic(0.1 * DL, 0.1 * DH, 0.1 * DW);
        Transform translation_cubic(Vecd::Zero());
        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
        //subtract<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_cubic), halfsize_cubic);
    }
};

class Cubic : public ComplexShape
{
  public:
    explicit Cubic(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_cubic(0.1 * DL, 0.1 * DH, 0.1 * DW);
        Transform translation_cubic(Vecd::Zero());
        add<TransformShape<GeometricShapeBox>>(Transform(translation_cubic), halfsize_cubic);
        
    }
};

/** create wall shape */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_out(0.5 * DL + BW, 0.5 * DH_T + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH_T, 0.5 * DW);
        Vecd translation_point(0.0, 100.0, 0.0);
        Transform translation_out(translation_point);
        Transform translation_inner(translation_point);
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_out), halfsize_out);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_inner), halfsize_inner);
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
         current_position[0]= previous_position[0] - 0.1 * run_time;
         current_position[1] = previous_position[1];
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
        Real theta = atan2(difference_0, difference_2);
        Real run_time = GlobalStaticVariables::physical_time_;
        Vecd current_position(0.0, 0.0, 0.0);
        current_position[0] = rotation_axis_[0] + cos(theta + rotation_v_ * run_time) * rho;
        current_position[1] = previous_position[1];
        current_position[2] = rotation_axis_[1] + sin(theta + rotation_v_ * run_time) * rho;
       
        return current_position;
    }

    virtual Vecd updateNormalForVector(Vecd previous_vector) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real magnitude = (previous_vector - findProjectionPoint(previous_vector)).norm();
        Real theta = atan2(previous_vector[0], previous_vector[2]);
        Vecd current_vector(0.0, 0.0, 0.0);
        current_vector[0] = magnitude * cos(theta + run_time * rotation_v_);
        current_vector[2] = magnitude * sin(theta + run_time * rotation_v_);
        current_vector[1] = previous_vector[1];

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

//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, dp_0);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody wall(system, makeShared<WallBoundary>("Wall"));
   // wall.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    wall.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall.generateParticles<ParticleGeneratorLattice>();
    wall.addBodyStateForRecording<Vecd>("NormalDirection");
    
    //RealBody imported_model(system, makeShared<Cubic>("cubic"));
    //// level set shape is used for particle relaxation
    //imported_model.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    //imported_model.defineParticlesAndMaterial();
    //imported_model.generateParticles<ParticleGeneratorLattice>();
  
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    BodyStatesRecordingToVtp body_states_recording(io_environment, system.real_bodies_);
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    //ComplexRelation water_block_complex(water_block, {&wall });
    /** Initialize particle acceleration. */
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> density_relaxation(water_block_inner);
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);

    //RotationMovement rotation_movement(rotation_axis_1, rotation_axis_2, 0.5 * Pi);
    ThreeDRotation rotation_movement(rotation_axis_1, rotation_axis_2, 0.5 * Pi);
    //CircleMovement circle_movement(square_center, Pi);
    //HorizontalMovement horizaontal_movement;
    NearShapeSurfaceTracing near_surface_circle(water_block, makeShared<InverseShape<Cubic>>("cubic"), rotation_movement);
    near_surface_circle.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::MovingConfinementGeneral confinement_condition_circle(near_surface_circle);

    /** Define the confinement condition for wall. */
    NearShapeSurface near_surface_wall(water_block, makeShared<WallBoundary>("Wall"));
    near_surface_wall.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinement confinement_condition_wall(near_surface_wall);
    /** Define the confinement condition for structure. */

   
    
    /*NearShapeSurface near_surface_triangle(water_block, makeShared<InverseShape<Triangle>>("Triangle"));
    near_surface_triangle.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinement confinement_condition_triangle(near_surface_triangle);*/
    /** Push back the static confinement conditiont to corresponding dynamics. */
    update_density_by_summation.post_processes_.push_back(&confinement_condition_wall.density_summation_);
    update_density_by_summation.post_processes_.push_back(&confinement_condition_circle.density_relaxation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_wall.pressure_relaxation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_circle.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_circle.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.surface_bounding_);
    density_relaxation.post_processes_.push_back(&confinement_condition_circle.surface_bounding_);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    //BodyStatesRecordingToVtp body_states_recording(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
    /*RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);*/
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 10.0;       /**< End time. */
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;              /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    //body_states_recording.writeToFile(0);
    write_water_mechanical_energy.writeToFile(0);
    //write_recorded_water_pressure.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force and gravity. */
            time_instance = TickCount::now();
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                //body_states_recording.writeToFile();


            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    //write_recorded_water_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_inner.updateConfiguration();
            //water_block_complex.updateConfiguration();
            //fluid_observer_contact.updateConfiguration();
            near_surface_circle.updateCellList();
            interval_updating_configuration += TickCount::now() - time_instance;
            
        }

        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}
