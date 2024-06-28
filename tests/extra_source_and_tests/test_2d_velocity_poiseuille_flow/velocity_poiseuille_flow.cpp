/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
#include "base_particle_dynamics.h"
#include "bidirectional_buffer.h"
#include "data_type.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "io_environment.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
#include <cmath>
#include <fstream> // Include fstream to handle file writing
#include <gtest/gtest.h>
#include <string>

using namespace SPH;
Real DH = 0.;             // prescribed value in main()
Real U_f = 0.;            // prescribed value in main()
Real outlet_pressure = 0; // prescribed value in main()

//----------------------------------------------------------------------
//	Circular buffer for checking convergence usage
//----------------------------------------------------------------------
template <typename T>
class CircularBuffer
{
  private:
    std::vector<T> buffer;
    size_t head = 0;
    size_t capacity;
    size_t count = 0;

  public:
    CircularBuffer(size_t size) : capacity(size)
    {
        buffer.resize(size);
    }

    void push(const T &item)
    {
        buffer[head] = item;
        head = (head + 1) % capacity;
        if (count < capacity)
            count++;
    }

    T &operator[](size_t index)
    {
        return buffer[(head - count + index) % capacity];
    }

    size_t size() const
    {
        return count;
    }

    bool is_full() const
    {
        return count == capacity;
    }
};

//----------------------------------------------------------------------
//	ConvergenceChecker
//----------------------------------------------------------------------
template <typename T>
class ConvergenceChecker
{
  private:
    CircularBuffer<T> buffer;
    T threshold;
    T percentage_difference_ = std::numeric_limits<T>::infinity();

    bool calculate_convergence()
    {
        if (!buffer.is_full())
        {
            return false; // Not enough data to determine convergence
        }

        T max_val = buffer[0];
        T min_val = buffer[0];
        for (size_t i = 1; i < buffer.size(); ++i)
        {
            if (buffer[i] > max_val)
                max_val = buffer[i];
            if (buffer[i] < min_val)
                min_val = buffer[i];
        }

        // Calculate the percentage difference
        T range = max_val - min_val;
        T average = (max_val + min_val) / 2;
        T percentage_difference = (range / average) * 100;
        percentage_difference_ = percentage_difference;
        std::cout << "converger percentage_difference_ :" << percentage_difference_ << "\n";
        return percentage_difference < threshold;
    }

  public:
    ConvergenceChecker(size_t size, T conv_threshold)
        : buffer(size), threshold(conv_threshold) {}

    bool update(T new_value)
    {
        buffer.push(new_value);
        return calculate_convergence();
    }

    T get_percentage_difference()
    {
        return percentage_difference_;
    }
};

//----------------------------------------------------------------------
//	Function to create observer locations
//----------------------------------------------------------------------
std::vector<Vecd> createObserverLocations(double x, double DH, int num_points)
{
    std::vector<Vecd> observer_location;
    observer_location.reserve(num_points);
    double dy = DH / (num_points - 1);

    for (int i = 0; i < num_points; ++i)
    {
        double y = -0.5 * DH + i * dy;
        observer_location.push_back(Vecd(x, y));
    }

    return observer_location;
}

//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

struct RightInflowPressure
{
    Real outlet_pressure_;

    RightInflowPressure(Real outlet_pressure) {}

    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) : outlet_pressure_(outlet_pressure) {}

    Real operator()(Real &p_)
    {
        return outlet_pressure_;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
Real parabolic_velocity(Real y, Real U_f, Real DH)
{
    return U_f * (pow(0.5 * DH, 2) - pow(y, 2)) / pow(0.5 * DH, 2);
}

struct InflowVelocity
{
    Real u_ref_, t_ref_, DH_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(5.0), DH_(DH)
    {
    }

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;

        target_velocity[0] = (1.0 - pow(position[1] / (0.5 * DH_), 2)) * u_ave;

        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name, const Real DH, const Real DL) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, -0.5 * DH));
        water_block_shape.push_back(Vecd(0.0, 0.5 * DH));
        water_block_shape.push_back(Vecd(DL, 0.5 * DH));
        water_block_shape.push_back(Vecd(DL, -0.5 * DH));
        water_block_shape.push_back(Vecd(0.0, -0.5 * DH));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        std::cout << "water done" << std::endl;
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name, const Real DH, const Real DL, const Real BW) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(0.0, -0.5 * DH - BW));
        outer_wall_shape.push_back(Vecd(0.0, 0.5 * DH + BW));
        outer_wall_shape.push_back(Vecd(DL, 0.5 * DH + BW));
        outer_wall_shape.push_back(Vecd(DL, -0.5 * DH - BW));
        outer_wall_shape.push_back(Vecd(0.0, -0.5 * DH - BW));
        std::cout << "wall outer done" << std::endl;
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(-BW, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(-BW, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(DL + BW, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(DL + BW, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(-BW, -0.5 * DH));
        std::cout << "wall inner done" << std::endl;

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        std::cout << "wall outer done" << std::endl;
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
        std::cout << "wall inner done" << std::endl;
    }
};
//----------------------------------------------------------------------
//	Basic geometry parameters for FDA.
//----------------------------------------------------------------------
//        ______(L_in)____    |    (L_throat)    |________(L_out)________
//       |                \                      |                       |
//       |                 \                     |                       |
//       |                  \                    |                       |
//       |                   \___________________|                       |
//       |                            |                                  |
//       |                            |                                  |
//      (DH)                       H_noozle                              |
//       |                            |                                  |
//       |           _________________|__________                        |
//       |             alpha /                   |                       |
//       |                \ /                  (DH/3)                    |
//       |                 /                     |                       |
//       |________________/                      |_______________________|
//       |                |   |
//       |              (L_slope)
//       |
//----------------------------------------------------------------------
//	Nozzle fluid body definition.
//----------------------------------------------------------------------
class NozzleWaterBlock : public MultiPolygonShape
{
  public:
    explicit NozzleWaterBlock(const std::string &shape_name, const Real DH, const Real L_in, const Real L_throat, const Real L_slope, const Real L_out) : MultiPolygonShape(shape_name)
    {
        const Real H_noozle = DH / 3;
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, -0.5 * DH));
        water_block_shape.push_back(Vecd(0.0, 0.5 * DH));
        water_block_shape.push_back(Vecd(L_in, 0.5 * DH));
        water_block_shape.push_back(Vecd(L_in + L_slope, 0.5 * DH - H_noozle));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH - H_noozle));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out, 0.5 * DH));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out, -0.5 * DH));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH));
        water_block_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH + H_noozle));
        water_block_shape.push_back(Vecd(L_in + L_slope, -0.5 * DH + H_noozle));
        water_block_shape.push_back(Vecd(L_in, -0.5 * DH));
        water_block_shape.push_back(Vecd(0.0, -0.5 * DH));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Nozzle wall boundary body definition.
//----------------------------------------------------------------------
class NozzleWallBoundary : public MultiPolygonShape
{
  public:
    explicit NozzleWallBoundary(const std::string &shape_name, const Real DH, const Real L_in, const Real L_throat, const Real L_slope, const Real L_out, const Real Bw) : MultiPolygonShape(shape_name)
    {
        const Real H_noozle = DH / 3;
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(0.0, -0.5 * DH - Bw));
        outer_wall_shape.push_back(Vecd(0.0, 0.5 * DH + Bw));
        outer_wall_shape.push_back(Vecd(L_in, 0.5 * DH + Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope, 0.5 * DH - H_noozle + Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH - H_noozle + Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH + Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out, 0.5 * DH + Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out, -0.5 * DH - Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH - Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH + H_noozle - Bw));
        outer_wall_shape.push_back(Vecd(L_in + L_slope, -0.5 * DH + H_noozle - Bw));
        outer_wall_shape.push_back(Vecd(L_in, -0.5 * DH - Bw));
        outer_wall_shape.push_back(Vecd(0.0, -0.5 * DH - Bw));
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(0.0 - Bw, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(0.0 - Bw, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in + L_slope, 0.5 * DH - H_noozle));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH - H_noozle));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out + Bw, 0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat + L_out + Bw, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(L_in + L_slope + L_throat, -0.5 * DH + H_noozle));
        inner_wall_shape.push_back(Vecd(L_in + L_slope, -0.5 * DH + H_noozle));
        inner_wall_shape.push_back(Vecd(L_in, -0.5 * DH));
        inner_wall_shape.push_back(Vecd(0.0 - Bw, -0.5 * DH));
        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
void channel_flow(int ac, char *av[], const Real length_to_height_ratio, const size_t number_of_particles, const bool is_FDA = false, const bool use_transport_correction = true, const bool use_linear_gradient_correction = true, const size_t number_of_observer = 50)
{
    //----------------------------------------------------------------------
    //	Basic geometry parameters and numerical setup.
    //----------------------------------------------------------------------
    DH = 2.0;                              /**< Channel height. */
    Real DL = DH * length_to_height_ratio; /**< Channel length. */
    const Real L_throat = DH * 4 / 1.2;    // Based on FDA geometry
    const Real L_slope = DH * 2.2685 / 1.2;
    const Real L_out = DH * 7;
    const Real L_in = DL - L_throat - L_slope - L_out;
    //----------------------------------------------------------------------
    //	Global parameters on the fluid properties
    //----------------------------------------------------------------------
    const Real rho0_f = 1.0;                  /**< Density. */
    U_f = 1.0;                                /**< Characteristic velocity. */
    const Real Re = 100.0;                    /**< Reynolds number. */
    const Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
    Real predicted_pressure = 8 * mu_f * DL * U_f / pow(DH, 2);
    Real maximum_pressure_fluctuation = predicted_pressure * 2.0;
    Real c_f = std::max(sqrt(maximum_pressure_fluctuation / 0.01 * (1 / rho0_f)), U_f * 10.); // sqrt(maximum_pressure_fluctuation / 0.01 * (1 / rho0_f)) means maximum_pressure_fluctuation will be 1% of p0
    std::cout << "predicted_pressure = " << predicted_pressure << ", accepted predicted_pressure fluctuation= " << maximum_pressure_fluctuation << std::endl;
    Real U_max = std::max(c_f / 10.0, U_f);
    outlet_pressure = 0.0;
    //----------------------------------------------------------------------
    Real resolution_ref = DH / number_of_particles; /**< Initial reference particle spacing. */
    Real BW = resolution_ref * 4.;                  /**< Extending width for BCs. */
    BoundingBox b_box(Vec2d(-BW, -BW - 0.5 * DH), Vec2d(DL + BW, 0.5 * DH + BW));
    BoundingBox system_domain_bounds(b_box);
    //----------------------------------------------------------------------
    //	Geometric shapes used in this case.
    //----------------------------------------------------------------------
    Real boundary_width = 3.0 * resolution_ref;
    Real half_boundary_width = 0.5 * boundary_width;
    Vec2d bidirectional_buffer_halfsize = Vec2d(half_boundary_width, DH * 0.5);
    Vec2d left_bidirectional_translation = Vec2d(half_boundary_width, 0.);
    Vec2d right_bidirectional_translation = Vec2d(DL - half_boundary_width, 0.);
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    // GlobalStaticVariables::physical_time_ is reset to zero to ensure clean starts for cases
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    std::string filename_suffix = "_DL" + std::to_string(int(length_to_height_ratio));
    filename_suffix += "_NumPart" + std::to_string(int(number_of_particles));
    filename_suffix += "_TransPort" + std::to_string(int(use_transport_correction));
    filename_suffix += "LGC" + std::to_string(int(use_linear_gradient_correction));
    if (is_FDA)
        filename_suffix += "FDA";
    auto &output_folder = sph_system.getIOEnvironment().output_folder_;
    output_folder += filename_suffix;
    std::cout << "sph_system.getIOEnvironment().output_folder_ : " << sph_system.getIOEnvironment().output_folder_ << std::endl;
    std::filesystem::remove_all(sph_system.getIOEnvironment().output_folder_);
    std::filesystem::create_directory(sph_system.getIOEnvironment().output_folder_);

    GlobalStaticVariables::physical_time_ = 0.;
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    std::shared_ptr<MultiPolygonShape> water_block_shape;
    std::shared_ptr<MultiPolygonShape> wall_block_shape;

    if (is_FDA)
    {
        water_block_shape = makeShared<NozzleWaterBlock>("NozzleWaterBody", DH, L_in, L_throat, L_slope, L_out);
        wall_block_shape = makeShared<NozzleWallBoundary>("NozzleWallBody", DH, L_in, L_throat, L_slope, L_out, BW);
    }
    else
    {
        water_block_shape = makeShared<WaterBlock>("WaterBody", DH, DL);
        wall_block_shape = makeShared<WallBoundary>("WallBoundary", DH, DL, BW);
    }

    FluidBody water_block(sph_system, water_block_shape);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_block_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    Real radial_observer_x = DL * 0.5;
    if (is_FDA)
        radial_observer_x = L_in * 0.5;
    velocity_observer.generateParticles<BaseParticles, Observer>(createObserverLocations(radial_observer_x, DH, 50));

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block, &wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxillary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell right_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_emitter_inflow_injection(right_emitter, in_outlet_particle_buffer, xAxis);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);

    Vecd inflow_translation(half_boundary_width, 0);
    Vecd inflow_halfsize(half_boundary_width, 0.5 * DH);
    // BodyAlignedBoxByCell inflow_region(
    //     water_block, makeShared<AlignedBoxShape>(Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    BodyAlignedBoxByCell inflow_region(
        water_block, makeShared<AlignedBoxShape>(Transform(inflow_translation), inflow_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inflow_region);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp waterblock_recording(sph_system);
    BodyStatesRecordingToVtp observer_recording(velocity_observer);

    waterblock_recording.addVariableRecording<Real>(water_block, "Pressure");
    waterblock_recording.addVariableRecording<int>(water_block, "Indicator");
    waterblock_recording.addVariableRecording<Real>(water_block, "Density");
    waterblock_recording.addVariableRecording<int>(water_block, "BufferParticleIndicator");
    ObservingAQuantity<Vecd> update_observer_velocity(velocity_observer_contact, "Velocity");
    observer_recording.addVariableRecording<Vecd>(velocity_observer, "Velocity");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 20.0;       /**< End time. */
    Real max_end_time = 1000.0; /**< End time. */
    Real Output_Time = 2;       /**< Time stamps for output of body states. */
    Real dt = 0.0;              /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Defined convergence checker
    //----------------------------------------------------------------------
    ConvergenceChecker<double> conv_checker(10, 0.5); // Buffer of 10 values, convergence threshold of 0.5 percent
    bool is_converged = false;
    size_t mid_index_of_observer = number_of_observer / 2.0;
    StdLargeVec<Vecd> &pos_radial = velocity_observer.getBaseParticles().ParticlePositions();
    StdLargeVec<Vecd> &vel_radial = *velocity_observer.getBaseParticles().getVariableByName<Vecd>("Velocity");
    auto &vel_of_mid_index_observer = vel_radial[mid_index_of_observer][0];
    auto &pos_y_of_mid_index_observer = pos_radial[mid_index_of_observer][1];

    int convergence_checker_output_interval = 50;
    std::ofstream file_mid_observer(sph_system.getIOEnvironment().output_folder_ + "/output_velocity_of_mid_observer_" + std::to_string(int(length_to_height_ratio)) + "DH_" + std::to_string(number_of_particles) + ".csv");
    file_mid_observer << "Time,Velocity X,Convergence Rate\n";

    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    waterblock_recording.writeToFile();
    observer_recording.writeToFile();
    update_observer_velocity.exec();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while ((GlobalStaticVariables::physical_time_ < end_time || !is_converged) && GlobalStaticVariables::physical_time_ <= max_end_time)
    {

        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            if (use_linear_gradient_correction)
                corrected_configuration_fluid.exec();
            viscous_acceleration.exec();
            if (use_transport_correction)
                transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "	convergence = " << conv_checker.get_percentage_difference() << "    is_converged = " << is_converged << "\n";
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();
            if (number_of_iterations % convergence_checker_output_interval == 0)
            {
                velocity_observer_contact.updateConfiguration();
                update_observer_velocity.exec();
                is_converged = conv_checker.update(vel_of_mid_index_observer);
                std::cout << "add value to analysis convergence :" << vel_of_mid_index_observer << ", analytical value" << parabolic_velocity(pos_y_of_mid_index_observer, U_f, DH) << std::endl;
                if (is_converged)
                {
                    std::cout << "Converged at iteration " << vel_of_mid_index_observer << std::endl;
                }
                if (std::isfinite(conv_checker.get_percentage_difference()))
                    file_mid_observer << GlobalStaticVariables::physical_time_ << "," << vel_of_mid_index_observer << "," << conv_checker.get_percentage_difference() << "\n";
                else
                    file_mid_observer << GlobalStaticVariables::physical_time_ << "," << vel_of_mid_index_observer << "," << 0 << "\n";
            }
        }
        TickCount t2 = TickCount::now();
        velocity_observer_contact.updateConfiguration();
        update_observer_velocity.exec();
        waterblock_recording.writeToFile();
        observer_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    file_mid_observer.close();

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

    // Create and open a CSV file
    std::ofstream file_radial_velocity(sph_system.getIOEnvironment().output_folder_ + "/output_velocity_" + std::to_string(int(length_to_height_ratio)) + "DH_" + std::to_string(number_of_particles) + ".csv");
    // Write the header row
    file_radial_velocity << "Position Y,Parabolic Velocity X,Velocity X\n";

    // Loop over all particles to calculate and write the required data
    for (size_t i = 0; i < pos_radial.size(); i++)
    {
        // Fetch parabolic velocity for the current position
        Real par_vel = parabolic_velocity(pos_radial[i][1], U_f, DH);

        // Write data to file
        file_radial_velocity << pos_radial[i][1] << "," << par_vel << "," << vel_radial[i][0] << "\n";

        // Existing test, assuming you still want to run it
        EXPECT_NEAR(parabolic_velocity(pos_radial[i][1], U_f, DH), vel_radial[i][0], U_f * 5e-2);
    }

    // Close the file
    file_radial_velocity.close();
}
int main(int ac, char *av[])
{
    channel_flow(ac, av, 25, 10, true);
    return 0;
}
