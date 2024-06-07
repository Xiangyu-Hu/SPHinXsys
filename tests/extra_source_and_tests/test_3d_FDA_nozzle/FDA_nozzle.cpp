//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include <cmath>
#include <functional>
#include <gtest/gtest.h>
#include <ostream>
#include <sphinxsys.h>

#include <numbers>
using namespace SPH;

constexpr Real Pi = Real(M_PI);

struct AxialVelocityProfile
{
    double z;
    std::vector<std::pair<double, double>> velocity_data; // Pair of position and velocity

    // Default constructor
    AxialVelocityProfile() : z(0)
    {
    }

    // Constructor with initial z-value
    AxialVelocityProfile(double z) : z(z)
    {
    }

    void addData(double position, double velocity) { velocity_data.emplace_back(position, velocity); }
};

void readDataFromFile(
    const std::string &filename,
    std::string string_pattern,
    std::map<double, AxialVelocityProfile> &velocity_profiles)
{
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Failed to open file" << std::endl;
        return;
    }

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;

        // Look for the specific pattern indicating the start of the data we care about
        if (line.find(string_pattern) != std::string::npos)
        {
            double z_value;
            int count;

            // Read and discard the first token ("plot-profile-axial-velocity-at-z")
            iss >> token;   // Consumes "plot-profile-axial-velocity-at-z"
            iss >> z_value; // Next token should be the z value

            // Debug output
            std::cout << "Found profile for Z-value: " << z_value << std::endl;

            // Move to the next line to read the count of data points
            if (std::getline(file, line))
            {
                std::istringstream iss_count(line);
                iss_count >> count; // First integer on the next line should be the count

                std::cout << "Data count: " << count << std::endl;

                AxialVelocityProfile profile(z_value);

                // Read the actual data points
                for (int i = 0; i < count; ++i)
                {
                    double position, velocity;
                    if (file >> position >> velocity)
                    {
                        profile.addData(position, velocity);
                    }
                }

                // Use emplace to avoid unnecessary object creation
                auto result = velocity_profiles.emplace(z_value, std::move(profile));
                if (!result.second)
                {
                    std::cerr << "Profile for z=" << z_value << " already exists." << std::endl;
                }
            }
            else
            {
                std::cerr << "Expected data count line after z-value but got nothing." << std::endl;
            }
        }
    }
};

StdVec<Vecd> ObersverAxialGenerator(double x_min, double x_max)
{
    StdVec<Vecd> observer_location;
    int ny = 51;
    double full_length = x_max - x_min;
    for (int i = 0; i < ny; ++i)
    {
        double x = full_length / (ny - 1) * i + x_min;
        observer_location.emplace_back(Vec3d(x, 0, 0));
    }
    return observer_location;
}

StdVec<Vecd> ObserverRadialGenerator(double x, double diameter)
{
    StdVec<Vecd> observer_location;

    int ny = 51;
    for (int i = 0; i < ny - 1; i++)
    {
        double z = diameter * i / (ny - 1);
        observer_location.emplace_back(Vec3d(x, 0, -diameter / 2.0 + z));
    }
    return observer_location;
}

double compute_pressure(double p)
{
    double run_time = GlobalStaticVariables::physical_time_;
    double pressure = run_time < 0.5 ? 0.5 * p * (1.0 - cos(M_PI * run_time / 0.5)) : p;

    // double pressure = p;
    return pressure;
}

struct FDA_nozzle_parameters
{ // using default blood parameters
    // resolution, particles per diameter
    double scale = 0.001;
    double inlet_diameter = 12 * scale; // Entering part diameter
    double thoat_diameter = 4 * scale;
    double outlet_diameter = 12 * scale;

    uint number_of_particles = 10;

    Vec3d inlet_normal = Vec3d::UnitX();
    Vec3d inlet_center = Vec3d(-120 * scale, 0, 0);

    Vec3d outlet_normal = Vec3d::UnitX() * -1;
    Vec3d outlet_center = Vec3d(120 * scale, 0, 0);

    // defualt inlet area
    double inlet_area = pow(inlet_diameter, 2) * M_PI * 0.25;
    double throat_area = pow(thoat_diameter, 2) * M_PI * 0.25;

    double Q_f = 5e-6; // for Re = 500

    // material paramaters, default: blood
    // Newtonian
    double rho0_f = 1056; //[kg/m3]
    double mu_f = 0.0035; //[N s / m2]

    // fluid flle path
    fs::path fluid_file_path = "./input/FDA_nozzle_fluid.stl";
    // wall flle path
    fs::path wall_file_path;
    // total length of fluid
    double fluid_length = 200 * scale; // defined in mesh file
    double outlet_pressure = 0;
    // end time
    double end_time = 15.0;
};

void setup_directory(const std::string &path)
{
    // Check if the directory exists
    if (fs::exists(path))
    {
        // Try to remove the directory and its contents
        if (!fs::remove_all(path))
        {
            std::cerr << "Failed to remove existing directory: " << path << std::endl;
            return;
        }
    }

    // Create the directory
    if (!fs::create_directory(path))
    {
        std::cerr << "Failed to create directory: " << path << std::endl;
        return;
    }

    std::cout << "Directory set up at: " << path << std::endl;
}

void print_profiles_data(const std::map<double, AxialVelocityProfile> &profiles, const std::string &description)
{
    std::cout << description << std::endl;
    for (const auto &entry : profiles)
    {
        double z = entry.first;
        const AxialVelocityProfile &profile = entry.second;

        std::cout << "Profile for Z-value: " << z << std::endl;
        for (const auto &data : profile.velocity_data)
        {
            std::cout << "    Position: " << data.first << ", Value: " << data.second << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
}

const AxialVelocityProfile get_profiles_data(const std::map<double, AxialVelocityProfile> &profiles, double z)
{
    auto it = profiles.find(z); // Use the find method to search for the z value directly
    return (it != profiles.end()) ? it->second
                                  : throw std::runtime_error("Profile not found for z-value: " + std::to_string(z));
}

/**
 * @brief 	Pressure boundary definition.
 */
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
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = 0;
        return pressure;
    }
};

Real U_f = 0.;
Real inlet_diameter = 0;

/**
 * @brief 	inflow velocity definition.
 */
struct InflowVelocity
{
    Real u_max_;
    Real radius_squared_; // radius at inlet

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_max_(0.0), radius_squared_(0.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vec3d pos = position;
        pos[0] = 0; // set x direction to 0
        double radius_square = pos.squaredNorm();

        u_max_ = U_f * 2;
        radius_squared_ = pow(inlet_diameter * 0.5, 2);
        double vel = u_max_ * (1 - radius_square / radius_squared_);
        Vec3d target_velocity = Vec3d::UnitX() * vel;

        return target_velocity;
    }

    void set_u_max(double u_max)
    {
        u_max_ = u_max;
    }

    void set_radius_squared_(double radius_squared)
    {
        radius_squared_ = radius_squared;
    }
};

void FDA_nozzle(int ac, char *av[], FDA_nozzle_parameters &params, const double Re = 500)
{

    const double scale = params.scale;
    inlet_diameter = params.inlet_diameter;
    const double Q_f = params.Q_f;
    U_f = Q_f / params.inlet_area;
    const double U_max = Q_f / params.throat_area * 2; // U max suppose to happend at throat
    const double c_f = 10 * U_max;                     // Speed of sound

    std::cout << "U_f (velocity at inlet): " << U_f << std::endl;
    std::cout << "Re: " << Re << std::endl;

    // Map to hold each profile accessed by the z-value as a key (obtaining radial velocity)
    std::map<double, AxialVelocityProfile> radial_velocity_profiles;
    std::map<double, AxialVelocityProfile> axial_velocity_profiles;
    std::map<double, AxialVelocityProfile> axial_pressure_profiles;

    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-profile-axial-velocity-at-z", radial_velocity_profiles);
    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-z-distribution-axial-velocity", axial_velocity_profiles);
    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-z-distribution-pressure", axial_pressure_profiles);

    // radial observer at x axis
    std::vector<double> radiao_observer_X;
    for (const auto &entry : radial_velocity_profiles)
    {
        double z = entry.first;
        // Print the Z-value
        std::cout << "Profile for Z-value: " << z << std::endl;
        radiao_observer_X.emplace_back(z);
    }

    // Inside your FDA_nozzle function or the main function after data has been loaded
    std::cout << "Printing Axial Velocity Profiles Data:" << std::endl;
    print_profiles_data(axial_velocity_profiles, "Axial Velocity Profiles:");

    std::cout << "Printing Axial Pressure Profiles Data:" << std::endl;
    print_profiles_data(axial_pressure_profiles, "Axial Pressure Profiles:");

    const uint number_of_particles = params.number_of_particles;
    const double resolution_ref = params.thoat_diameter / number_of_particles;
    std::cout << "resolution_ref: " << resolution_ref << std::endl;

    // const double resolution_wall = resolution_ref;
    // const uint simtk_resolution = 10;

    const double rho0_f = params.rho0_f;
    const double mu_f = params.mu_f;

    // GEOMETRY
    auto fluid_shape = makeShared<ComplexShape>("fluid");
    fluid_shape->add<TriangleMeshShapeSTL>(params.fluid_file_path, Vec3d::Zero(), scale);
    const double x_min_domain = fluid_shape->getBounds().first_.x();
    const double x_max_domain = fluid_shape->getBounds().second_.x();

    Real boundary_width = 5 * resolution_ref;
    Real half_boundary_width = boundary_width * 0.5;
    Vecd bidirectional_buffer_halfsize = Vec3d(half_boundary_width, params.inlet_diameter * 0.5, params.inlet_diameter * 0.5);
    Vec3d left_disposer_translation = Vec3d(x_min_domain + half_boundary_width, 0, 0);
    Vec3d right_disposer_translation = Vec3d(x_max_domain - half_boundary_width, 0, 0);

    Vecd left_bidirectional_translation = Vec3d(x_min_domain + half_boundary_width, 0, 0);
    Vecd right_bidirectional_translation = Vec3d(x_max_domain - half_boundary_width, 0, 0);

    auto wall_shape = makeShared<ComplexShape>("wall");
    wall_shape->add<TriangleMeshShapeSTL>(params.wall_file_path, Vec3d::Zero(), scale);

    // SYSTEM
    SPHSystem sph_system(wall_shape->getBounds(), resolution_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    IOEnvironment in_output(sph_system);

    setup_directory(in_output.output_folder_);

    // Check if the folder exists
    std::filesystem::path folderPath = in_output.output_folder_;
    if (std::filesystem::exists(folderPath) && std::filesystem::is_directory(folderPath))
    {
        // Delete the folder if it exists
        std::filesystem::remove_all(folderPath);
    }

    // Create the folder
    std::filesystem::create_directory(folderPath);
    // FLUID
    FluidBody water_block(sph_system, fluid_shape);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<Lattice>();

    InnerRelation wall_boundary_inner(wall_boundary);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    /** check whether run particle relaxation for body fitted particle distribution. */
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        // Write the particle reload files.
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        // A  Physics relaxation step.
        RelaxationStepInner relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation loop.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        // Output particles position for reload.
        write_particle_reload_files.writeToFile(0);
    }

    ObserverBody observer_axial(sph_system, "fluid_observer_axial");
    observer_axial.generateParticles<Observer>(ObersverAxialGenerator(x_min_domain, x_max_domain));

    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation axial_velocity_observer_contact(observer_axial, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** delete outflow particles */
    BodyAlignedBoxByCell left_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(left_disposer_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vecd(right_disposer_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        boundary_indicator(water_block_inner, water_block_contact);
    /** bidirectional buffer */
    BodyAlignedBoxByCell left_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize));

    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure, SequencedPolicy> left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell right_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(right_bidirectional_translation)), bidirectional_buffer_halfsize));

    fluid_dynamics::BidirectionalBuffer<RightInflowPressure, SequencedPolicy> right_emitter_inflow_injection(right_emitter, in_outlet_particle_buffer, xAxis); /** output parameters */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    // water_block.addBodyStateForRecording<Vec3d>("KernelSummation");

    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** momentum equation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** mass equation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** pressure boundary condition. */
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);

    // Set up velocity condition

    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>
        viscous_acceleration(water_block_inner, water_block_contact);

    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact);

    /**
     * @brief Output.
     */
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.sph_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_centerline_velocity("Velocity", axial_velocity_observer_contact);
    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();

    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = params.end_time; /**< End time. */
    Real Output_Time = 0.25;         /**< Time stamps for output of body states. */
    Real dt = 0.0;                   /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    write_centerline_velocity.writeToFile(number_of_iterations);
    /**
     * @brief 	Main loop starts here.
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
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
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
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
        }
        TickCount t2 = TickCount::now();
        axial_velocity_observer_contact.updateConfiguration();
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
}

int main(int argc, char *argv[])
{
    const double Re = 500;
    FDA_nozzle_parameters params;
    params.number_of_particles = 5;
    std::stringstream wall_file;
    wall_file << "./input/FDA_nozzle_wall_N" << params.number_of_particles << ".stl";
    params.wall_file_path = wall_file.str();
    FDA_nozzle(argc, argv, params, Re);
    return 0;
}
