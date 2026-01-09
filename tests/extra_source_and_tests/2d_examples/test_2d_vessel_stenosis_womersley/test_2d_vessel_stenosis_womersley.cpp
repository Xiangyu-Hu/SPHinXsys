/**
 * @file 	vessel_stenosis_womersley.cpp
 * @brief 	2D vessel stenosis with Womersley inflow example.
 * @details This is a test case for Womersley flow in a 2D vessel with stenosis
 * @author 	Minhui Zhou, Dong Wu
 */
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
// Select stenosis case (30% / 50% / 70%)
// By default we use the 30% stenosis CSVs.
// If you want 50% or 70% stenosis, change the filename suffix from "_0.3"
// to "_0.5" or "_0.7" accordingly (0.3=30%, 0.5=50%, 0.7=70%).
//----------------------------------------------------------------------
std::string womersley_velocity_profile_csv = "./input/womersley_velocity_profile_0.3.csv";
std::string outlet_pressure_csv = "./input/outlet_pressure_0.3.csv";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 0.0058;
Real DL1 = 4 * DH;
Real DL2 = 2 * DH;
Real DL3 = 16 * DH;
Real DL = DL1 + DL2 + DL3;
Real resolution_ref = DH / 12.0; 
Real BW = resolution_ref * 4.0; 
Real max_narrowing = 0.3; // chage to 0.5 or 0.7 for other stenosis cases
Real interpolationNum = 100;
BoundingBoxd system_domain_bounds(Vec2d(-DL1 - 0.5 * DL2 - BW, -0.5 * DH - BW), Vec2d(0.5 * DL2 + DL3, 0.5 * DH + BW));
//----------------------------------------------------------------------
//	Observation points
//----------------------------------------------------------------------
StdVec<Vecd> createAxialObservationPoints(Real full_length = DL, Vecd translation = Vecd(-DL1 - 0.5 * DL2, 0.0))
{
    StdVec<Vecd> observation_points;
    const int n_pts = 101;
    for (int i = 0; i < n_pts; ++i)
    {
        Real x = full_length * i / (n_pts - 1);
        Vecd point_coordinate(x, 0.0);
        observation_points.emplace_back(point_coordinate + translation);
    }
    return observation_points;
}
//----------------------------------------------------------------------
//	Material parameters of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1040.0;
Real mu_f = 0.004;
Real U_f = 4.44;
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Material parameters of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1060.0;
Real Youngs_modulus1 = 0.4e6;
Real poisson = 0.499;
//----------------------------------------------------------------------
//	Buffer parameters
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = Vec2d(-DL1 - 0.5 * DL2 + 2 * resolution_ref, 0);
Vec2d right_bidirectional_translation = Vec2d(0.5 * DL2 + DL3 - 2 * resolution_ref, 0);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

class OutletPressureCSV
{
  public:
    std::vector<Real> times;
    std::vector<Real> pressures;
    Real period{0.0};

    explicit OutletPressureCSV(const std::string &csv_file, Real period_override = -1.0)
    {
        std::ifstream file(csv_file);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, tok_t, tok_p;

        while (std::getline(file, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);

            if (!std::getline(ss, tok_t, ','))
                continue;
            if (!std::getline(ss, tok_p, ','))
                continue;
            try
            {
                Real t = static_cast<Real>(std::stod(tok_t));
                Real p = static_cast<Real>(std::stod(tok_p));
                times.push_back(t);
                pressures.push_back(p);
            }
            catch (...)
            {
                continue;
            }
        }

        if (!times.empty())
        {
            if (period_override > Real(0))
            {
                period = period_override; 
            }
            else
            {
                period = times.back() - times.front(); 
            }
        }

        std::cout << "[INFO] Loaded " << times.size()
                  << " pressure samples. period = " << period << " s" << std::endl;
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin())
            return 0;
        if (it == v.end())
            return v.size() - 1;
        size_t hi = static_cast<size_t>(std::distance(v.begin(), it));
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real t) const
    {
        if (times.empty())
            return Real(0);
        if (period > Real(0))
        {
            Real t0 = times.front();
            Real dt = t - t0;
            Real m = std::fmod(dt, period);
            if (m < Real(0))
                m += period;
            t = t0 + m;
        }

        size_t it = nearestIndex(times, t);
        return pressures[it];
    }
};

struct RightInflowPressure
{
    OutletPressureCSV pressure;

    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition,
                        const std::string &csv_path = outlet_pressure_csv,
                        Real T_override = 0.8)
        : pressure(csv_path, T_override) {}

    Real operator()(Real p, Real physical_time)
    {
        return 1000.0 * pressure.value(physical_time);
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
class WomersleyProfileCSV
{
  public:
    std::vector<Real> times;                   
    std::vector<Real> radii;                   
    std::vector<std::vector<Real>> velocities; 
    Real period{0.0};

    explicit WomersleyProfileCSV(const std::string &csv_file, Real period_override = -1.0)
    {
        std::ifstream file(csv_file);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, token;

        std::getline(file, line);
        {
            std::stringstream ss(line);
            std::getline(ss, token, ','); 
            while (std::getline(ss, token, ','))
                radii.push_back(static_cast<Real>(std::stod(token)));
        }

        while (std::getline(file, line))
        {
            if (line.empty())
                continue;
            std::stringstream ls(line);
            std::getline(ls, token, ','); // time
            if (token == "NaN" || token.empty())
                continue;
            times.push_back(static_cast<Real>(std::stod(token)));

            std::vector<Real> row;
            while (std::getline(ls, token, ','))
                row.push_back(static_cast<Real>(std::stod(token)));
            velocities.push_back(std::move(row));
        }

        if (period_override > 0.0)
        {
            period = period_override; 
        }
        else if (!times.empty())
        {
            period = times.back() - times.front();
        }

        std::cout << "[INFO] Loaded " << times.size()
                  << " time steps and " << radii.size()
                  << " radii. period = " << period << " s\n";
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin())
            return 0;
        if (it == v.end())
            return v.size() - 1;
        size_t hi = size_t(std::distance(v.begin(), it));
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real r, Real t) const
    {
        // clamp r
        if (r <= radii.front())
            r = radii.front();
        else if (r >= radii.back())
            r = radii.back();
        if (period > 0.0 && !times.empty())
        {
            Real t0 = times.front();
            Real dt = t - t0;
            Real m = std::fmod(dt, period);
            if (m < Real(0))
                m += period;
            t = t0 + m;
        }

        size_t it = nearestIndex(times, t);
        size_t ir = nearestIndex(radii, r);
        return velocities[it][ir];
    }
};

struct InflowVelocity
{
    WomersleyProfileCSV profile;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &,
                   const std::string &csv_path = womersley_velocity_profile_csv,
                   Real T_override = 0.8) 
        : profile(csv_path, T_override) {}

    Vecd operator()(Vecd &position, Vecd &, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        Real r = std::fabs(position[1]);
        target_velocity[0] = profile.value(r, current_time);
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	stenosis definition
//  a - maximaum narrowing, X0 - halflength of stenosis
//----------------------------------------------------------------------
Real outline(Real x_rel, Real a, Real X0)
{
    if (std::abs(x_rel) <= X0)
    {
        return 1.0 - 0.5 * a * (1.0 + std::cos(Pi * x_rel / X0));
    }
    else
    {
        return 1.0;
    }
}

std::vector<Vecd> createStenosisUpper(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N; 
    std::vector<Vecd> stenosis_upper;
    stenosis_upper.reserve(N + 1);

    for (int i = N; i >= 0; --i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = (D * 0.5) * y_rel;
        stenosis_upper.push_back(Vecd(x_rel, +y));
    }
    stenosis_upper.push_back(stenosis_upper.front());
    return stenosis_upper;
}

std::vector<Vecd> createStenosisLower(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N;
    std::vector<Vecd> stenosis_lower;
    stenosis_lower.reserve(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = -(D * 0.5) * y_rel;
        stenosis_lower.push_back(Vecd(x_rel, y));
    }
    stenosis_lower.push_back(stenosis_lower.front());
    return stenosis_lower;
}
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createBloodShape()
{
    // geometry
    std::vector<Vecd> blood_shape;
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    return blood_shape;
}
class Blood : public ComplexShape
{
  public:
    explicit Blood(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon blood_outer_boundary(createBloodShape());
        add<MultiPolygonShape>(blood_outer_boundary, "BloodOuterBoundary");

        MultiPolygon stenosisUpper(createStenosisUpper(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisUpper);
        MultiPolygon stenosisLower(createStenosisLower(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisLower);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    return outer_wall_shape;
}
std::vector<Vecd> createCompleteInnerWallShape(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    std::vector<Vecd> lumen;
    lumen.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * D));
    lumen.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * D));
    lumen.push_back(Vecd(-X0, 0.5 * D));
    Real dx = (2.0 * X0) / N;
    std::vector<Vecd> stenosis_lower1;
    stenosis_lower1.reserve(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = (D * 0.5) * y_rel;
        lumen.push_back(Vecd(x_rel, +y));
    }
    lumen.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * D));
    lumen.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * D));
    lumen.push_back(Vecd(X0, -0.5 * D));

    for (int i = N; i >= 0; --i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = -(D * 0.5) * y_rel;
        lumen.push_back(Vecd(x_rel, y));
    }

    lumen.push_back(lumen.front());
    return lumen;
}
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon boundary_outer_wall_shape(createOuterWallShape());
        add<MultiPolygonShape>(boundary_outer_wall_shape, "OuterWallBoundary");
        MultiPolygon boundary_inner_wall_shape(createCompleteInnerWallShape(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(boundary_inner_wall_shape, "OuterWallBoundary");
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // sph_system.setRunParticleRelaxation(true);
    // sph_system.setReloadParticles(false);

    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.setGenerateRegressionData(true);

    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody blood(sph_system, makeShared<Blood>("WaterBody"));
    blood.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet();
    blood.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(6.0);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? blood.generateParticlesWithReserve<BaseParticles, Reload>(in_outlet_particle_buffer, blood.getName())
        : blood.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet();
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody velocity_axial_observer(sph_system, "VelocityAxialObserver");
    velocity_axial_observer.defineAdaptationRatios(0.25, 1.0); 
    velocity_axial_observer.generateParticles<ObserverParticles>(createAxialObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation blood_inner(blood);
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation blood_contact(blood, {&wall_boundary});
    ContactRelation wall_boundary_contact(wall_boundary, {&blood});
    ContactRelation velocity_observer_contact_axial(velocity_axial_observer, {&blood});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation blood_complex(blood_inner, blood_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_blood_particles(blood);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_blood(blood_inner);
        BodyStatesRecordingToVtp write_wall_boundary_and_blood(sph_system);
        ReloadParticleIO write_particle_reload_files({&wall_boundary, &blood});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_wall_boundary_particles.exec(0.25);
        random_blood_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_inner_blood.SurfaceBounding().exec();
        write_wall_boundary_and_blood.writeToFile();

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            relaxation_step_inner_blood.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_wall_boundary_and_blood.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        write_particle_reload_files.writeToFile();
        return 0;
    }

    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NablaWVComplex> kernel_summation(blood_inner, blood_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_boundary_corrected_configuration(wall_boundary_inner);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(blood_inner, blood_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(blood_inner, blood_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(blood_inner, blood_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(blood_inner, blood_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(blood_inner, blood_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(blood, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(blood);
    //----------------------------------------------------------------------
    //  Buffer
    //----------------------------------------------------------------------
    AlignedBox left_emitter_shape(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell left_emitter(blood, left_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);

    AlignedBox right_emitter_shape(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell right_emitter(blood, right_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(blood_inner, blood_contact);

    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(blood);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(blood, "Pressure");
    body_states_recording.addToWrite<int>(blood, "Indicator");
    body_states_recording.addToWrite<Real>(blood, "Density");
    body_states_recording.addToWrite<int>(blood, "BufferIndicator");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_centerline_velocity_axial("Velocity", velocity_observer_contact_axial);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    wall_boundary_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    if (sph_system.RestartStep() != 0)
    {
        physical_time = restart_io.readRestartFiles(sph_system.RestartStep());
        blood.updateCellLinkedList();
        blood_complex.updateConfiguration();
        velocity_observer_contact_axial.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control   
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 1.61;  
    Real Output_Time = 0.01; 
    Real dt = 0.0;       
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_centerline_velocity_axial.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();

            viscous_acceleration.exec();
            transport_velocity_correction.exec();

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
                physical_time += dt;
            }
            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            blood.updateCellLinkedList();
            blood_complex.updateConfiguration();
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        velocity_observer_contact_axial.updateConfiguration();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;

        /** Update observer and write output of observer. */
        write_centerline_velocity_axial.writeToFile(number_of_iterations);
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_centerline_velocity_axial.generateDataBase(0.05); 
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_centerline_velocity_axial.testResult();
    }

    return 0;
}
