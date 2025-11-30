/**
 * @file 	Three ring impact.cpp
 * @brief 	This is the case file for the test of dynamic contacts between shell and solid.
 * @author  Weiyi Kong, Xiangyu Hu
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

namespace SPH
{
class Wall;
template <>
class ParticleGenerator<SurfaceParticles, Wall> : public ParticleGenerator<SurfaceParticles>
{
    const Vec3d center_;
    const Real length_;
    const Real width_;
    const Real dp_;

  public:
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, const Vec3d &center, Real length, Real width, Real dp)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          center_(center),
          length_(length),
          width_(width),
          dp_(dp) {};
    void prepareGeometricData() override
    {
        Real x = center_.x() - 0.5 * length_;
        while (x < center_.x() + 0.5 * length_)
        {
            Real z = -0.5 * width_;
            while (z < 0.5 * width_)
            {
                addPositionAndVolumetricMeasure(Vec3d(x, center_.y(), z), dp_ * dp_);
                addSurfaceProperties(Vec3d::UnitY(), dp_);
                z += dp_;
            }
            x += dp_;
        }
    }
};
} // namespace SPH

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

BoundingBoxd union_bounding_box(const BoundingBoxd &a, const BoundingBoxd &b)
{
    BoundingBoxd out = a;
    out.lower_[0] = std::min(a.lower_[0], b.lower_[0]);
    out.lower_[1] = std::min(a.lower_[1], b.lower_[1]);
    out.lower_[2] = std::min(a.lower_[2], b.lower_[2]);
    out.upper_[0] = std::max(a.upper_[0], b.upper_[0]);
    out.upper_[1] = std::max(a.upper_[1], b.upper_[1]);
    out.upper_[2] = std::max(a.upper_[2], b.upper_[2]);
    return out;
}

void block_sliding(
    int resolution_factor_cube,
    int resolution_factor_slope)
{
    // Global parameters
    constexpr Real threshold = 5e-2; // 5% error

    const Real end_time = 3.0;
    const Real scale = 1.0;

    // geometry
    const Real cube_length = 1.0 * scale;
    const Real slope_length = 25.0 * scale;
    const Real slope_width = 2.0 * scale;
    const Real slope_angle = 30.0 / 180.0 * Pi;
    const auto rotation = Rotation3d(slope_angle, Vec3d::UnitZ()).toRotationMatrix();
    auto rotation_inverse = rotation.transpose();

    // resolutions
    const Real resolution_cube = cube_length / (resolution_factor_cube * 5.0);
    const Real resolution_slope = cube_length / (resolution_factor_slope * 5.0);

    // Material properties
    const Real rho0_s = 2e3;
    const Real Youngs_modulus = 1e7;
    const Real poisson = 0.45;

    // Import meshes
    auto cube_translation = Vec3d(0.5 * cube_length + 2 * resolution_cube, 0.5 * cube_length, 0.0);
    auto mesh_cube = std::make_shared<TriangleMeshShapeBrick>(0.5 * cube_length * Vec3d::Ones(), 5, cube_translation, "Cube");
    auto slope_translation = Vec3d(0.5 * slope_length, -(0.65 * resolution_cube + 1.15 * resolution_slope), 0);
    auto mesh_slope = std::make_shared<TriangleMeshShapeBrick>(0.5 * Vec3d(slope_length, resolution_slope, slope_width), 5, slope_translation, "Slope");

    // Material models
    auto material_cube = makeShared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);

    // System bounding box
    BoundingBoxd bb_system = union_bounding_box(mesh_cube->getBounds(), mesh_slope->getBounds());

    // System
    SPHSystem system(bb_system, resolution_cube);

    // Create objects
    SolidBody cube_body(system, mesh_cube);
    cube_body.defineMaterial<SaintVenantKirchhoffSolid>(*material_cube.get());
    cube_body.generateParticles<BaseParticles, Lattice>();

    SolidBody slope_body(system, mesh_slope);
    slope_body.defineAdaptationRatios(1.15, resolution_cube / resolution_slope);
    slope_body.defineMaterial<SaintVenantKirchhoffSolid>(*material_cube.get());
    slope_body.generateParticles<SurfaceParticles, Wall>(slope_translation, slope_length, slope_width, resolution_slope);

    // Inner relation
    InnerRelation cube_inner(cube_body);

    // Methods
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(cube_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(cube_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(cube_inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        velocity_damping(0.5, cube_inner, "Velocity", get_physical_viscosity_general(rho0_s, Youngs_modulus, resolution_cube));
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(cube_body);

    // Contact relation
    SurfaceContactRelation contact_cube_to_slope(cube_body, {&slope_body}, {true});
    // Contact density
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor(contact_cube_to_slope);
    // Contact Force
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> contact_force(contact_cube_to_slope);

    // gravity
    const Real g = 9.81;
    Gravity gravity(rotation * (-g * Vec3d::UnitY()));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(cube_body, gravity);

    // analytical solution
    Real sin_theta = sin(slope_angle);
    Real cos_theta = cos(slope_angle);
    // analytical displacement under gravity
    auto get_analytical_displacement = [&](Real time)
    {
        Real a = 0.5 * g * sin_theta * time * time;
        Real u_x = a * cos_theta;
        Real u_y = a * sin_theta;
        return Vec3d(u_x, u_y, 0);
    };
    const Real max_error = threshold * get_analytical_displacement(end_time).norm();

    // Output
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.writeToFile(0);

    // Observer
    ObserverBody cube_observer(system, "CubeObserver");
    cube_observer.generateParticles<ObserverParticles>(StdVec<Vecd>{cube_translation});
    ContactRelation cube_observer_contact(cube_observer, {&cube_body});
    ObservedQuantityRecording<Vecd>
        write_cube_displacement("Position", cube_observer_contact);

    // initialize
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();
    constant_gravity.exec();

    // simulation
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 20.0;
    Real dt = 0.0;

    auto check_disp = [&]()
    {
        const Vec3d analytical_disp = get_analytical_displacement(physical_time);
        const Vec3d pos_observer = write_cube_displacement.getObservedQuantity()[0];
        const Vec3d rotated_disp = rotation_inverse * (pos_observer - cube_translation);
        for (int n = 0; n < 3; n++)
            ASSERT_NEAR(abs(rotated_disp[n]), abs(analytical_disp[n]), max_error);
    };

    TickCount t1 = TickCount::now();
    const Real dt_ref = computing_time_step_size.exec();
    auto run_simulation = [&]()
    {
        while (physical_time < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                contact_factor.exec();
                contact_force.exec();

                dt = computing_time_step_size.exec();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                stress_relaxation_first_half.exec(dt);
                velocity_damping.exec(dt);
                stress_relaxation_second_half.exec(dt);

                cube_body.updateCellLinkedList();
                contact_cube_to_slope.updateConfiguration();

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            write_cube_displacement.writeToFile(ite);
            check_disp();
            vtp_output.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        vtp_output.writeToFile(ite_output + 1);
        exit(0);
    }
}
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
TEST(sliding_3d, factor_2x_2x)
{
    block_sliding(2, 2);
}

int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}