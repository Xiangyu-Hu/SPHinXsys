/**
 * @file 	2d_gradient.cpp
 * @brief 	Accuracy test of the linearly corrected gradient and the second-order
 * 			(Hessian-based) gradient on disordered particles.
 * @details	This test is derived from unit_test_gradient_ck/2d_gradient.cpp.
 * 			The fluid particles are shifted by independent, seeded random displacements.
 * 			The gradient and the Hessian of a quadratic and of a sinusoidal field are compared
 * 			with the analytical values at all fluid particles. The errors are evaluated in the
 * 			interior, which is at least the kernel support radius away from all boundaries,
 * 			and in a fixed window, which is the same for all resolutions.
 * 			The following optional environment variables allow convergence studies:
 * 			GRADIENT_DP: particle spacing, default 0.01, at most 0.04,
 * 			GRADIENT_JITTER: amplitude of the random shift relative to the spacing, default 0.25.
 * 			The lines starting with "[summary]" report the RMS and the maximum errors.
 * @author 	Xiangyu Hu (original test)
 */
#include "sphinxsys.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <type_traits>
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
// Amplitude of the random particle shift relative to the spacing, and the seed of the generator.
Real jitter_ratio = 0.25;
constexpr unsigned jitter_seed = 1234;
// Margin of the fixed evaluation window. It equals the kernel support radius (2.6 spacings)
// of the coarsest resolution, so the window lies inside the interior for all resolutions.
constexpr Real window_margin = 0.104;
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
//	Test fields and their analytical gradients and Hessians.
//	The Hessian is stored as [phi_xx, phi_yy, 2 phi_xy], the same layout as in the operators.
//----------------------------------------------------------------------
// Quadratic field with fixed and distinct coefficients, so that a wrong order of the
// Hessian components is detected. The reference Hessian is [1.0, -0.8, 1.6].
Vec2d first_coefficient(0.3, -0.7);
VecMat2d second_coefficient{0.5, -0.4, 0.8};
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
Vec2d quadraticGradient(const Vec2d &position)
{
    return first_coefficient +
           Vec2d(2.0 * second_coefficient[0] * position[0] + second_coefficient[2] * position[1],
                 2.0 * second_coefficient[1] * position[1] + second_coefficient[2] * position[0]);
}
VecMat2d quadraticHessian()
{
    return 2.0 * second_coefficient;
}

// Smooth, non-polynomial field: phi = sin(2 pi x) cos(pi y).
class SinusoidalProfile : public ReturnFunction<Real>
{
  public:
    Real operator()(const Vec2d &position)
    {
        return std::sin(2.0 * Pi * position[0]) * std::cos(Pi * position[1]);
    }
};
Vec2d sinusoidalGradient(const Vec2d &position)
{
    return Vec2d(2.0 * Pi * std::cos(2.0 * Pi * position[0]) * std::cos(Pi * position[1]),
                 -Pi * std::sin(2.0 * Pi * position[0]) * std::sin(Pi * position[1]));
}
VecMat2d sinusoidalHessian(const Vec2d &position)
{
    return VecMat2d(-4.0 * Pi * Pi * std::sin(2.0 * Pi * position[0]) * std::cos(Pi * position[1]),
                    -Pi * Pi * std::sin(2.0 * Pi * position[0]) * std::cos(Pi * position[1]),
                    -4.0 * Pi * Pi * std::cos(2.0 * Pi * position[0]) * std::sin(Pi * position[1]));
}
//----------------------------------------------------------------------
//	Google test items: point-wise checks at one observation point, as in the original test.
//----------------------------------------------------------------------
const Vec2d observation_point(0.37, 0.21);

Matd approximated_gradient(Matd::Zero());
Matd reference_gradient(Matd::Identity());
TEST(LinearGradient, Error)
{
    EXPECT_LT((reference_gradient - approximated_gradient).norm(), 1.0e-6);
    std::cout << "Reference Gradient: " << reference_gradient << " and "
              << "Predicted Gradient: " << approximated_gradient << std::endl;
};

VecMat2d approximated_hessian = VecMat2d::Zero();
VecMat2d reference_hessian = quadraticHessian();
TEST(Hessian, Error)
{
    EXPECT_LT((reference_hessian - approximated_hessian).norm(), 1.0e-6);
    std::cout << "Reference Hessian: " << reference_hessian << " and "
              << "Predicted Hessian: " << approximated_hessian << std::endl;
};

Vec2d approximated_2nd_order_gradient = Vec2d::Zero();
Vec2d reference_2nd_order_gradient = quadraticGradient(observation_point);
TEST(SecondOrderGradient, Error)
{
    EXPECT_LT((reference_2nd_order_gradient - approximated_2nd_order_gradient).norm(), 1.0e-6);
    std::cout << "Reference Second Order Gradient: " << reference_2nd_order_gradient << " and "
              << "Predicted Second Order Gradient: " << approximated_2nd_order_gradient << std::endl;
};

// The analytical references are checked against central finite differences, independent of SPH.
// The tolerance is dominated by the truncation error of the finite differences, about 3e-4.
template <class Field>
void expectAnalyticalReference(Field field, const Vec2d &gradient, const VecMat2d &hessian, const Vec2d &position)
{
    const Real step = 1.0e-3;
    const Vec2d ex(step, 0.0), ey(0.0, step);
    const Real f0 = field(position);
    Vec2d fd_gradient((field(position + ex) - field(position - ex)) / (2.0 * step),
                      (field(position + ey) - field(position - ey)) / (2.0 * step));
    Real fxx = (field(position + ex) - 2.0 * f0 + field(position - ex)) / (step * step);
    Real fyy = (field(position + ey) - 2.0 * f0 + field(position - ey)) / (step * step);
    Real fxy = (field(position + ex + ey) - field(position + ex - ey) -
                field(position - ex + ey) + field(position - ex - ey)) /
               (4.0 * step * step);
    VecMat2d fd_hessian(fxx, fyy, 2.0 * fxy);
    EXPECT_LT((gradient - fd_gradient).norm(), 2.0e-3);
    EXPECT_LT((hessian - fd_hessian).norm(), 2.0e-3);
}
TEST(AnalyticalReference, FiniteDifference)
{
    // the finite difference tolerance assumes double precision
    if (!std::is_same<Real, double>::value)
        return;
    for (const Vec2d &position : {Vec2d(0.37, 0.21), Vec2d(0.81, 0.07), Vec2d(0.05, 0.48)})
    {
        expectAnalyticalReference(ParabolicProfile(), quadraticGradient(position),
                                  quadraticHessian(), position);
        expectAnalyticalReference(SinusoidalProfile(), sinusoidalGradient(position),
                                  sinusoidalHessian(position), position);
    }
}
//----------------------------------------------------------------------
//	Particle disorder: independent displacements from a seeded generator.
//	RandomizeParticlePosition is not used since it draws a clock-seeded generator for
//	every component, which can give strongly correlated displacements.
//----------------------------------------------------------------------
void applySeededJitter(Vecd *position, size_t total_particles, Real amplitude, unsigned seed)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<Real> distribution(-1.0, 1.0);
    for (size_t i = 0; i != total_particles; ++i)
        for (int k = 0; k != position[i].size(); ++k)
            position[i][k] += amplitude * distribution(generator);
}
// Pearson correlation between a[i] and b[i + lag].
Real correlation(const StdVec<Real> &a, const StdVec<Real> &b, size_t lag)
{
    const size_t count = a.size() - lag;
    Real mean_a = 0.0, mean_b = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        mean_a += a[i];
        mean_b += b[i + lag];
    }
    mean_a /= Real(count);
    mean_b /= Real(count);
    Real covariance = 0.0, variance_a = 0.0, variance_b = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        covariance += (a[i] - mean_a) * (b[i + lag] - mean_b);
        variance_a += (a[i] - mean_a) * (a[i] - mean_a);
        variance_b += (b[i + lag] - mean_b) * (b[i + lag] - mean_b);
    }
    return covariance / std::sqrt(variance_a * variance_b);
}
//----------------------------------------------------------------------
//	Operator results and error statistics.
//----------------------------------------------------------------------
// Results of the three operators, one entry for every fluid particle. SecondOrderGradient
// overwrites the variable written by LinearGradient, so each result is copied right after
// its operator has been executed.
struct GradientResults
{
    StdVec<Vecd> linear_gradient;
    StdVec<VecMat2d> hessian;
    StdVec<Vecd> second_order_gradient;
};
template <typename DataType>
StdVec<DataType> copyParticleVariable(BaseParticles &particles, const std::string &name, size_t total_particles)
{
    DataType *data = particles.getVariableDataByName<DataType>(name);
    return StdVec<DataType>(data, data + total_particles);
}

struct ErrorStatistics
{
    size_t count = 0;
    Real sum_square = 0.0;
    Real max_error = 0.0;
    void add(Real error)
    {
        ++count;
        sum_square += error * error;
        max_error = SMAX(max_error, error);
    }
    Real rms() const { return count == 0 ? Real(0.0) : Real(std::sqrt(sum_square / Real(count))); }
};
struct OperatorErrors
{
    ErrorStatistics linear_gradient, second_order_gradient, hessian;
};
// The interior is at least the kernel support radius away from all four boundaries.
// The window is a fixed region, the same for all resolutions.
struct FieldErrors
{
    OperatorErrors interior, window;
};
bool isInInterior(const Vecd &position, Real support_radius)
{
    return position[0] >= support_radius && position[0] <= width - support_radius &&
           position[1] >= support_radius && position[1] <= height - support_radius;
}
bool isInWindow(const Vecd &position)
{
    return position[0] >= window_margin && position[0] <= width - window_margin &&
           position[1] >= window_margin && position[1] <= height - window_margin;
}
template <class GradientReference, class HessianReference>
FieldErrors computeFieldErrors(const GradientResults &results, const StdVec<Vecd> &position, Real support_radius,
                               GradientReference reference_gradient, HessianReference reference_hessian)
{
    FieldErrors field_errors;
    for (size_t i = 0; i != position.size(); ++i)
    {
        const Vecd gradient_reference = reference_gradient(position[i]);
        const VecMat2d hessian_reference = reference_hessian(position[i]);
        const Real linear_error = (results.linear_gradient[i] - gradient_reference).norm();
        const Real second_order_error = (results.second_order_gradient[i] - gradient_reference).norm();
        const Real hessian_error = (results.hessian[i] - hessian_reference).norm();
        if (isInInterior(position[i], support_radius))
        {
            field_errors.interior.linear_gradient.add(linear_error);
            field_errors.interior.second_order_gradient.add(second_order_error);
            field_errors.interior.hessian.add(hessian_error);
        }
        if (isInWindow(position[i]))
        {
            field_errors.window.linear_gradient.add(linear_error);
            field_errors.window.second_order_gradient.add(second_order_error);
            field_errors.window.hessian.add(hessian_error);
        }
    }
    return field_errors;
}
void printSummary(const std::string &field_name, const FieldErrors &field_errors, Real spacing, Real jitter)
{
    const OperatorErrors *regions[2] = {&field_errors.window, &field_errors.interior};
    const char *region_names[2] = {"window", "interior"};
    char line[512];
    for (int k = 0; k != 2; ++k)
    {
        const OperatorErrors &e = *regions[k];
        std::snprintf(line, sizeof(line),
                      "[summary] dp=%g jitter=%g field=%s region=%s n=%zu "
                      "linear_rms=%.6e second_rms=%.6e hessian_rms=%.6e "
                      "linear_linf=%.6e second_linf=%.6e hessian_linf=%.6e",
                      double(spacing), double(jitter), field_name.c_str(), region_names[k],
                      e.linear_gradient.count, double(e.linear_gradient.rms()),
                      double(e.second_order_gradient.rms()), double(e.hessian.rms()),
                      double(e.linear_gradient.max_error), double(e.second_order_gradient.max_error),
                      double(e.hessian.max_error));
        std::cout << line << std::endl;
    }
}
//----------------------------------------------------------------------
//	Google test items: accuracy on all fluid particles. The results are computed in main.
//	The tolerances assume double precision, since the correction matrices are regularized
//	with a constant of about 1e-8. The assertions on the benefit of the second-order gradient
//	apply to a substantial particle disorder (jitter of at least 0.2 spacings).
//----------------------------------------------------------------------
struct AccuracyResults
{
    Real jitter_correlation_xy = 0.0;       // correlation of the x and y displacements
    Real jitter_correlation_neighbor = 0.0; // correlation of the x displacements of consecutive particles
    size_t jitter_sample_size = 0;
    FieldErrors phi; // quadratic field
    FieldErrors psi; // sinusoidal field
} accuracy;

TEST(ParticleDisorder, IndependentDisplacements)
{
    if (jitter_ratio <= 0.0)
        return;
    // five standard deviations of the sampling error of a correlation coefficient
    const Real tolerance = 5.0 / std::sqrt(Real(accuracy.jitter_sample_size));
    EXPECT_LT(std::abs(accuracy.jitter_correlation_xy), tolerance);
    EXPECT_LT(std::abs(accuracy.jitter_correlation_neighbor), tolerance);
}
TEST(QuadraticField, SecondOrderGradientAndHessianAreExact)
{
    if (!std::is_same<Real, double>::value)
        return;
    const OperatorErrors &interior = accuracy.phi.interior;
    EXPECT_GT(interior.linear_gradient.count, size_t(0));
    EXPECT_LT(interior.second_order_gradient.max_error, 1.0e-6);
    EXPECT_LT(interior.hessian.max_error, 1.0e-4);
    // The linearly corrected gradient is only first-order consistent, so it is not exact on
    // disordered particles. This also shows that the particle disorder is effective.
    if (jitter_ratio >= 0.2)
    {
        EXPECT_GT(interior.linear_gradient.rms(), 100.0 * interior.second_order_gradient.rms());
    }
}
TEST(SinusoidalField, SecondOrderGradientIsMoreAccurateOnDisorderedParticles)
{
    if (jitter_ratio < 0.2)
        return;
    const OperatorErrors &interior = accuracy.psi.interior;
    EXPECT_LT(interior.second_order_gradient.rms(), interior.linear_gradient.rms());
    // at the default and finer resolutions the error is at least halved
    if (particle_spacing <= 0.0105)
    {
        EXPECT_LT(interior.second_order_gradient.rms(), 0.5 * interior.linear_gradient.rms());
    }
}
//----------------------------------------------------------------------
//	Reads a real number from an environment variable. An invalid or out-of-range value stops the test.
//----------------------------------------------------------------------
Real readEnvironmentReal(const char *name, Real default_value, Real lower, Real upper)
{
    const char *text = std::getenv(name);
    if (text == nullptr)
        return default_value;
    char *end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || !(value >= lower && value <= upper))
    {
        std::cerr << "Invalid value \"" << text << "\" of " << name
                  << ": expected a number in [" << lower << ", " << upper << "]." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return Real(value);
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    particle_spacing = readEnvironmentReal("GRADIENT_DP", particle_spacing, 0.001, 0.04);
    boundary_width = particle_spacing * 4;
    jitter_ratio = readEnvironmentReal("GRADIENT_JITTER", jitter_ratio, 0.0, 0.4);
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
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    water_block.getBaseParticles().registerStateVariable<Real>("Phi");
    water_block.getBaseParticles().registerStateVariable<Real>("Psi");
    // shift the fluid particles randomly to avoid the symmetry of the lattice
    BaseParticles &water_particles = water_block.getBaseParticles();
    const size_t total_water_particles = water_particles.TotalRealParticles();
    Vecd *water_position = water_particles.getVariableDataByName<Vecd>("Position");
    const Real jitter_amplitude = jitter_ratio * particle_spacing;
    if (jitter_amplitude > 0.0)
    {
        const StdVec<Vecd> lattice_position(water_position, water_position + total_water_particles);
        applySeededJitter(water_position, total_water_particles, jitter_amplitude, jitter_seed);
        StdVec<Real> dx(total_water_particles), dy(total_water_particles);
        for (size_t i = 0; i != total_water_particles; ++i)
        {
            dx[i] = (water_position[i][0] - lattice_position[i][0]) / jitter_amplitude;
            dy[i] = (water_position[i][1] - lattice_position[i][1]) / jitter_amplitude;
        }
        accuracy.jitter_sample_size = total_water_particles;
        accuracy.jitter_correlation_xy = correlation(dx, dy, 0);
        accuracy.jitter_correlation_neighbor = correlation(dx, dx, 1);
        std::cout << "[disorder] shift amplitude " << jitter_ratio << " spacings, correlation of dx and dy: "
                  << accuracy.jitter_correlation_xy << ", of dx of consecutive particles: "
                  << accuracy.jitter_correlation_neighbor << std::endl;
    }
    else
        std::cout << "[disorder] disabled, the particles stay on the lattice" << std::endl;

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();
    wall.getBaseParticles().registerStateVariable<Real>("Phi");
    wall.getBaseParticles().registerStateVariable<Real>("Psi");

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vec2d> observation_location = {observation_point};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, wall);
    Contact<> fluid_observer_contact(fluid_observer, water_block);
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody>
        water_cell_linked_list(water_block);
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
        observed_position_gradient(fluid_observer_contact, "PositionGradient");

    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>, Contact<>>>
        displacement_matrix_gradient(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        hessian_correction_matrix(DynamicsArgs(water_block_inner, 0.0), water_wall_contact);

    // the quadratic field "Phi" and the sinusoidal field "Psi", assigned to the fluid and the wall particles
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ParabolicProfile>>>
        water_block_initial_condition(water_block, "Phi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ParabolicProfile>>>
        wall_initial_condition(wall, "Phi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<SinusoidalProfile>>>
        water_block_sinusoidal_condition(water_block, "Psi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<SinusoidalProfile>>>
        wall_sinusoidal_condition(wall, "Psi");

    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        variable_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>, Contact<Real>>>
        variable_hessian(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian(fluid_observer_contact, "PhiHessian");

    InteractionDynamicsCK<MainExecutionPolicy, SecondOrderGradient<Inner<Real>, Contact<Real>>>
        variable_2nd_order_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    ObservedQuantityRecording<MainExecutionPolicy, Vec2d, RestoringCorrection>
        observed_2nd_order_gradient(fluid_observer_contact, "PhiGradient");

    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        psi_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Psi")),
            DynamicsArgs(water_wall_contact, std::string("Psi")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>, Contact<Real>>>
        psi_hessian(
            DynamicsArgs(water_block_inner, std::string("Psi")),
            DynamicsArgs(water_wall_contact, std::string("Psi")));
    InteractionDynamicsCK<MainExecutionPolicy, SecondOrderGradient<Inner<Real>, Contact<Real>>>
        psi_2nd_order_gradient(
            DynamicsArgs(water_block_inner, std::string("Psi")),
            DynamicsArgs(water_wall_contact, std::string("Psi")));
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
    water_block_sinusoidal_condition.exec();
    wall_sinusoidal_condition.exec();

    linear_correction_matrix.exec();
    position_linear_gradient.exec();
    observed_position_gradient.writeToFile(0);
    approximated_gradient = *observed_position_gradient.getObservedQuantity();

    displacement_matrix_gradient.exec();
    hessian_correction_matrix.exec();

    // The results of every operator are copied right after its execution.
    GradientResults phi_results, psi_results;
    variable_linear_gradient.exec();
    phi_results.linear_gradient = copyParticleVariable<Vecd>(water_particles, "PhiGradient", total_water_particles);
    variable_hessian.exec();
    phi_results.hessian = copyParticleVariable<VecMat2d>(water_particles, "PhiHessian", total_water_particles);
    observed_hessian.writeToFile(0);
    approximated_hessian = *observed_hessian.getObservedQuantity();

    variable_2nd_order_gradient.exec();
    phi_results.second_order_gradient = copyParticleVariable<Vecd>(water_particles, "PhiGradient", total_water_particles);
    observed_2nd_order_gradient.writeToFile(0);
    approximated_2nd_order_gradient = *observed_2nd_order_gradient.getObservedQuantity();

    psi_linear_gradient.exec();
    psi_results.linear_gradient = copyParticleVariable<Vecd>(water_particles, "PsiGradient", total_water_particles);
    psi_hessian.exec();
    psi_results.hessian = copyParticleVariable<VecMat2d>(water_particles, "PsiHessian", total_water_particles);
    psi_2nd_order_gradient.exec();
    psi_results.second_order_gradient = copyParticleVariable<Vecd>(water_particles, "PsiGradient", total_water_particles);

    // errors against the analytical references at all fluid particles
    const Real support_radius = water_block.getSPHAdaptation().getKernel()->CutOffRadius();
    const StdVec<Vecd> fluid_position(water_position, water_position + total_water_particles);
    accuracy.phi = computeFieldErrors(
        phi_results, fluid_position, support_radius,
        [](const Vecd &position) { return quadraticGradient(position); },
        [](const Vecd &) { return quadraticHessian(); });
    accuracy.psi = computeFieldErrors(
        psi_results, fluid_position, support_radius,
        [](const Vecd &position) { return sinusoidalGradient(position); },
        [](const Vecd &position) { return sinusoidalHessian(position); });
    printSummary("Phi", accuracy.phi, particle_spacing, jitter_ratio);
    printSummary("Psi", accuracy.psi, particle_spacing, jitter_ratio);

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
