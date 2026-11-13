/**
 * @file 	relaxation_consistency_circle.cpp
 * @brief 	Standard particle relaxation in a circle, its zero-order consistency, and the accuracy of the
 * 			LinearGradient and SecondOrderGradient operators on the relaxed particles.
 * @details	A circle of radius 1 is filled with lattice particles. They are shifted by independent, seeded
 * 			random displacements and relaxed with the standard relaxation of SPHinXsys
 * 			(RelaxationStepLevelSetCorrectionInner), until the RMS of the dimensionless zero-order residual
 * 			h |R_P| in the interior falls below a target, or a step limit is reached. The interior is at least
 * 			the kernel support radius away from the boundary. The residuals are
 * 			R_P = sum_j V_j grad W_ij               of the standard kernel gradient,
 * 			R_B = sum_j V_j (B_i + B_j) grad W_ij   of the corrected form, with the matrix B of
 * 			                                        LinearCorrectionMatrix.
 * 			Both are evaluated at the same particle positions. The standard relaxation reduces R_P. R_B is
 * 			only reported, it is not tested.
 * 			On the relaxed particles, the gradient of a constant, a linear and a quadratic field is computed
 * 			with LinearGradient and SecondOrderGradient and compared with the analytical gradient in a fixed
 * 			disk |x| <= 0.45 and in the interior. Expected accuracy: both operators are exact for constant and
 * 			linear fields. For the quadratic field LinearGradient is only first-order consistent, while
 * 			SecondOrderGradient is exact.
 * 			Optional environment variables: RELAX_DP (particle spacing, default 0.05),
 * 			RELAX_TARGET (target residual, default 1e-3) and RELAX_MAX_STEPS (step limit, default 3000).
 * 			The lines starting with "[setup]", "[relax]", "[result]" and "[gradient]" are also written to
 * 			relaxation_consistency_dp<spacing>.txt in the working directory. For plotting, the same results
 * 			are written to relaxation_consistency_dp<spacing>_history.csv (residuals at every check),
 * 			_gradient.csv (gradient errors) and _particles.csv (initial and final particle positions).
 */
#include "sphinxsys.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <type_traits>
using namespace SPH;
//----------------------------------------------------------------------
//	Basic parameters.
//----------------------------------------------------------------------
const Real circle_radius = 1.0;
Real particle_spacing = 0.05;
Real target_residual = 1.0e-3;     // target of the RMS of h |R_P| in the interior
int max_relaxation_steps = 3000;   // step limit of the relaxation
constexpr int check_interval = 50; // the residual is evaluated every check_interval steps
constexpr Real jitter_ratio = 0.25;
constexpr unsigned jitter_seed = 1234;
// The errors of the gradients are also evaluated in a fixed disk, the same for all resolutions.
// All neighbors of its particles have their full kernel support inside the circle.
constexpr Real core_radius = 0.45;
using MainExecutionPolicy = execution::ParallelPolicy;
//----------------------------------------------------------------------
//	Shape of the body.
//----------------------------------------------------------------------
class CircleBody : public ComplexShape
{
  public:
    explicit CircleBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(Vecd::Zero(), circle_radius);
    }
};
//----------------------------------------------------------------------
//	Test fields and their analytical gradients. The coefficients are those of the
//	quadratic field in unit_test_gradient_accuracy_ck.
//----------------------------------------------------------------------
Vec2d first_coefficient(0.3, -0.7);
VecMat2d second_coefficient{0.5, -0.4, 0.8};
class ConstantProfile : public ReturnFunction<Real>
{
  public:
    Real operator()(const Vec2d &position) { return 1.0; }
};
class LinearProfile : public ReturnFunction<Real>
{
  public:
    Real operator()(const Vec2d &position) { return first_coefficient.dot(position); }
};
class ParabolicProfile : public ReturnFunction<Real>
{
  public:
    Real operator()(const Vec2d &position)
    {
        return first_coefficient.dot(position) + second_coefficient.dot(vectorizeTensorSquare(position));
    }
};
Vec2d constantGradient(const Vec2d &position) { return Vec2d::Zero(); }
Vec2d linearGradient(const Vec2d &position) { return first_coefficient; }
Vec2d quadraticGradient(const Vec2d &position)
{
    return first_coefficient +
           Vec2d(2.0 * second_coefficient[0] * position[0] + second_coefficient[2] * position[1],
                 2.0 * second_coefficient[1] * position[1] + second_coefficient[2] * position[0]);
}
//----------------------------------------------------------------------
//	Seeded random displacements. RandomizeParticlePosition is not used since it is clock-seeded,
//	so that the particle configuration, and with it the results, would differ from run to run.
//----------------------------------------------------------------------
void applySeededJitter(Vecd *position, size_t total_particles, Real amplitude, unsigned seed)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<Real> distribution(-1.0, 1.0);
    for (size_t i = 0; i != total_particles; ++i)
        for (int k = 0; k != position[i].size(); ++k)
            position[i][k] += amplitude * distribution(generator);
}
//----------------------------------------------------------------------
//	Diagnostics of the particle configuration at the current positions.
//----------------------------------------------------------------------
// The zero-order residual of the kernel gradient, the kernel sum, the nearest neighbor distance
// and the number of neighbors.
class ZeroOrderConsistency : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ZeroOrderConsistency(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          kernel_gradient_sum_(particles_->registerStateVariableData<Vecd>("KernelGradientSum")),
          kernel_sum_(particles_->registerStateVariableData<Real>("KernelSum")),
          neighbor_number_(particles_->registerStateVariableData<Real>("NeighborNumber")),
          nearest_distance_(particles_->registerStateVariableData<Real>("NearestNeighborDistance")),
          W0_(getSPHAdaptation().getKernel()->W0(ZeroVecd)) {};
    virtual ~ZeroOrderConsistency() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd gradient_sum = Vecd::Zero();
        Real kernel_sum = Vol_[index_i] * W0_; // the particle itself
        Real nearest_distance = MaxReal;
        const Neighborhood &neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            const size_t index_j = neighborhood.j_[n];
            gradient_sum += neighborhood.dW_ij_[n] * Vol_[index_j] * neighborhood.e_ij_[n];
            kernel_sum += neighborhood.W_ij_[n] * Vol_[index_j];
            nearest_distance = SMIN(nearest_distance, neighborhood.r_ij_[n]);
        }
        kernel_gradient_sum_[index_i] = gradient_sum;
        kernel_sum_[index_i] = kernel_sum;
        neighbor_number_[index_i] = Real(neighborhood.current_size_);
        nearest_distance_[index_i] = nearest_distance;
    };

  protected:
    Real *Vol_;
    Vecd *kernel_gradient_sum_;
    Real *kernel_sum_;
    Real *neighbor_number_;
    Real *nearest_distance_;
    Real W0_;
};
// The zero-order residual of the corrected form, R_B = sum_j V_j (B_i + B_j) grad W_ij.
// The matrix B is the one of LinearCorrectionMatrix, which must be executed before.
class CorrectionMatrixResidual : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit CorrectionMatrixResidual(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          B_(particles_->getVariableDataByName<Matd>("LinearCorrectionMatrix")),
          residual_(particles_->registerStateVariableData<Vecd>("CorrectionMatrixResidual")) {};
    virtual ~CorrectionMatrixResidual() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd residual = Vecd::Zero();
        const Neighborhood &neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            const size_t index_j = neighborhood.j_[n];
            residual += (B_[index_i] + B_[index_j]) * neighborhood.e_ij_[n] *
                        (neighborhood.dW_ij_[n] * Vol_[index_j]);
        }
        residual_[index_i] = residual;
    };

  protected:
    Real *Vol_;
    Matd *B_;
    Vecd *residual_;
};
//----------------------------------------------------------------------
//	Statistics and results.
//----------------------------------------------------------------------
struct Statistics
{
    size_t count = 0;
    Real sum_square = 0.0;
    Real max_value = 0.0;
    void add(Real value)
    {
        ++count;
        sum_square += value * value;
        max_value = SMAX(max_value, value);
    }
    Real rms() const { return count == 0 ? Real(0.0) : Real(std::sqrt(sum_square / Real(count))); }
};
struct ConsistencyMetrics
{
    Statistics residual_p_interior, residual_b_interior; // h |R_P| and h |R_B| in the interior
    Statistics residual_p_band, residual_b_band;         // in the boundary band
    Statistics kernel_sum_interior;                      // |kernel sum - 1| in the interior
    Real nearest_min = 0.0;                              // minimum nearest neighbor distance in spacings
    Real neighbors_min = 0.0;                            // minimum number of neighbors
};
struct OperatorErrors
{
    Statistics linear_gradient, second_order_gradient;
};
struct FieldErrors
{
    OperatorErrors core, interior;
};
struct RunResults
{
    int steps = 0;
    bool converged = false;
    ConsistencyMetrics initial, final;
    FieldErrors constant, linear, quadratic;
} results;
struct HistoryRow
{
    int step;
    ConsistencyMetrics metrics;
};
StdVec<HistoryRow> history; // the diagnostics at step 0 and at every check
//----------------------------------------------------------------------
//	Output: every line is printed and also collected for the result file.
//----------------------------------------------------------------------
std::string output_text;
void emit(const char *line)
{
    std::cout << line << std::endl;
    output_text += line;
    output_text += "\n";
}
void printMetrics(const char *tag, int step, const ConsistencyMetrics &m)
{
    char line[640];
    std::snprintf(line, sizeof(line),
                  "[%s] step=%d interior_n=%zu hR_P_rms=%.3e hR_P_max=%.3e hR_B_rms=%.3e hR_B_max=%.3e "
                  "kernel_sum_dev_rms=%.3e kernel_sum_dev_max=%.3e band_hR_P_rms=%.3e band_hR_B_rms=%.3e "
                  "nearest_min=%.3f neighbors_min=%.0f",
                  tag, step, m.residual_p_interior.count, double(m.residual_p_interior.rms()),
                  double(m.residual_p_interior.max_value), double(m.residual_b_interior.rms()),
                  double(m.residual_b_interior.max_value), double(m.kernel_sum_interior.rms()),
                  double(m.kernel_sum_interior.max_value), double(m.residual_p_band.rms()),
                  double(m.residual_b_band.rms()), double(m.nearest_min), double(m.neighbors_min));
    emit(line);
}
void printFieldErrors(const char *field, const FieldErrors &e)
{
    const OperatorErrors *regions[2] = {&e.core, &e.interior};
    const char *region_names[2] = {"core", "interior"};
    char line[512];
    for (int k = 0; k != 2; ++k)
    {
        const OperatorErrors &r = *regions[k];
        std::snprintf(line, sizeof(line),
                      "[gradient] field=%s region=%s n=%zu linear_rms=%.6e linear_max=%.6e "
                      "second_order_rms=%.6e second_order_max=%.6e",
                      field, region_names[k], r.linear_gradient.count, double(r.linear_gradient.rms()),
                      double(r.linear_gradient.max_value), double(r.second_order_gradient.rms()),
                      double(r.second_order_gradient.max_value));
        emit(line);
    }
}
//----------------------------------------------------------------------
//	Google test items. The results are computed in main.
//	The tolerances assume double precision, since the correction matrices are regularized with a
//	constant of about 1e-8.
//----------------------------------------------------------------------
TEST(Relaxation, ReachesTheTargetResidual)
{
    EXPECT_TRUE(results.converged);
    EXPECT_GT(results.final.residual_p_interior.count, size_t(0));
    EXPECT_LT(results.final.residual_p_interior.rms(), target_residual);
    // the relaxation removes most of the residual of the randomly shifted particles
    EXPECT_LT(results.final.residual_p_interior.rms(), 0.05 * results.initial.residual_p_interior.rms());
}
TEST(Relaxation, ParticleDistributionIsRegular)
{
    // the kernel sum is close to one, and no particles come close to each other
    EXPECT_LT(results.final.kernel_sum_interior.max_value, 0.1);
    EXPECT_GT(results.final.nearest_min, 0.8);
    EXPECT_GT(results.final.neighbors_min, 10.0);
}
TEST(Relaxation, ResidualsAreReported)
{
    // R_B is only reported: the standard relaxation reduces R_P, not R_B
    EXPECT_TRUE(std::isfinite(results.final.residual_b_interior.rms()));
    EXPECT_TRUE(std::isfinite(results.final.residual_b_interior.max_value));
    EXPECT_GT(results.final.residual_b_interior.count, size_t(0));
}
TEST(ConstantField, GradientIsZero)
{
    if (!std::is_same<Real, double>::value)
        return;
    EXPECT_GT(results.constant.core.linear_gradient.count, size_t(0));
    EXPECT_LT(results.constant.interior.linear_gradient.max_value, 1.0e-8);
    EXPECT_LT(results.constant.interior.second_order_gradient.max_value, 1.0e-8);
}
TEST(LinearField, GradientIsExact)
{
    if (!std::is_same<Real, double>::value)
        return;
    EXPECT_LT(results.linear.interior.linear_gradient.max_value, 1.0e-6);
    EXPECT_LT(results.linear.interior.second_order_gradient.max_value, 1.0e-6);
}
TEST(QuadraticField, SecondOrderGradientIsExactLinearGradientIsNot)
{
    if (!std::is_same<Real, double>::value)
        return;
    const OperatorErrors &core = results.quadratic.core;
    EXPECT_GT(core.linear_gradient.count, size_t(0));
    EXPECT_LT(core.second_order_gradient.max_value, 1.0e-6);
    EXPECT_LT(results.quadratic.interior.second_order_gradient.max_value, 1.0e-6);
    // LinearGradient is only first-order consistent: a small, but clearly non-zero error
    EXPECT_GT(core.linear_gradient.rms(), 1.0e-5);
    EXPECT_LT(core.linear_gradient.rms(), 1.0e-2);
    EXPECT_GT(core.linear_gradient.rms(), 100.0 * core.second_order_gradient.rms());
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
template <typename DataType>
StdVec<DataType> copyParticleVariable(BaseParticles &particles, const std::string &name, size_t total_particles)
{
    DataType *data = particles.getVariableDataByName<DataType>(name);
    return StdVec<DataType>(data, data + total_particles);
}
//----------------------------------------------------------------------
//	Gradient errors of one field at the current particle positions.
//	SecondOrderGradient overwrites the variable written by LinearGradient, so each result is
//	copied right after its operator has been executed.
//----------------------------------------------------------------------
template <class GradientReference>
FieldErrors computeFieldErrors(Inner<> &inner_relation, BaseParticles &particles, const std::string &name,
                               GradientReference reference_gradient, Real support_radius)
{
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        linear_gradient(DynamicsArgs(inner_relation, name));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        hessian(DynamicsArgs(inner_relation, name));
    InteractionDynamicsCK<MainExecutionPolicy, SecondOrderGradient<Inner<Real>>>
        second_order_gradient(DynamicsArgs(inner_relation, name));

    const size_t total_particles = particles.TotalRealParticles();
    linear_gradient.exec();
    const StdVec<Vecd> linear_result = copyParticleVariable<Vecd>(particles, name + "Gradient", total_particles);
    hessian.exec();
    second_order_gradient.exec();
    const StdVec<Vecd> second_order_result = copyParticleVariable<Vecd>(particles, name + "Gradient", total_particles);

    const Vecd *position = particles.getVariableDataByName<Vecd>("Position");
    FieldErrors errors;
    for (size_t i = 0; i != total_particles; ++i)
    {
        const Vecd reference = reference_gradient(position[i]);
        const Real linear_error = (linear_result[i] - reference).norm();
        const Real second_order_error = (second_order_result[i] - reference).norm();
        if (circle_radius - position[i].norm() >= support_radius)
        {
            errors.interior.linear_gradient.add(linear_error);
            errors.interior.second_order_gradient.add(second_order_error);
        }
        if (position[i].norm() <= core_radius)
        {
            errors.core.linear_gradient.add(linear_error);
            errors.core.second_order_gradient.add(second_order_error);
        }
    }
    return errors;
}
//----------------------------------------------------------------------
//	CSV files for plotting. They only write what has been computed above.
//----------------------------------------------------------------------
void writeHistoryCsv(const std::string &file_name)
{
    std::ofstream file(file_name);
    file << "step,interior_n,band_n,hR_P_interior_rms,hR_B_interior_rms,hR_P_band_rms,hR_B_band_rms\n";
    char line[256];
    for (const HistoryRow &row : history)
    {
        const ConsistencyMetrics &m = row.metrics;
        std::snprintf(line, sizeof(line), "%d,%zu,%zu,%.9e,%.9e,%.9e,%.9e\n", row.step,
                      m.residual_p_interior.count, m.residual_p_band.count,
                      double(m.residual_p_interior.rms()), double(m.residual_b_interior.rms()),
                      double(m.residual_p_band.rms()), double(m.residual_b_band.rms()));
        file << line;
    }
}
void writeGradientCsv(const std::string &file_name)
{
    std::ofstream file(file_name);
    file << "field,region,n,linear_rms,linear_max,second_order_rms,second_order_max\n";
    const char *field_names[3] = {"constant", "linear", "quadratic"};
    const FieldErrors *fields[3] = {&results.constant, &results.linear, &results.quadratic};
    const char *region_names[2] = {"core", "interior"};
    char line[256];
    for (int f = 0; f != 3; ++f)
    {
        const OperatorErrors *regions[2] = {&fields[f]->core, &fields[f]->interior};
        for (int k = 0; k != 2; ++k)
        {
            const OperatorErrors &r = *regions[k];
            std::snprintf(line, sizeof(line), "%s,%s,%zu,%.9e,%.9e,%.9e,%.9e\n", field_names[f],
                          region_names[k], r.linear_gradient.count, double(r.linear_gradient.rms()),
                          double(r.linear_gradient.max_value), double(r.second_order_gradient.rms()),
                          double(r.second_order_gradient.max_value));
            file << line;
        }
    }
}
// interior = 1 for the particles at least the support radius away from the boundary
void writeParticlesCsv(const std::string &file_name, const StdVec<Vecd> &initial_position,
                       const Vecd *final_position, Real support_radius)
{
    std::ofstream file(file_name);
    file << "index,x_initial,y_initial,x_final,y_final,interior\n";
    char line[256];
    for (size_t i = 0; i != initial_position.size(); ++i)
    {
        std::snprintf(line, sizeof(line), "%zu,%.12e,%.12e,%.12e,%.12e,%d\n", i,
                      double(initial_position[i][0]), double(initial_position[i][1]),
                      double(final_position[i][0]), double(final_position[i][1]),
                      circle_radius - final_position[i].norm() >= support_radius ? 1 : 0);
        file << line;
    }
}
//----------------------------------------------------------------------
//	The main program.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    particle_spacing = readEnvironmentReal("RELAX_DP", particle_spacing, 0.005, 0.2);
    target_residual = readEnvironmentReal("RELAX_TARGET", target_residual, 1.0e-6, 1.0);
    max_relaxation_steps = int(readEnvironmentReal("RELAX_MAX_STEPS", Real(max_relaxation_steps), 1, 100000));
    //----------------------------------------------------------------------
    //	Build up the SPHSystem, the body and the particles.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-1.5 * circle_radius, -1.5 * circle_radius),
                                      Vecd(1.5 * circle_radius, 1.5 * circle_radius));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av);

    RealBody circle(sph_system, makeShared<CircleBody>("Circle"));
    circle.defineBodyLevelSetShape(2.0);
    circle.generateParticles<BaseParticles, Lattice>();
    BaseParticles &particles = circle.getBaseParticles();
    particles.registerStateVariable<Real>("Constant");
    particles.registerStateVariable<Real>("Linear");
    particles.registerStateVariable<Real>("Quadratic");
    //----------------------------------------------------------------------
    //	Relations and the dynamics. The relaxation and the residual of the kernel gradient use the
    //	standard relation, the correction matrix and the gradient operators use the relation of
    //	the CK dynamics. Both are built from the same particle positions and kernel.
    //----------------------------------------------------------------------
    InnerRelation circle_inner(circle);
    Inner<> circle_inner_ck(circle);
    using namespace relax_dynamics;
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(circle_inner);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> cell_linked_list_ck(circle);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_ck(circle_inner_ck);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        linear_correction_matrix(DynamicsArgs(circle_inner_ck, 0.0));
    InteractionDynamics<ZeroOrderConsistency> consistency(circle_inner);
    InteractionDynamics<CorrectionMatrixResidual> correction_residual(circle_inner);
    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>>>
        displacement_matrix_gradient(circle_inner_ck);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>>>
        hessian_correction_matrix(DynamicsArgs(circle_inner_ck, 0.0));

    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ConstantProfile>>>
        constant_assignment(circle, "Constant");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<LinearProfile>>>
        linear_assignment(circle, "Linear");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ParabolicProfile>>>
        quadratic_assignment(circle, "Quadratic");
    //----------------------------------------------------------------------
    //	Random shift of the lattice with a seeded generator.
    //----------------------------------------------------------------------
    const size_t total_particles = particles.TotalRealParticles();
    Vecd *position = particles.getVariableDataByName<Vecd>("Position");
    applySeededJitter(position, total_particles, jitter_ratio * particle_spacing, jitter_seed);
    relaxation_step_inner.SurfaceBounding().exec();
    circle.updateCellLinkedList();
    const StdVec<Vecd> initial_position(position, position + total_particles);
    //----------------------------------------------------------------------
    //	Quantities needed by the diagnostics.
    //----------------------------------------------------------------------
    const Vecd *kernel_gradient_sum = particles.getVariableDataByName<Vecd>("KernelGradientSum");
    const Real *kernel_sum = particles.getVariableDataByName<Real>("KernelSum");
    const Real *neighbor_number = particles.getVariableDataByName<Real>("NeighborNumber");
    const Real *nearest_distance = particles.getVariableDataByName<Real>("NearestNeighborDistance");
    const Vecd *correction_matrix_residual = particles.getVariableDataByName<Vecd>("CorrectionMatrixResidual");
    const Real smoothing_length = circle.getSPHAdaptation().ReferenceSmoothingLength();
    const Real support_radius = circle.getSPHAdaptation().getKernel()->CutOffRadius();
    char line[512];
    std::snprintf(line, sizeof(line),
                  "[setup] dp=%g seed=%u jitter=%g target_hR_P_rms=%.1e max_steps=%d check_interval=%d "
                  "particles=%zu h=%g support_radius=%g core_radius=%g precision=%s",
                  double(particle_spacing), jitter_seed, double(jitter_ratio), double(target_residual),
                  max_relaxation_steps, check_interval, total_particles, double(smoothing_length),
                  double(support_radius), double(core_radius), std::is_same<Real, double>::value ? "double" : "float");
    emit(line);

    // Evaluates the diagnostics at the current particle positions.
    auto measure = [&]() -> ConsistencyMetrics
    {
        circle.updateCellLinkedList();
        circle_inner.updateConfiguration();
        cell_linked_list_ck.exec();
        update_inner_ck.exec();
        linear_correction_matrix.exec();
        consistency.exec();
        correction_residual.exec();
        ConsistencyMetrics metrics;
        metrics.nearest_min = MaxReal;
        metrics.neighbors_min = MaxReal;
        for (size_t i = 0; i != total_particles; ++i)
        {
            const Real residual_p = smoothing_length * kernel_gradient_sum[i].norm();
            const Real residual_b = smoothing_length * correction_matrix_residual[i].norm();
            if (circle_radius - position[i].norm() >= support_radius)
            {
                metrics.residual_p_interior.add(residual_p);
                metrics.residual_b_interior.add(residual_b);
                metrics.kernel_sum_interior.add(std::abs(kernel_sum[i] - 1.0));
                metrics.nearest_min = SMIN(metrics.nearest_min, nearest_distance[i] / particle_spacing);
                metrics.neighbors_min = SMIN(metrics.neighbors_min, neighbor_number[i]);
            }
            else
            {
                metrics.residual_p_band.add(residual_p);
                metrics.residual_b_band.add(residual_b);
            }
        }
        return metrics;
    };
    //----------------------------------------------------------------------
    //	Relaxation until the target residual is reached or the step limit.
    //----------------------------------------------------------------------
    results.initial = measure();
    printMetrics("relax", 0, results.initial);
    history.push_back({0, results.initial});
    results.final = results.initial;
    int step = 0;
    while (step != max_relaxation_steps && !results.converged)
    {
        relaxation_step_inner.exec();
        ++step;
        if (step % check_interval == 0 || step == max_relaxation_steps)
        {
            results.final = measure();
            printMetrics("relax", step, results.final);
            history.push_back({step, results.final});
            results.converged = results.final.residual_p_interior.rms() < target_residual;
        }
    }
    results.steps = step;
    std::snprintf(line, sizeof(line), "[result] relaxation_steps=%d converged=%d final_hR_P_rms=%.3e initial_hR_P_rms=%.3e",
                  results.steps, results.converged ? 1 : 0, double(results.final.residual_p_interior.rms()),
                  double(results.initial.residual_p_interior.rms()));
    emit(line);
    printMetrics("result", results.steps, results.final);
    //----------------------------------------------------------------------
    //	Gradient accuracy on the relaxed particles. The correction matrices are computed at the
    //	final positions by the last call of measure().
    //----------------------------------------------------------------------
    displacement_matrix_gradient.exec();
    hessian_correction_matrix.exec();
    constant_assignment.exec();
    linear_assignment.exec();
    quadratic_assignment.exec();
    results.constant = computeFieldErrors(circle_inner_ck, particles, "Constant", constantGradient, support_radius);
    results.linear = computeFieldErrors(circle_inner_ck, particles, "Linear", linearGradient, support_radius);
    results.quadratic = computeFieldErrors(circle_inner_ck, particles, "Quadratic", quadraticGradient, support_radius);
    printFieldErrors("constant", results.constant);
    printFieldErrors("linear", results.linear);
    printFieldErrors("quadratic", results.quadratic);

    char file_stem[128];
    std::snprintf(file_stem, sizeof(file_stem), "relaxation_consistency_dp%g", double(particle_spacing));
    std::ofstream file(std::string(file_stem) + ".txt");
    file << output_text;
    writeHistoryCsv(std::string(file_stem) + "_history.csv");
    writeGradientCsv(std::string(file_stem) + "_gradient.csv");
    writeParticlesCsv(std::string(file_stem) + "_particles.csv", initial_position, position, support_radius);
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
