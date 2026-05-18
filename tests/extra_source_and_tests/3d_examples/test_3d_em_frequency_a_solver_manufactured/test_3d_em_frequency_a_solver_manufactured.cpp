/**
 * @file test_3d_em_frequency_a_solver_manufactured.cpp
 * @brief Manufactured-solution verification for frequency-domain A-component solver bodies.
 *
 * This case isolates the missing link after the operator smoke tests:
 *
 *   solve A_real with exact A_imag and exact grad(phi_real) frozen,
 *   solve A_imag with exact A_real and exact grad(phi_imag) frozen,
 *
 * using the actual frequency-domain inner A-equation updater
 * `VectorPotentialFrequencyEquationInner`.
 *
 * Analytical fields:
 *   A_real = 0.5 * (B_real x r)
 *   A_imag = 0.5 * (B_imag x r)
 *   grad(phi_real), grad(phi_imag) = prescribed constants
 *   curl(nu curl A_real/im)) = 0  (constant reluctivity, linear A)
 *
 * Manufactured source currents are chosen so that the continuous frequency-domain
 * conductor equations are exactly satisfied:
 *   Js_real = sigma * grad(phi_real) - sigma * omega * A_imag
 *   Js_imag = sigma * grad(phi_imag) + sigma * omega * A_real
 *
 * The case therefore checks whether the numerical A-component solver converges to
 * the analytical linear A field under fixed auxiliary fields.
 */

#include "sphinxsys.h"
#include "aphi_case_support/electromagnetic_team7_aphi_dynamics.hpp"
#include "aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

using namespace SPH;

namespace
{
Real get_env_real_local(const std::string &name, Real default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    char *end_ptr = nullptr;
    Real parsed = static_cast<Real>(std::strtod(value, &end_ptr));
    if (end_ptr == value || !std::isfinite(parsed))
    {
        return default_value;
    }
    return parsed;
}

bool get_env_bool_local(const std::string &name, bool default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    std::string token(value);
    std::transform(token.begin(), token.end(), token.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (token == "1" || token == "true" || token == "on" || token == "yes")
    {
        return true;
    }
    if (token == "0" || token == "false" || token == "off" || token == "no")
    {
        return false;
    }
    return default_value;
}

std::string get_env_string_local(const std::string &name,
                                 const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

const Real dp_0 = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_DP", 1.0);
const Real body_length = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_LENGTH", 20.0);
const Real body_height = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BOUNDARY_WIDTH", 3.0 * dp_0);
const Real boundary_shell_thickness =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BOUNDARY_SHELL_THICKNESS", 2.5 * dp_0);
const Real conductivity = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_SIGMA", 3.0e6);
const Real rho_cp = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_RHO_CP", 1.0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Real frequency_hz = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_FREQUENCY_HZ", 50.0);
const Real omega = 2.0 * Pi * frequency_hz;
const Real dt_pseudo = get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_DT", 1.0e-6);
const Real reference_dt_pseudo =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_REFERENCE_DT", dt_pseudo);
const int solve_iterations =
    static_cast<int>(get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_ITERATIONS", 200.0));
const Real sigma_relaxation_scaling =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_SIGMA_RELAX_SCALE", 1.0);
const Real sigma_relaxation_floor =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_SIGMA_RELAX_FLOOR", TinyReal);
const Real magnetic_diagonal_scaling =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_MAG_DIAG_SCALE", 1.0);
const Real relaxation_scaling =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_RELAXATION_SCALING", 1.0);
const Real max_change_rate =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_MAX_CHANGE_RATE", 1.0e6);
const Real curl_scaling =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_CURL_SCALING", 1.0);
const Real curl_nu_b_scaling =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_CURL_NUB_SCALING", 1.0);
const Real initial_guess_scale_real =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_INITIAL_GUESS_SCALE_REAL", 0.0);
const Real initial_guess_scale_imag =
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_INITIAL_GUESS_SCALE_IMAG", 0.0);
const bool write_vtp = get_env_bool_local("EM_FREQ_A_SOLVER_VERIFY_WRITE_VTP", true);
const bool write_particles = get_env_bool_local("EM_FREQ_A_SOLVER_VERIFY_WRITE_PARTICLES", true);
const bool use_boundary_constraint =
    get_env_bool_local("EM_FREQ_A_SOLVER_VERIFY_USE_BOUNDARY_CONSTRAINT", true);
const std::string solver_mode_raw =
    get_env_string_local("EM_FREQ_A_SOLVER_VERIFY_SOLVER_MODE", "single");
const int history_interval =
    static_cast<int>(get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_HISTORY_INTERVAL", 50.0));

const Vec3d target_b_real(get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BREAL_X", 0.0),
                          get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BREAL_Y", 0.0),
                          get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BREAL_Z", 1.0));
const Vec3d target_b_imag(get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BIMAG_X", 0.0),
                          get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BIMAG_Y", 0.0),
                          get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_BIMAG_Z", 0.5));
const Vec3d target_grad_phi_real(
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_REAL_X", 0.4),
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_REAL_Y", -0.2),
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_REAL_Z", 0.1));
const Vec3d target_grad_phi_imag(
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_IMAG_X", -0.3),
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_IMAG_Y", 0.1),
    get_env_real_local("EM_FREQ_A_SOLVER_VERIFY_GRADPHI_IMAG_Z", 0.2));

const Vec3d body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d gauge_origin = body_center;
BoundingBoxd system_domain_bounds(
    Vec3d(-boundary_width, -boundary_width, -boundary_width),
    Vec3d(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

class ConductorShape : public ComplexShape
{
  public:
    explicit ConductorShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(body_center), body_halfsize);
    }
};

class OuterBoundaryShellShape : public ComplexShape
{
  public:
    explicit OuterBoundaryShellShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d inner_halfsize = body_halfsize - Vec3d::Ones() * boundary_shell_thickness;
        add<GeometricShapeBox>(Transform(body_center), body_halfsize, "OuterBody");
        if (inner_halfsize.minCoeff() > TinyReal)
        {
            subtract<GeometricShapeBox>(Transform(body_center), inner_halfsize, "InnerCore");
        }
    }
};

struct LinearVectorPotentialField
{
    Vec3d target_b;
    Vec3d origin;

    Vec3d evaluate(const Vec3d &position) const
    {
        Vec3d relative_position = position - origin;
        return 0.5 * target_b.cross(relative_position);
    }
};

const LinearVectorPotentialField a_real_model{target_b_real, gauge_origin};
const LinearVectorPotentialField a_imag_model{target_b_imag, gauge_origin};

std::string to_lower_copy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

bool solver_mode_contains(const std::string &solver_mode, const std::string &token)
{
    return solver_mode == token || solver_mode == "both";
}

class AssignManufacturedAuxiliaryFields : public LocalDynamics
{
  public:
    explicit AssignManufacturedAuxiliaryFields(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          a_real_(particles_->getVariableDataByName<Vecd>("VectorPotentialReal")),
          a_imag_(particles_->getVariableDataByName<Vecd>("VectorPotentialImag")),
          grad_phi_real_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientReal")),
          grad_phi_imag_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientImag")),
          source_current_density_real_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityReal")),
          source_current_density_imag_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityImag")),
          curl_real_(particles_->getVariableDataByName<AngularVecd>("VectorPotentialCurlReal")),
          curl_imag_(particles_->getVariableDataByName<AngularVecd>("VectorPotentialCurlImag")),
          curl_nu_b_real_(particles_->getVariableDataByName<Vecd>("CurlNuBReal")),
          curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>("CurlNuBImag")),
          change_rate_real_(particles_->getVariableDataByName<Vecd>("VectorPotentialChangeRateReal")),
          change_rate_imag_(particles_->getVariableDataByName<Vecd>("VectorPotentialChangeRateImag")),
          conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        const Vec3d position = positions_[index_i];
        const Vec3d a_real_exact = a_real_model.evaluate(position);
        const Vec3d a_imag_exact = a_imag_model.evaluate(position);
        const Real sigma_i = conductivity_[index_i];

        a_real_[index_i] = a_real_exact;
        a_imag_[index_i] = a_imag_exact;
        grad_phi_real_[index_i] = target_grad_phi_real;
        grad_phi_imag_[index_i] = target_grad_phi_imag;
        source_current_density_real_[index_i] =
            sigma_i * target_grad_phi_real - sigma_i * omega * a_imag_exact;
        source_current_density_imag_[index_i] =
            sigma_i * target_grad_phi_imag + sigma_i * omega * a_real_exact;
        curl_real_[index_i] = ZeroData<AngularVecd>::value;
        curl_imag_[index_i] = ZeroData<AngularVecd>::value;
        curl_nu_b_real_[index_i] = ZeroData<Vecd>::value;
        curl_nu_b_imag_[index_i] = ZeroData<Vecd>::value;
        change_rate_real_[index_i] = ZeroData<Vecd>::value;
        change_rate_imag_[index_i] = ZeroData<Vecd>::value;
    }

  protected:
    Vecd *positions_;
    Vecd *a_real_, *a_imag_;
    Vecd *grad_phi_real_, *grad_phi_imag_;
    Vecd *source_current_density_real_, *source_current_density_imag_;
    AngularVecd *curl_real_, *curl_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    Vecd *change_rate_real_, *change_rate_imag_;
    Real *conductivity_;
};

class SetManufacturedInitialGuess : public LocalDynamics
{
  public:
    explicit SetManufacturedInitialGuess(SPHBody &sph_body,
                                         const std::string &field_name,
                                         const LinearVectorPotentialField &model,
                                         Real scale)
        : LocalDynamics(sph_body),
          field_name_(field_name),
          model_(model),
          scale_(scale),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          vector_field_(particles_->getVariableDataByName<Vecd>(field_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        vector_field_[index_i] = scale_ * model_.evaluate(positions_[index_i]);
    }

  protected:
    std::string field_name_;
    LinearVectorPotentialField model_;
    Real scale_;
    Vecd *positions_;
    Vecd *vector_field_;
};

class ConstrainManufacturedVectorPotential : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainManufacturedVectorPotential(BodyPartByParticle &body_part,
                                                  const std::string &field_name,
                                                  const LinearVectorPotentialField &model)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          model_(model),
          positions_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          vector_field_(this->particles_->template getVariableDataByName<Vecd>(field_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        vector_field_[index_i] = model_.evaluate(positions_[index_i]);
    }

  protected:
    LinearVectorPotentialField model_;
    Vecd *positions_;
    Vecd *vector_field_;
};

Vec3d angular_to_vec(const AngularVecd &value)
{
    return Vec3d(value[0], value[1], value[2]);
}

struct ASolverSummary
{
    std::string component_name;
    std::string solver_mode;
    std::string solve_phase;
    size_t total_particles = 0;
    int iterations = 0;
    Real initial_guess_scale = 0.0;
    Real mean_a_error = 0.0;
    Real max_a_error = 0.0;
    Real mean_b_error = 0.0;
    Real max_b_error = 0.0;
    Real mean_curl_nu_b_norm = 0.0;
    Real max_curl_nu_b_norm = 0.0;
    Real mean_residual_norm = 0.0;
    Real max_residual_norm = 0.0;
    Real mean_relative_residual = 0.0;
    Real max_relative_residual = 0.0;
    Real mean_change_rate_norm = 0.0;
    Real max_change_rate_norm = 0.0;
};

struct ComponentIterationMetrics
{
    Real mean_curl_nu_b_norm = 0.0;
    Real max_curl_nu_b_norm = 0.0;
    Real mean_residual_norm = 0.0;
    Real max_residual_norm = 0.0;
    Real mean_relative_residual = 0.0;
    Real max_relative_residual = 0.0;
    Real mean_change_rate_norm = 0.0;
    Real max_change_rate_norm = 0.0;
};

ComponentIterationMetrics evaluate_real_iteration_metrics(BaseParticles &particles)
{
    ComponentIterationMetrics metrics;
    Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
    Vecd *grad_phi_real = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
    Vecd *source_current_density_real = particles.getVariableDataByName<Vecd>("SourceCurrentDensityReal");
    Vecd *curl_nu_b_real = particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *change_rate_real = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
    Real *relative_residual_real =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualReal");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_curl_nu_b_norm = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d residual = source_current_density_real[i] - curl_nu_b_real[i] -
                         sigma[i] * grad_phi_real[i] + sigma[i] * omega * a_imag[i];
        Real curl_nu_b_norm = curl_nu_b_real[i].norm();
        Real residual_norm = residual.norm();
        Real relative_residual = relative_residual_real[i];
        Real change_rate_norm = change_rate_real[i].norm();
        sum_curl_nu_b_norm += curl_nu_b_norm;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual;
        sum_change_rate_norm += change_rate_norm;
        metrics.max_curl_nu_b_norm = SMAX(metrics.max_curl_nu_b_norm, curl_nu_b_norm);
        metrics.max_residual_norm = SMAX(metrics.max_residual_norm, residual_norm);
        metrics.max_relative_residual = SMAX(metrics.max_relative_residual, relative_residual);
        metrics.max_change_rate_norm = SMAX(metrics.max_change_rate_norm, change_rate_norm);
    }

    metrics.mean_curl_nu_b_norm =
        sum_curl_nu_b_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_residual_norm =
        sum_residual_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_relative_residual =
        sum_relative_residual / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_change_rate_norm =
        sum_change_rate_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    return metrics;
}

ComponentIterationMetrics evaluate_imag_iteration_metrics(BaseParticles &particles)
{
    ComponentIterationMetrics metrics;
    Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
    Vecd *grad_phi_imag = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
    Vecd *source_current_density_imag = particles.getVariableDataByName<Vecd>("SourceCurrentDensityImag");
    Vecd *curl_nu_b_imag = particles.getVariableDataByName<Vecd>("CurlNuBImag");
    Vecd *change_rate_imag = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
    Real *relative_residual_imag =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualImag");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_curl_nu_b_norm = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d residual = source_current_density_imag[i] - curl_nu_b_imag[i] -
                         sigma[i] * grad_phi_imag[i] - sigma[i] * omega * a_real[i];
        Real curl_nu_b_norm = curl_nu_b_imag[i].norm();
        Real residual_norm = residual.norm();
        Real relative_residual = relative_residual_imag[i];
        Real change_rate_norm = change_rate_imag[i].norm();
        sum_curl_nu_b_norm += curl_nu_b_norm;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual;
        sum_change_rate_norm += change_rate_norm;
        metrics.max_curl_nu_b_norm = SMAX(metrics.max_curl_nu_b_norm, curl_nu_b_norm);
        metrics.max_residual_norm = SMAX(metrics.max_residual_norm, residual_norm);
        metrics.max_relative_residual = SMAX(metrics.max_relative_residual, relative_residual);
        metrics.max_change_rate_norm = SMAX(metrics.max_change_rate_norm, change_rate_norm);
    }

    metrics.mean_curl_nu_b_norm =
        sum_curl_nu_b_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_residual_norm =
        sum_residual_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_relative_residual =
        sum_relative_residual / (static_cast<Real>(total_real_particles) + TinyReal);
    metrics.mean_change_rate_norm =
        sum_change_rate_norm / (static_cast<Real>(total_real_particles) + TinyReal);
    return metrics;
}

ASolverSummary evaluate_real_component_summary(BaseParticles &particles)
{
    ASolverSummary summary;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
    Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
    Vecd *grad_phi_real = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
    Vecd *source_current_density_real = particles.getVariableDataByName<Vecd>("SourceCurrentDensityReal");
    AngularVecd *curl_real = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
    Vecd *curl_nu_b_real = particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *change_rate_real = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
    Real *relative_residual_real =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualReal");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_a_error = 0.0;
    Real sum_b_error = 0.0;
    Real sum_curl_nu_b_norm = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = a_real_model.evaluate(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_real[i]);
        Vec3d residual = source_current_density_real[i] - curl_nu_b_real[i] -
                         sigma[i] * grad_phi_real[i] + sigma[i] * omega * a_imag[i];
        Real a_error = (a_real[i] - a_exact).norm();
        Real b_error = (b_discrete - target_b_real).norm();
        Real curl_nu_b_norm = curl_nu_b_real[i].norm();
        Real residual_norm = residual.norm();
        Real relative_residual = relative_residual_real[i];
        Real change_rate_norm = change_rate_real[i].norm();
        summary.total_particles++;
        sum_a_error += a_error;
        sum_b_error += b_error;
        sum_curl_nu_b_norm += curl_nu_b_norm;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual;
        sum_change_rate_norm += change_rate_norm;
        summary.max_a_error = SMAX(summary.max_a_error, a_error);
        summary.max_b_error = SMAX(summary.max_b_error, b_error);
        summary.max_curl_nu_b_norm = SMAX(summary.max_curl_nu_b_norm, curl_nu_b_norm);
        summary.max_residual_norm = SMAX(summary.max_residual_norm, residual_norm);
        summary.max_relative_residual = SMAX(summary.max_relative_residual, relative_residual);
        summary.max_change_rate_norm = SMAX(summary.max_change_rate_norm, change_rate_norm);
    }

    summary.mean_a_error = sum_a_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_b_error = sum_b_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_curl_nu_b_norm =
        sum_curl_nu_b_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_residual_norm =
        sum_residual_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_relative_residual =
        sum_relative_residual / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_change_rate_norm =
        sum_change_rate_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    return summary;
}

ASolverSummary evaluate_imag_component_summary(BaseParticles &particles)
{
    ASolverSummary summary;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
    Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
    Vecd *grad_phi_imag = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
    Vecd *source_current_density_imag = particles.getVariableDataByName<Vecd>("SourceCurrentDensityImag");
    AngularVecd *curl_imag = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
    Vecd *curl_nu_b_imag = particles.getVariableDataByName<Vecd>("CurlNuBImag");
    Vecd *change_rate_imag = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
    Real *relative_residual_imag =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualImag");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_a_error = 0.0;
    Real sum_b_error = 0.0;
    Real sum_curl_nu_b_norm = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = a_imag_model.evaluate(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_imag[i]);
        Vec3d residual = source_current_density_imag[i] - curl_nu_b_imag[i] -
                         sigma[i] * grad_phi_imag[i] - sigma[i] * omega * a_real[i];
        Real a_error = (a_imag[i] - a_exact).norm();
        Real b_error = (b_discrete - target_b_imag).norm();
        Real curl_nu_b_norm = curl_nu_b_imag[i].norm();
        Real residual_norm = residual.norm();
        Real relative_residual = relative_residual_imag[i];
        Real change_rate_norm = change_rate_imag[i].norm();
        summary.total_particles++;
        sum_a_error += a_error;
        sum_b_error += b_error;
        sum_curl_nu_b_norm += curl_nu_b_norm;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual;
        sum_change_rate_norm += change_rate_norm;
        summary.max_a_error = SMAX(summary.max_a_error, a_error);
        summary.max_b_error = SMAX(summary.max_b_error, b_error);
        summary.max_curl_nu_b_norm = SMAX(summary.max_curl_nu_b_norm, curl_nu_b_norm);
        summary.max_residual_norm = SMAX(summary.max_residual_norm, residual_norm);
        summary.max_relative_residual = SMAX(summary.max_relative_residual, relative_residual);
        summary.max_change_rate_norm = SMAX(summary.max_change_rate_norm, change_rate_norm);
    }

    summary.mean_a_error = sum_a_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_b_error = sum_b_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_curl_nu_b_norm =
        sum_curl_nu_b_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_residual_norm =
        sum_residual_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_relative_residual =
        sum_relative_residual / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_change_rate_norm =
        sum_change_rate_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    return summary;
}

void write_particle_diagnostics_real(const std::string &file_path, BaseParticles &particles,
                                     const std::string &component_label)
{
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    file << std::setprecision(12);
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
    Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
    Vecd *grad_phi_real = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
    Vecd *source_current_density_real = particles.getVariableDataByName<Vecd>("SourceCurrentDensityReal");
    AngularVecd *curl_real = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
    Vecd *curl_nu_b_real = particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *change_rate_real = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
    Real *relative_residual_real =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualReal");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");
    size_t total_real_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = a_real_model.evaluate(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_real[i]);
        Vec3d residual = source_current_density_real[i] - curl_nu_b_real[i] -
                         sigma[i] * grad_phi_real[i] + sigma[i] * omega * a_imag[i];
        file << component_label << ","
             << i << ","
             << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << ","
             << a_real[i][0] << "," << a_real[i][1] << "," << a_real[i][2] << ","
             << a_exact[0] << "," << a_exact[1] << "," << a_exact[2] << ","
             << (a_real[i] - a_exact).norm() << ","
             << b_discrete[0] << "," << b_discrete[1] << "," << b_discrete[2] << ","
             << (b_discrete - target_b_real).norm() << ","
             << curl_nu_b_real[i][0] << "," << curl_nu_b_real[i][1] << "," << curl_nu_b_real[i][2] << ","
             << curl_nu_b_real[i].norm() << ","
             << residual[0] << "," << residual[1] << "," << residual[2] << ","
             << residual.norm() << ","
             << relative_residual_real[i] << ","
             << change_rate_real[i][0] << "," << change_rate_real[i][1] << "," << change_rate_real[i][2] << ","
             << change_rate_real[i].norm() << "\n";
    }
}

void write_particle_diagnostics_imag(const std::string &file_path, BaseParticles &particles,
                                     const std::string &component_label)
{
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    file << std::setprecision(12);
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
    Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
    Vecd *grad_phi_imag = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
    Vecd *source_current_density_imag = particles.getVariableDataByName<Vecd>("SourceCurrentDensityImag");
    AngularVecd *curl_imag = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
    Vecd *curl_nu_b_imag = particles.getVariableDataByName<Vecd>("CurlNuBImag");
    Vecd *change_rate_imag = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
    Real *relative_residual_imag =
        particles.getVariableDataByName<Real>("AEquationRelativeResidualImag");
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");
    size_t total_real_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = a_imag_model.evaluate(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_imag[i]);
        Vec3d residual = source_current_density_imag[i] - curl_nu_b_imag[i] -
                         sigma[i] * grad_phi_imag[i] - sigma[i] * omega * a_real[i];
        file << component_label << ","
             << i << ","
             << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << ","
             << a_imag[i][0] << "," << a_imag[i][1] << "," << a_imag[i][2] << ","
             << a_exact[0] << "," << a_exact[1] << "," << a_exact[2] << ","
             << (a_imag[i] - a_exact).norm() << ","
             << b_discrete[0] << "," << b_discrete[1] << "," << b_discrete[2] << ","
             << (b_discrete - target_b_imag).norm() << ","
             << curl_nu_b_imag[i][0] << "," << curl_nu_b_imag[i][1] << "," << curl_nu_b_imag[i][2] << ","
             << curl_nu_b_imag[i].norm() << ","
             << residual[0] << "," << residual[1] << "," << residual[2] << ","
             << residual.norm() << ","
             << relative_residual_imag[i] << ","
             << change_rate_imag[i][0] << "," << change_rate_imag[i][1] << "," << change_rate_imag[i][2] << ","
             << change_rate_imag[i].norm() << "\n";
    }
}

} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_FREQ_A_SOLVER_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::string solver_mode = to_lower_copy(solver_mode_raw);
    if (solver_mode != "single" && solver_mode != "coupled" && solver_mode != "both")
    {
        std::cout << "[em-freq-a-solver-verify-warning] unsupported solver mode '"
                  << solver_mode_raw << "', fallback to 'single'." << std::endl;
        solver_mode = "single";
    }
    const int effective_history_interval = SMAX(1, history_interval);
    const bool run_single = solver_mode_contains(solver_mode, "single");
    const bool run_coupled = solver_mode_contains(solver_mode, "coupled");

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-a-solver-verify-config] dp=" << dp_0
              << ", sigma=" << conductivity
              << ", frequency_hz=" << frequency_hz
              << ", omega=" << omega
              << ", iterations=" << solve_iterations
              << ", dt=" << dt_pseudo
              << ", reference_dt=" << reference_dt_pseudo
              << ", initial_guess_scale_real=" << initial_guess_scale_real
              << ", initial_guess_scale_imag=" << initial_guess_scale_imag
              << ", solver_mode=" << solver_mode
              << ", history_interval=" << effective_history_interval
              << ", use_boundary_constraint=" << use_boundary_constraint
              << ", B_real=(" << target_b_real[0] << "," << target_b_real[1] << "," << target_b_real[2] << ")"
              << ", B_imag=(" << target_b_imag[0] << "," << target_b_imag[1] << "," << target_b_imag[2] << ")"
              << std::fixed << std::setprecision(6) << std::endl;

    SolidBody conductor_body(sph_system, makeShared<ConductorShape>("Conductor"));
    conductor_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    conductor_body.defineMaterial<Solid>();
    conductor_body.defineBodyLevelSetShape();
    conductor_body.generateParticles<BaseParticles, Lattice>();

    BodyRegionByParticle boundary_region(
        conductor_body,
        makeShared<OuterBoundaryShellShape>("ConductorBoundaryShell"));
    InnerRelation conductor_inner(conductor_body);

    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_frequency_em(conductor_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<AssignManufacturedAuxiliaryFields>
        assign_manufactured_auxiliary_fields(conductor_body);
    SimpleDynamics<SetManufacturedInitialGuess>
        set_initial_guess_real(conductor_body, "VectorPotentialReal", a_real_model, initial_guess_scale_real);
    SimpleDynamics<SetManufacturedInitialGuess>
        set_initial_guess_imag(conductor_body, "VectorPotentialImag", a_imag_model, initial_guess_scale_imag);
    SimpleDynamics<ConstrainManufacturedVectorPotential>
        constrain_a_real_boundary(boundary_region, "VectorPotentialReal", a_real_model);
    SimpleDynamics<ConstrainManufacturedVectorPotential>
        constrain_a_imag_boundary(boundary_region, "VectorPotentialImag", a_imag_model);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_real_inner(conductor_inner, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_imag_inner(conductor_inner, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_inner(conductor_inner, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_inner(conductor_inner, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling);

    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationInner>
        solve_a_real_inner(conductor_inner, omega, 1.0,
                           "VectorPotentialReal", "VectorPotentialImag",
                           "SourceCurrentDensityReal", "ElectricPotentialGradientReal",
                           "CurlNuBReal", "VectorPotentialChangeRateReal",
                           sigma_relaxation_scaling, sigma_relaxation_floor,
                           magnetic_diagonal_scaling, reference_dt_pseudo,
                           relaxation_scaling, max_change_rate);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationInner>
        solve_a_imag_inner(conductor_inner, omega, -1.0,
                           "VectorPotentialImag", "VectorPotentialReal",
                           "SourceCurrentDensityImag", "ElectricPotentialGradientImag",
                           "CurlNuBImag", "VectorPotentialChangeRateImag",
                           sigma_relaxation_scaling, sigma_relaxation_floor,
                           magnetic_diagonal_scaling, reference_dt_pseudo,
                           relaxation_scaling, max_change_rate);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyCoupledEquationInner>
        solve_a_coupled_inner(conductor_inner, omega,
                              "VectorPotentialReal", "VectorPotentialImag",
                              "SourceCurrentDensityReal", "SourceCurrentDensityImag",
                              "ElectricPotentialGradientReal", "ElectricPotentialGradientImag",
                              "CurlNuBReal", "CurlNuBImag",
                              "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
                              sigma_relaxation_scaling, sigma_relaxation_floor,
                              magnetic_diagonal_scaling, reference_dt_pseudo,
                              relaxation_scaling, max_change_rate);
    SimpleDynamics<electromagnetics::FrequencyAEquationResidualDiagnostic>
        diagnose_a_real(conductor_body, omega, 1.0,
                        "VectorPotentialReal", "VectorPotentialImag",
                        "SourceCurrentDensityReal", "ElectricPotentialGradientReal",
                        "CurlNuBReal",
                        "AEquationResidualVectorReal", "AEquationRelativeResidualReal");
    SimpleDynamics<electromagnetics::FrequencyAEquationResidualDiagnostic>
        diagnose_a_imag(conductor_body, omega, -1.0,
                        "VectorPotentialImag", "VectorPotentialReal",
                        "SourceCurrentDensityImag", "ElectricPotentialGradientImag",
                        "CurlNuBImag",
                        "AEquationResidualVectorImag", "AEquationRelativeResidualImag");

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(conductor_body, "VectorPotentialReal");
    write_states.addToWrite<Vecd>(conductor_body, "VectorPotentialImag");
    write_states.addToWrite<AngularVecd>(conductor_body, "VectorPotentialCurlReal");
    write_states.addToWrite<AngularVecd>(conductor_body, "VectorPotentialCurlImag");
    write_states.addToWrite<Vecd>(conductor_body, "CurlNuBReal");
    write_states.addToWrite<Vecd>(conductor_body, "CurlNuBImag");
    write_states.addToWrite<Vecd>(conductor_body, "SourceCurrentDensityReal");
    write_states.addToWrite<Vecd>(conductor_body, "SourceCurrentDensityImag");
    write_states.addToWrite<Vecd>(conductor_body, "AEquationResidualVectorReal");
    write_states.addToWrite<Vecd>(conductor_body, "AEquationResidualVectorImag");
    write_states.addToWrite<Real>(conductor_body, "AEquationRelativeResidualReal");
    write_states.addToWrite<Real>(conductor_body, "AEquationRelativeResidualImag");

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_frequency_em.exec();
    assign_manufactured_auxiliary_fields.exec();

    const std::string history_path =
        io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_history.csv";
    std::ofstream history_file(history_path, std::ios::out | std::ios::trunc);
    history_file << std::setprecision(12);
    history_file << "solver_mode,solve_phase,iteration,component,"
                 << "mean_curl_nu_b_norm,max_curl_nu_b_norm,"
                 << "mean_residual_norm,max_residual_norm,"
                 << "mean_relative_residual,max_relative_residual,"
                 << "mean_change_rate_norm,max_change_rate_norm\n";
    history_file.close();

    auto write_history_row = [&](const std::string &mode,
                                 const std::string &phase,
                                 int iteration,
                                 const std::string &component,
                                 const ComponentIterationMetrics &metrics)
    {
        std::ofstream append_file(history_path, std::ios::out | std::ios::app);
        append_file << std::setprecision(12);
        append_file << mode << ","
                    << phase << ","
                    << iteration << ","
                    << component << ","
                    << metrics.mean_curl_nu_b_norm << ","
                    << metrics.max_curl_nu_b_norm << ","
                    << metrics.mean_residual_norm << ","
                    << metrics.max_residual_norm << ","
                    << metrics.mean_relative_residual << ","
                    << metrics.max_relative_residual << ","
                    << metrics.mean_change_rate_norm << ","
                    << metrics.max_change_rate_norm << "\n";
    };
    auto refresh_equation_residual_diagnostics = [&]()
    {
        diagnose_a_real.exec();
        diagnose_a_imag.exec();
    };

    auto write_particle_header = [&](const std::string &path)
    {
        std::ofstream particle_file(path, std::ios::out | std::ios::trunc);
        particle_file << std::setprecision(12);
        particle_file << "component,particle_id,x,y,z,"
                      << "a_x,a_y,a_z,a_exact_x,a_exact_y,a_exact_z,a_error,"
                      << "curl_a_x,curl_a_y,curl_a_z,b_error,"
                      << "curl_nu_b_x,curl_nu_b_y,curl_nu_b_z,curl_nu_b_norm,"
                      << "residual_x,residual_y,residual_z,residual_norm,relative_residual,"
                      << "change_rate_x,change_rate_y,change_rate_z,change_rate_norm\n";
    };

    auto set_summary_meta = [&](ASolverSummary &summary,
                                const std::string &component,
                                const std::string &mode,
                                const std::string &phase,
                                Real initial_scale)
    {
        summary.component_name = component;
        summary.solver_mode = mode;
        summary.solve_phase = phase;
        summary.iterations = solve_iterations;
        summary.initial_guess_scale = initial_scale;
    };

    std::vector<ASolverSummary> summaries;

    if (run_single)
    {
        // Phase 1: solve A_real with exact A_imag and exact grad(phi_real) fixed.
        set_initial_guess_real.exec();
        if (use_boundary_constraint)
        {
            constrain_a_real_boundary.exec();
        }
        curl_a_real_inner.exec();
        curl_nu_b_real_inner.exec();
        refresh_equation_residual_diagnostics();
        write_history_row("single", "phase1_real", 0, "A_real",
                          evaluate_real_iteration_metrics(conductor_body.getBaseParticles()));
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            curl_a_real_inner.exec();
            curl_nu_b_real_inner.exec();
            solve_a_real_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_real_boundary.exec();
            }
            bool sample_this_iter =
                ((iter + 1) % effective_history_interval == 0) || (iter + 1 == solve_iterations);
            if (sample_this_iter)
            {
                curl_a_real_inner.exec();
                curl_nu_b_real_inner.exec();
                refresh_equation_residual_diagnostics();
                write_history_row("single", "phase1_real", iter + 1, "A_real",
                                  evaluate_real_iteration_metrics(conductor_body.getBaseParticles()));
            }
        }
        curl_a_real_inner.exec();
        curl_nu_b_real_inner.exec();
        refresh_equation_residual_diagnostics();
        ASolverSummary real_summary =
            evaluate_real_component_summary(conductor_body.getBaseParticles());
        set_summary_meta(real_summary, "A_real", "single", "phase1_real", initial_guess_scale_real);
        summaries.push_back(real_summary);
        if (write_particles)
        {
            const std::string particles_phase1_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_particles_phase1_real.csv";
            write_particle_header(particles_phase1_path);
            write_particle_diagnostics_real(
                particles_phase1_path, conductor_body.getBaseParticles(), "A_real");
        }

        // Reset exact auxiliaries before the A_imag phase.
        assign_manufactured_auxiliary_fields.exec();

        // Phase 2: solve A_imag with exact A_real and exact grad(phi_imag) fixed.
        set_initial_guess_imag.exec();
        if (use_boundary_constraint)
        {
            constrain_a_imag_boundary.exec();
        }
        curl_a_imag_inner.exec();
        curl_nu_b_imag_inner.exec();
        refresh_equation_residual_diagnostics();
        write_history_row("single", "phase2_imag", 0, "A_imag",
                          evaluate_imag_iteration_metrics(conductor_body.getBaseParticles()));
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            curl_a_imag_inner.exec();
            curl_nu_b_imag_inner.exec();
            solve_a_imag_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_imag_boundary.exec();
            }
            bool sample_this_iter =
                ((iter + 1) % effective_history_interval == 0) || (iter + 1 == solve_iterations);
            if (sample_this_iter)
            {
                curl_a_imag_inner.exec();
                curl_nu_b_imag_inner.exec();
                refresh_equation_residual_diagnostics();
                write_history_row("single", "phase2_imag", iter + 1, "A_imag",
                                  evaluate_imag_iteration_metrics(conductor_body.getBaseParticles()));
            }
        }
        curl_a_imag_inner.exec();
        curl_nu_b_imag_inner.exec();
        refresh_equation_residual_diagnostics();
        ASolverSummary imag_summary =
            evaluate_imag_component_summary(conductor_body.getBaseParticles());
        set_summary_meta(imag_summary, "A_imag", "single", "phase2_imag", initial_guess_scale_imag);
        summaries.push_back(imag_summary);
        if (write_particles)
        {
            const std::string particles_phase2_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_particles_phase2_imag.csv";
            write_particle_header(particles_phase2_path);
            write_particle_diagnostics_imag(
                particles_phase2_path, conductor_body.getBaseParticles(), "A_imag");
        }
    }

    if (run_coupled)
    {
        assign_manufactured_auxiliary_fields.exec();
        set_initial_guess_real.exec();
        set_initial_guess_imag.exec();
        if (use_boundary_constraint)
        {
            constrain_a_real_boundary.exec();
            constrain_a_imag_boundary.exec();
        }
        curl_a_real_inner.exec();
        curl_nu_b_real_inner.exec();
        curl_a_imag_inner.exec();
        curl_nu_b_imag_inner.exec();
        refresh_equation_residual_diagnostics();
        write_history_row("coupled", "coupled", 0, "A_real",
                          evaluate_real_iteration_metrics(conductor_body.getBaseParticles()));
        write_history_row("coupled", "coupled", 0, "A_imag",
                          evaluate_imag_iteration_metrics(conductor_body.getBaseParticles()));
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            curl_a_real_inner.exec();
            curl_nu_b_real_inner.exec();
            curl_a_imag_inner.exec();
            curl_nu_b_imag_inner.exec();
            solve_a_coupled_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_real_boundary.exec();
                constrain_a_imag_boundary.exec();
            }
            bool sample_this_iter =
                ((iter + 1) % effective_history_interval == 0) || (iter + 1 == solve_iterations);
            if (sample_this_iter)
            {
                curl_a_real_inner.exec();
                curl_nu_b_real_inner.exec();
                curl_a_imag_inner.exec();
                curl_nu_b_imag_inner.exec();
                refresh_equation_residual_diagnostics();
                write_history_row("coupled", "coupled", iter + 1, "A_real",
                                  evaluate_real_iteration_metrics(conductor_body.getBaseParticles()));
                write_history_row("coupled", "coupled", iter + 1, "A_imag",
                                  evaluate_imag_iteration_metrics(conductor_body.getBaseParticles()));
            }
        }
        curl_a_real_inner.exec();
        curl_nu_b_real_inner.exec();
        curl_a_imag_inner.exec();
        curl_nu_b_imag_inner.exec();
        refresh_equation_residual_diagnostics();
        ASolverSummary real_summary =
            evaluate_real_component_summary(conductor_body.getBaseParticles());
        set_summary_meta(real_summary, "A_real", "coupled", "coupled", initial_guess_scale_real);
        summaries.push_back(real_summary);

        ASolverSummary imag_summary =
            evaluate_imag_component_summary(conductor_body.getBaseParticles());
        set_summary_meta(imag_summary, "A_imag", "coupled", "coupled", initial_guess_scale_imag);
        summaries.push_back(imag_summary);

        if (write_particles)
        {
            const std::string particles_coupled_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_particles_coupled.csv";
            write_particle_header(particles_coupled_path);
            write_particle_diagnostics_real(
                particles_coupled_path, conductor_body.getBaseParticles(), "A_real");
            write_particle_diagnostics_imag(
                particles_coupled_path, conductor_body.getBaseParticles(), "A_imag");
        }
    }

    if (write_vtp)
    {
        write_states.writeToFile(0.0);
    }

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "component,solver_mode,solve_phase,total_particles,iterations,initial_guess_scale,"
                 << "mean_a_error,max_a_error,"
                 << "mean_b_error,max_b_error,"
                 << "mean_curl_nu_b_norm,max_curl_nu_b_norm,"
                 << "mean_residual_norm,max_residual_norm,"
                 << "mean_relative_residual,max_relative_residual,"
                 << "mean_change_rate_norm,max_change_rate_norm\n";

    auto write_summary = [&](const ASolverSummary &summary)
    {
        summary_file << summary.component_name << ","
                     << summary.solver_mode << ","
                     << summary.solve_phase << ","
                     << summary.total_particles << ","
                     << summary.iterations << ","
                     << summary.initial_guess_scale << ","
                     << summary.mean_a_error << ","
                     << summary.max_a_error << ","
                     << summary.mean_b_error << ","
                     << summary.max_b_error << ","
                     << summary.mean_curl_nu_b_norm << ","
                     << summary.max_curl_nu_b_norm << ","
                     << summary.mean_residual_norm << ","
                     << summary.max_residual_norm << ","
                     << summary.mean_relative_residual << ","
                     << summary.max_relative_residual << ","
                     << summary.mean_change_rate_norm << ","
                     << summary.max_change_rate_norm << "\n";
    };
    for (const ASolverSummary &summary : summaries)
    {
        write_summary(summary);
    }
    summary_file.flush();

    for (const ASolverSummary &summary : summaries)
    {
        std::cout << std::scientific << std::setprecision(6)
                  << "[em-freq-a-solver-verify-summary] mode=" << summary.solver_mode
                  << ", phase=" << summary.solve_phase
                  << ", component=" << summary.component_name
                  << ", mean|A-A_exact|=" << summary.mean_a_error
                  << ", mean|curlA-B|=" << summary.mean_b_error
                  << ", mean|curlNuB|=" << summary.mean_curl_nu_b_norm
                  << ", mean|residual|=" << summary.mean_residual_norm
                  << ", mean(relative residual)=" << summary.mean_relative_residual
                  << ", mean|dA/dtau|=" << summary.mean_change_rate_norm << std::endl;
    }
    std::cout << "[em-freq-a-solver-verify] summary file: "
              << io_environment.OutputFolder() + "/em_frequency_a_solver_manufactured_summary.csv"
              << std::endl;
    std::cout << "[em-freq-a-solver-verify] history file: "
              << history_path << std::endl;

    return 0;
}
