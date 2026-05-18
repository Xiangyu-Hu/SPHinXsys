/**
 * @file test_3d_em_frequency_phi_only_linear.cpp
 * @brief Verify frequency-domain scalar-potential operators with fixed A.
 *
 * This case keeps A_real/A_imag fixed (no A-equation solve) and verifies only
 * phi-path operators for both real/imag parts:
 *   1) grad(phi) against analytical gradient (linear/quadratic selectable)
 *   2) sigma*grad(phi) against analytical vector
 *   3) dphi/dtau from scalar relaxation (with consistent analytical source term)
 */

#include "sphinxsys.h"
#include "aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.hpp"
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

const Real dp_0 = get_env_real_local("EM_FREQ_PHI_VERIFY_DP", 1.0);
const Real total_length = get_env_real_local("EM_FREQ_PHI_VERIFY_TOTAL_LENGTH", 16.0);
const Real body_height = get_env_real_local("EM_FREQ_PHI_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_FREQ_PHI_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_FREQ_PHI_VERIFY_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real interface_shell_thickness =
    get_env_real_local("EM_FREQ_PHI_VERIFY_INTERFACE_SHELL_THICKNESS", 2.5 * dp_0);
const Real conductivity = get_env_real_local("EM_FREQ_PHI_VERIFY_SIGMA", 1.0);
const Real rho_cp = get_env_real_local("EM_FREQ_PHI_VERIFY_RHO_CP", 1.0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);

const Vec3d fixed_a_real(get_env_real_local("EM_FREQ_PHI_VERIFY_AREAL_X", 0.0),
                         get_env_real_local("EM_FREQ_PHI_VERIFY_AREAL_Y", 1.0e-4),
                         get_env_real_local("EM_FREQ_PHI_VERIFY_AREAL_Z", 0.0));
const Vec3d fixed_a_imag(get_env_real_local("EM_FREQ_PHI_VERIFY_AIMAG_X", 0.0),
                         get_env_real_local("EM_FREQ_PHI_VERIFY_AIMAG_Y", -2.0e-4),
                         get_env_real_local("EM_FREQ_PHI_VERIFY_AIMAG_Z", 0.0));

const Vec3d target_grad_phi_real(get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_REAL_X", 0.4),
                                 get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_REAL_Y", -0.2),
                                 get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_REAL_Z", 0.1));
const Vec3d target_grad_phi_imag(get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_IMAG_X", -0.3),
                                 get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_IMAG_Y", 0.1),
                                 get_env_real_local("EM_FREQ_PHI_VERIFY_GRADPHI_IMAG_Z", 0.2));
const bool use_quadratic_phi = get_env_bool_local("EM_FREQ_PHI_VERIFY_USE_QUADRATIC", false);
const Vec3d quadratic_coeff_phi_real(
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_REAL_X", 0.0),
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_REAL_Y", 0.0),
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_REAL_Z", 0.0));
const Vec3d quadratic_coeff_phi_imag(
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_IMAG_X", 0.0),
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_IMAG_Y", 0.0),
    get_env_real_local("EM_FREQ_PHI_VERIFY_K_IMAG_Z", 0.0));
const Real gradient_scaling =
    get_env_real_local("EM_FREQ_PHI_VERIFY_GRADIENT_SCALING", 1.0);
const Real laplacian_scaling =
    get_env_real_local("EM_FREQ_PHI_VERIFY_LAPLACIAN_SCALING", 1.0);
const bool use_inner_operator = get_env_bool_local("EM_FREQ_PHI_VERIFY_USE_INNER", true);
const bool use_contact_operator = get_env_bool_local("EM_FREQ_PHI_VERIFY_USE_CONTACT", true);

const Real half_total_length = 0.5 * total_length;
const Real half_body_length = 0.25 * total_length;
const Vec3d body_halfsize(half_body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d left_body_center(-half_body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d right_body_center(half_body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d gauge_origin(0.0, 0.5 * body_height, 0.5 * body_width);
BoundingBoxd system_domain_bounds(
    Vec3d(-half_total_length - boundary_width, -boundary_width, -boundary_width),
    Vec3d(half_total_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

class LeftBodyShape : public ComplexShape
{
  public:
    explicit LeftBodyShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(left_body_center), body_halfsize);
    }
};

class RightBodyShape : public ComplexShape
{
  public:
    explicit RightBodyShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(right_body_center), body_halfsize);
    }
};

class MergedBodyShape : public ComplexShape
{
  public:
    explicit MergedBodyShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(gauge_origin),
                               Vec3d(half_total_length, 0.5 * body_height, 0.5 * body_width));
    }
};

struct AnalyticalPhiField
{
    Vec3d linear_grad;
    Vec3d quadratic_coeff;

    Real evaluatePotential(const Vec3d &relative_position) const
    {
        Real x = relative_position[0];
        Real y = relative_position[1];
        Real z = relative_position[2];
        return linear_grad.dot(relative_position) +
               0.5 * (quadratic_coeff[0] * x * x +
                      quadratic_coeff[1] * y * y +
                      quadratic_coeff[2] * z * z);
    }

    Vec3d evaluateGradient(const Vec3d &relative_position) const
    {
        return linear_grad + Vec3d(quadratic_coeff[0] * relative_position[0],
                                   quadratic_coeff[1] * relative_position[1],
                                   quadratic_coeff[2] * relative_position[2]);
    }

    Real evaluateSource(Real sigma) const
    {
        // Source is chosen to balance div(sigma * grad(phi)) in the continuous equation.
        return -sigma * (quadratic_coeff[0] + quadratic_coeff[1] + quadratic_coeff[2]);
    }
};

const AnalyticalPhiField phi_real_model{
    target_grad_phi_real,
    use_quadratic_phi ? quadratic_coeff_phi_real : ZeroData<Vec3d>::value};
const AnalyticalPhiField phi_imag_model{
    target_grad_phi_imag,
    use_quadratic_phi ? quadratic_coeff_phi_imag : ZeroData<Vec3d>::value};

class AssignFixedAAndAnalyticalPhi : public LocalDynamics
{
  public:
    explicit AssignFixedAAndAnalyticalPhi(SPHBody &sph_body,
                                          const Vec3d &fixed_a_real,
                                          const Vec3d &fixed_a_imag,
                                          const AnalyticalPhiField &phi_real_model,
                                          const AnalyticalPhiField &phi_imag_model,
                                          Real sigma,
                                          const Vec3d &origin)
        : LocalDynamics(sph_body),
          fixed_a_real_(fixed_a_real),
          fixed_a_imag_(fixed_a_imag),
          phi_real_model_(phi_real_model),
          phi_imag_model_(phi_imag_model),
          source_real_value_(phi_real_model.evaluateSource(sigma)),
          source_imag_value_(phi_imag_model.evaluateSource(sigma)),
          origin_(origin),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          a_real_(particles_->getVariableDataByName<Vecd>("VectorPotentialReal")),
          a_imag_(particles_->getVariableDataByName<Vecd>("VectorPotentialImag")),
          phi_real_(particles_->getVariableDataByName<Real>("ElectricPotentialReal")),
          phi_imag_(particles_->getVariableDataByName<Real>("ElectricPotentialImag")),
          phi_source_real_(particles_->getVariableDataByName<Real>("ElectricPotentialSourceReal")),
          phi_source_imag_(particles_->getVariableDataByName<Real>("ElectricPotentialSourceImag")),
          grad_phi_real_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientReal")),
          grad_phi_imag_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientImag")),
          phi_rate_real_(particles_->getVariableDataByName<Real>("ElectricPotentialChangeRateReal")),
          phi_rate_imag_(particles_->getVariableDataByName<Real>("ElectricPotentialChangeRateImag"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d relative_position = positions_[index_i] - origin_;
        a_real_[index_i] = fixed_a_real_;
        a_imag_[index_i] = fixed_a_imag_;
        phi_real_[index_i] = phi_real_model_.evaluatePotential(relative_position);
        phi_imag_[index_i] = phi_imag_model_.evaluatePotential(relative_position);
        phi_source_real_[index_i] = source_real_value_;
        phi_source_imag_[index_i] = source_imag_value_;
        grad_phi_real_[index_i] = ZeroData<Vecd>::value;
        grad_phi_imag_[index_i] = ZeroData<Vecd>::value;
        phi_rate_real_[index_i] = 0.0;
        phi_rate_imag_[index_i] = 0.0;
    }

  protected:
    Vec3d fixed_a_real_;
    Vec3d fixed_a_imag_;
    AnalyticalPhiField phi_real_model_;
    AnalyticalPhiField phi_imag_model_;
    Real source_real_value_;
    Real source_imag_value_;
    Vec3d origin_;
    Vecd *positions_;
    Vecd *a_real_, *a_imag_;
    Real *phi_real_, *phi_imag_;
    Real *phi_source_real_, *phi_source_imag_;
    Vecd *grad_phi_real_, *grad_phi_imag_;
    Real *phi_rate_real_, *phi_rate_imag_;
};

class SetConstantConductivity : public LocalDynamics
{
  public:
    explicit SetConstantConductivity(SPHBody &sph_body, Real sigma)
        : LocalDynamics(sph_body),
          sigma_(sigma),
          electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        electrical_conductivity_[index_i] = sigma_;
    }

  protected:
    Real sigma_;
    Real *electrical_conductivity_;
};

struct PhiFrequencySummary
{
    std::string body_name;
    size_t total_particles = 0;
    size_t interface_particles = 0;
    size_t core_particles = 0;
    Real mean_grad_real_error = 0.0;
    Real mean_grad_imag_error = 0.0;
    Real max_grad_real_error = 0.0;
    Real max_grad_imag_error = 0.0;
    Real mean_sigma_grad_real_error = 0.0;
    Real mean_sigma_grad_imag_error = 0.0;
    Real max_sigma_grad_real_error = 0.0;
    Real max_sigma_grad_imag_error = 0.0;
    Real mean_phi_rate_real_abs = 0.0;
    Real mean_phi_rate_imag_abs = 0.0;
    Real max_phi_rate_real_abs = 0.0;
    Real max_phi_rate_imag_abs = 0.0;
    Real mean_phi_rate_real_inner_abs = 0.0;
    Real mean_phi_rate_real_contact_abs = 0.0;
    Real mean_phi_rate_real_balance_abs = 0.0;
    Real mean_phi_rate_imag_inner_abs = 0.0;
    Real mean_phi_rate_imag_contact_abs = 0.0;
    Real mean_phi_rate_imag_balance_abs = 0.0;
};

bool is_interface_particle(const Vec3d &position)
{
    return fabs(position[0]) <= interface_shell_thickness;
}

PhiFrequencySummary evaluate_phi_frequency_summary(const std::string &body_name,
                                                   BaseParticles &particles,
                                                   const AnalyticalPhiField &phi_real_model,
                                                   const AnalyticalPhiField &phi_imag_model,
                                                   Real sigma,
                                                   const std::vector<Real> &phi_rate_real_inner,
                                                   const std::vector<Real> &phi_rate_real_contact,
                                                   const std::vector<Real> &phi_rate_imag_inner,
                                                   const std::vector<Real> &phi_rate_imag_contact)
{
    PhiFrequencySummary summary;
    summary.body_name = body_name;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *grad_phi_real = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
    Vecd *grad_phi_imag = particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
    Real *phi_rate_real = particles.getVariableDataByName<Real>("ElectricPotentialChangeRateReal");
    Real *phi_rate_imag = particles.getVariableDataByName<Real>("ElectricPotentialChangeRateImag");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_grad_real_error = 0.0;
    Real sum_grad_imag_error = 0.0;
    Real sum_sigma_grad_real_error = 0.0;
    Real sum_sigma_grad_imag_error = 0.0;
    Real sum_phi_rate_real_abs = 0.0;
    Real sum_phi_rate_imag_abs = 0.0;
    Real sum_phi_rate_real_inner_abs = 0.0;
    Real sum_phi_rate_real_contact_abs = 0.0;
    Real sum_phi_rate_real_balance_abs = 0.0;
    Real sum_phi_rate_imag_inner_abs = 0.0;
    Real sum_phi_rate_imag_contact_abs = 0.0;
    Real sum_phi_rate_imag_balance_abs = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d relative_position = positions[i] - gauge_origin;
        Vec3d expected_grad_real = phi_real_model.evaluateGradient(relative_position);
        Vec3d expected_grad_imag = phi_imag_model.evaluateGradient(relative_position);
        Vec3d grad_real_error_vec = grad_phi_real[i] - expected_grad_real;
        Vec3d grad_imag_error_vec = grad_phi_imag[i] - expected_grad_imag;
        Real grad_real_error = grad_real_error_vec.norm();
        Real grad_imag_error = grad_imag_error_vec.norm();
        Real sigma_grad_real_error = (sigma * grad_real_error_vec).norm();
        Real sigma_grad_imag_error = (sigma * grad_imag_error_vec).norm();
        Real phi_rate_real_abs = fabs(phi_rate_real[i]);
        Real phi_rate_imag_abs = fabs(phi_rate_imag[i]);
        Real phi_rate_real_inner_value =
            phi_rate_real_inner.size() == total_real_particles ? phi_rate_real_inner[i] : 0.0;
        Real phi_rate_real_contact_value =
            phi_rate_real_contact.size() == total_real_particles ? phi_rate_real_contact[i] : 0.0;
        Real phi_rate_imag_inner_value =
            phi_rate_imag_inner.size() == total_real_particles ? phi_rate_imag_inner[i] : 0.0;
        Real phi_rate_imag_contact_value =
            phi_rate_imag_contact.size() == total_real_particles ? phi_rate_imag_contact[i] : 0.0;
        Real phi_rate_real_balance_abs =
            fabs(phi_rate_real[i] - (phi_rate_real_inner_value + phi_rate_real_contact_value));
        Real phi_rate_imag_balance_abs =
            fabs(phi_rate_imag[i] - (phi_rate_imag_inner_value + phi_rate_imag_contact_value));
        (void)is_interface_particle(positions[i]);

        summary.total_particles++;
        sum_grad_real_error += grad_real_error;
        sum_grad_imag_error += grad_imag_error;
        sum_sigma_grad_real_error += sigma_grad_real_error;
        sum_sigma_grad_imag_error += sigma_grad_imag_error;
        sum_phi_rate_real_abs += phi_rate_real_abs;
        sum_phi_rate_imag_abs += phi_rate_imag_abs;
        sum_phi_rate_real_inner_abs += fabs(phi_rate_real_inner_value);
        sum_phi_rate_real_contact_abs += fabs(phi_rate_real_contact_value);
        sum_phi_rate_real_balance_abs += phi_rate_real_balance_abs;
        sum_phi_rate_imag_inner_abs += fabs(phi_rate_imag_inner_value);
        sum_phi_rate_imag_contact_abs += fabs(phi_rate_imag_contact_value);
        sum_phi_rate_imag_balance_abs += phi_rate_imag_balance_abs;
        summary.max_grad_real_error = SMAX(summary.max_grad_real_error, grad_real_error);
        summary.max_grad_imag_error = SMAX(summary.max_grad_imag_error, grad_imag_error);
        summary.max_sigma_grad_real_error = SMAX(summary.max_sigma_grad_real_error, sigma_grad_real_error);
        summary.max_sigma_grad_imag_error = SMAX(summary.max_sigma_grad_imag_error, sigma_grad_imag_error);
        summary.max_phi_rate_real_abs = SMAX(summary.max_phi_rate_real_abs, phi_rate_real_abs);
        summary.max_phi_rate_imag_abs = SMAX(summary.max_phi_rate_imag_abs, phi_rate_imag_abs);
    }

    summary.mean_grad_real_error =
        sum_grad_real_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_grad_imag_error =
        sum_grad_imag_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_sigma_grad_real_error =
        sum_sigma_grad_real_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_sigma_grad_imag_error =
        sum_sigma_grad_imag_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_real_abs =
        sum_phi_rate_real_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_imag_abs =
        sum_phi_rate_imag_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_real_inner_abs =
        sum_phi_rate_real_inner_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_real_contact_abs =
        sum_phi_rate_real_contact_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_real_balance_abs =
        sum_phi_rate_real_balance_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_imag_inner_abs =
        sum_phi_rate_imag_inner_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_imag_contact_abs =
        sum_phi_rate_imag_contact_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_imag_balance_abs =
        sum_phi_rate_imag_balance_abs / (static_cast<Real>(summary.total_particles) + TinyReal);

    // Keep shell counters for CSV compatibility with other verification cases.
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (is_interface_particle(positions[i]))
        {
            summary.interface_particles++;
        }
        else
        {
            summary.core_particles++;
        }
    }
    return summary;
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_FREQ_PHI_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-phi-verify-config] dp=" << dp_0
              << ", sigma=" << conductivity
              << ", use_quadratic_phi=" << use_quadratic_phi
              << ", gradient_scaling=" << gradient_scaling
              << ", laplacian_scaling=" << laplacian_scaling
              << ", use_inner=" << use_inner_operator
              << ", use_contact=" << use_contact_operator
              << ", grad_phi_real=(" << target_grad_phi_real[0] << ","
              << target_grad_phi_real[1] << "," << target_grad_phi_real[2] << ")"
              << ", grad_phi_imag=(" << target_grad_phi_imag[0] << ","
              << target_grad_phi_imag[1] << "," << target_grad_phi_imag[2] << ")"
              << ", source_real=" << phi_real_model.evaluateSource(conductivity)
              << ", source_imag=" << phi_imag_model.evaluateSource(conductivity)
              << std::fixed << std::setprecision(6)
              << std::endl;

    SolidBody left_body(sph_system, makeShared<LeftBodyShape>("LeftBody"));
    left_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    left_body.defineMaterial<Solid>();
    left_body.defineBodyLevelSetShape();
    left_body.generateParticles<BaseParticles, Lattice>();

    SolidBody right_body(sph_system, makeShared<RightBodyShape>("RightBody"));
    right_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    right_body.defineMaterial<Solid>();
    right_body.defineBodyLevelSetShape();
    right_body.generateParticles<BaseParticles, Lattice>();

    SolidBody merged_body(sph_system, makeShared<MergedBodyShape>("MergedBody"));
    merged_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    merged_body.defineMaterial<Solid>();
    merged_body.defineBodyLevelSetShape();
    merged_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation left_inner(left_body);
    InnerRelation right_inner(right_body);
    InnerRelation merged_inner(merged_body);
    ContactRelation left_contact(left_body, {&right_body});
    ContactRelation right_contact(right_body, {&left_body});

    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_left_em(left_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_right_em(right_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_merged_em(merged_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<SetConstantConductivity> set_left_sigma(left_body, conductivity);
    SimpleDynamics<SetConstantConductivity> set_right_sigma(right_body, conductivity);
    SimpleDynamics<SetConstantConductivity> set_merged_sigma(merged_body, conductivity);
    SimpleDynamics<AssignFixedAAndAnalyticalPhi> assign_left_fields(
        left_body, fixed_a_real, fixed_a_imag, phi_real_model, phi_imag_model, conductivity, gauge_origin);
    SimpleDynamics<AssignFixedAAndAnalyticalPhi> assign_right_fields(
        right_body, fixed_a_real, fixed_a_imag, phi_real_model, phi_imag_model, conductivity, gauge_origin);
    SimpleDynamics<AssignFixedAAndAnalyticalPhi> assign_merged_fields(
        merged_body, fixed_a_real, fixed_a_imag, phi_real_model, phi_imag_model, conductivity, gauge_origin);

    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_real_left_inner(left_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_real_right_inner(right_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_real_merged_inner(merged_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        grad_phi_real_left_contact(left_contact, "ElectricPotentialReal", "ElectricPotentialGradientReal", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        grad_phi_real_right_contact(right_contact, "ElectricPotentialReal", "ElectricPotentialGradientReal", gradient_scaling);

    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_imag_left_inner(left_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_imag_right_inner(right_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        grad_phi_imag_merged_inner(merged_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        grad_phi_imag_left_contact(left_contact, "ElectricPotentialImag", "ElectricPotentialGradientImag", gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        grad_phi_imag_right_contact(right_contact, "ElectricPotentialImag", "ElectricPotentialGradientImag", gradient_scaling);

    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_real_relax_left_inner(left_inner, "ElectricPotentialReal", "ElectricPotentialSourceReal", "ElectricPotentialChangeRateReal", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_real_relax_right_inner(right_inner, "ElectricPotentialReal", "ElectricPotentialSourceReal", "ElectricPotentialChangeRateReal", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_real_relax_merged_inner(merged_inner, "ElectricPotentialReal", "ElectricPotentialSourceReal", "ElectricPotentialChangeRateReal", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationContactByName>
        phi_real_relax_left_contact(left_contact, "ElectricPotentialReal", "ElectricPotentialChangeRateReal", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationContactByName>
        phi_real_relax_right_contact(right_contact, "ElectricPotentialReal", "ElectricPotentialChangeRateReal", laplacian_scaling);

    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_imag_relax_left_inner(left_inner, "ElectricPotentialImag", "ElectricPotentialSourceImag", "ElectricPotentialChangeRateImag", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_imag_relax_right_inner(right_inner, "ElectricPotentialImag", "ElectricPotentialSourceImag", "ElectricPotentialChangeRateImag", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationInnerByName>
        phi_imag_relax_merged_inner(merged_inner, "ElectricPotentialImag", "ElectricPotentialSourceImag", "ElectricPotentialChangeRateImag", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationContactByName>
        phi_imag_relax_left_contact(left_contact, "ElectricPotentialImag", "ElectricPotentialChangeRateImag", laplacian_scaling);
    InteractionDynamics<electromagnetics::ScalarRelaxationContactByName>
        phi_imag_relax_right_contact(right_contact, "ElectricPotentialImag", "ElectricPotentialChangeRateImag", laplacian_scaling);

    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_left_em.exec();
    initialize_right_em.exec();
    initialize_merged_em.exec();
    set_left_sigma.exec();
    set_right_sigma.exec();
    set_merged_sigma.exec();
    assign_left_fields.exec();
    assign_right_fields.exec();
    assign_merged_fields.exec();

    auto reset_phi_rates = [](BaseParticles &particles)
    {
        Real *phi_rate_real = particles.getVariableDataByName<Real>("ElectricPotentialChangeRateReal");
        Real *phi_rate_imag = particles.getVariableDataByName<Real>("ElectricPotentialChangeRateImag");
        size_t total_real_particles = particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            phi_rate_real[i] = 0.0;
            phi_rate_imag[i] = 0.0;
        }
    };
    auto snapshot_phi_rates = [](BaseParticles &particles, const std::string &field_name)
    {
        Real *field = particles.getVariableDataByName<Real>(field_name);
        size_t total_real_particles = particles.TotalRealParticles();
        std::vector<Real> snapshot(total_real_particles, 0.0);
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            snapshot[i] = field[i];
        }
        return snapshot;
    };

    if (use_inner_operator)
    {
        grad_phi_real_left_inner.exec();
        grad_phi_real_right_inner.exec();
        grad_phi_real_merged_inner.exec();
        grad_phi_imag_left_inner.exec();
        grad_phi_imag_right_inner.exec();
        grad_phi_imag_merged_inner.exec();
    }
    if (use_contact_operator)
    {
        grad_phi_real_left_contact.exec();
        grad_phi_real_right_contact.exec();
        grad_phi_imag_left_contact.exec();
        grad_phi_imag_right_contact.exec();
    }

    reset_phi_rates(left_body.getBaseParticles());
    reset_phi_rates(right_body.getBaseParticles());
    reset_phi_rates(merged_body.getBaseParticles());
    if (use_inner_operator)
    {
        phi_real_relax_left_inner.exec();
        phi_real_relax_right_inner.exec();
        phi_real_relax_merged_inner.exec();
        phi_imag_relax_left_inner.exec();
        phi_imag_relax_right_inner.exec();
        phi_imag_relax_merged_inner.exec();
    }
    std::vector<Real> left_phi_rate_real_inner =
        snapshot_phi_rates(left_body.getBaseParticles(), "ElectricPotentialChangeRateReal");
    std::vector<Real> left_phi_rate_imag_inner =
        snapshot_phi_rates(left_body.getBaseParticles(), "ElectricPotentialChangeRateImag");
    std::vector<Real> right_phi_rate_real_inner =
        snapshot_phi_rates(right_body.getBaseParticles(), "ElectricPotentialChangeRateReal");
    std::vector<Real> right_phi_rate_imag_inner =
        snapshot_phi_rates(right_body.getBaseParticles(), "ElectricPotentialChangeRateImag");
    std::vector<Real> merged_phi_rate_real_inner =
        snapshot_phi_rates(merged_body.getBaseParticles(), "ElectricPotentialChangeRateReal");
    std::vector<Real> merged_phi_rate_imag_inner =
        snapshot_phi_rates(merged_body.getBaseParticles(), "ElectricPotentialChangeRateImag");

    reset_phi_rates(left_body.getBaseParticles());
    reset_phi_rates(right_body.getBaseParticles());
    reset_phi_rates(merged_body.getBaseParticles());
    if (use_contact_operator)
    {
        phi_real_relax_left_contact.exec();
        phi_real_relax_right_contact.exec();
        phi_imag_relax_left_contact.exec();
        phi_imag_relax_right_contact.exec();
    }
    std::vector<Real> left_phi_rate_real_contact =
        snapshot_phi_rates(left_body.getBaseParticles(), "ElectricPotentialChangeRateReal");
    std::vector<Real> left_phi_rate_imag_contact =
        snapshot_phi_rates(left_body.getBaseParticles(), "ElectricPotentialChangeRateImag");
    std::vector<Real> right_phi_rate_real_contact =
        snapshot_phi_rates(right_body.getBaseParticles(), "ElectricPotentialChangeRateReal");
    std::vector<Real> right_phi_rate_imag_contact =
        snapshot_phi_rates(right_body.getBaseParticles(), "ElectricPotentialChangeRateImag");
    std::vector<Real> merged_phi_rate_real_contact(
        merged_body.getBaseParticles().TotalRealParticles(), 0.0);
    std::vector<Real> merged_phi_rate_imag_contact(
        merged_body.getBaseParticles().TotalRealParticles(), 0.0);

    reset_phi_rates(left_body.getBaseParticles());
    reset_phi_rates(right_body.getBaseParticles());
    reset_phi_rates(merged_body.getBaseParticles());
    if (use_inner_operator)
    {
        phi_real_relax_left_inner.exec();
        phi_real_relax_right_inner.exec();
        phi_real_relax_merged_inner.exec();
        phi_imag_relax_left_inner.exec();
        phi_imag_relax_right_inner.exec();
        phi_imag_relax_merged_inner.exec();
    }
    if (use_contact_operator)
    {
        phi_real_relax_left_contact.exec();
        phi_real_relax_right_contact.exec();
        phi_imag_relax_left_contact.exec();
        phi_imag_relax_right_contact.exec();
    }

    write_states.writeToFile(0.0);

    PhiFrequencySummary left_summary = evaluate_phi_frequency_summary(
        "LeftBody", left_body.getBaseParticles(), phi_real_model, phi_imag_model, conductivity,
        left_phi_rate_real_inner, left_phi_rate_real_contact,
        left_phi_rate_imag_inner, left_phi_rate_imag_contact);
    PhiFrequencySummary right_summary = evaluate_phi_frequency_summary(
        "RightBody", right_body.getBaseParticles(), phi_real_model, phi_imag_model, conductivity,
        right_phi_rate_real_inner, right_phi_rate_real_contact,
        right_phi_rate_imag_inner, right_phi_rate_imag_contact);
    PhiFrequencySummary merged_summary = evaluate_phi_frequency_summary(
        "MergedBody", merged_body.getBaseParticles(), phi_real_model, phi_imag_model, conductivity,
        merged_phi_rate_real_inner, merged_phi_rate_real_contact,
        merged_phi_rate_imag_inner, merged_phi_rate_imag_contact);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_frequency_phi_only_linear_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file
        << "body,total_particles,interface_particles,core_particles,"
        << "mean_grad_real_error,max_grad_real_error,mean_grad_imag_error,max_grad_imag_error,"
        << "mean_sigma_grad_real_error,max_sigma_grad_real_error,"
        << "mean_sigma_grad_imag_error,max_sigma_grad_imag_error,"
        << "mean_phi_rate_real_abs,max_phi_rate_real_abs,"
        << "mean_phi_rate_imag_abs,max_phi_rate_imag_abs,"
        << "mean_phi_rate_real_inner_abs,mean_phi_rate_real_contact_abs,mean_phi_rate_real_balance_abs,"
        << "mean_phi_rate_imag_inner_abs,mean_phi_rate_imag_contact_abs,mean_phi_rate_imag_balance_abs\n";

    auto write_summary = [&](const PhiFrequencySummary &summary)
    {
        summary_file << summary.body_name << ","
                     << summary.total_particles << ","
                     << summary.interface_particles << ","
                     << summary.core_particles << ","
                     << summary.mean_grad_real_error << ","
                     << summary.max_grad_real_error << ","
                     << summary.mean_grad_imag_error << ","
                     << summary.max_grad_imag_error << ","
                     << summary.mean_sigma_grad_real_error << ","
                     << summary.max_sigma_grad_real_error << ","
                     << summary.mean_sigma_grad_imag_error << ","
                     << summary.max_sigma_grad_imag_error << ","
                     << summary.mean_phi_rate_real_abs << ","
                     << summary.max_phi_rate_real_abs << ","
                     << summary.mean_phi_rate_imag_abs << ","
                     << summary.max_phi_rate_imag_abs << ","
                     << summary.mean_phi_rate_real_inner_abs << ","
                     << summary.mean_phi_rate_real_contact_abs << ","
                     << summary.mean_phi_rate_real_balance_abs << ","
                     << summary.mean_phi_rate_imag_inner_abs << ","
                     << summary.mean_phi_rate_imag_contact_abs << ","
                     << summary.mean_phi_rate_imag_balance_abs << "\n";
    };
    write_summary(left_summary);
    write_summary(right_summary);
    write_summary(merged_summary);
    summary_file.flush();

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-phi-verify-summary] LeftBody mean|grad(phi_re)-ref|="
              << left_summary.mean_grad_real_error
              << ", mean|sigma*grad(phi_re)-ref|=" << left_summary.mean_sigma_grad_real_error
              << ", mean|dphi_re/dtau|=" << left_summary.mean_phi_rate_real_abs
              << " (inner=" << left_summary.mean_phi_rate_real_inner_abs
              << ", contact=" << left_summary.mean_phi_rate_real_contact_abs
              << ", balance=" << left_summary.mean_phi_rate_real_balance_abs << ")"
              << ", mean|grad(phi_im)-ref|=" << left_summary.mean_grad_imag_error
              << ", mean|dphi_im/dtau|=" << left_summary.mean_phi_rate_imag_abs
              << " (inner=" << left_summary.mean_phi_rate_imag_inner_abs
              << ", contact=" << left_summary.mean_phi_rate_imag_contact_abs
              << ", balance=" << left_summary.mean_phi_rate_imag_balance_abs << ")"
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-phi-verify-summary] RightBody mean|grad(phi_re)-ref|="
              << right_summary.mean_grad_real_error
              << ", mean|sigma*grad(phi_re)-ref|=" << right_summary.mean_sigma_grad_real_error
              << ", mean|dphi_re/dtau|=" << right_summary.mean_phi_rate_real_abs
              << " (inner=" << right_summary.mean_phi_rate_real_inner_abs
              << ", contact=" << right_summary.mean_phi_rate_real_contact_abs
              << ", balance=" << right_summary.mean_phi_rate_real_balance_abs << ")"
              << ", mean|grad(phi_im)-ref|=" << right_summary.mean_grad_imag_error
              << ", mean|dphi_im/dtau|=" << right_summary.mean_phi_rate_imag_abs
              << " (inner=" << right_summary.mean_phi_rate_imag_inner_abs
              << ", contact=" << right_summary.mean_phi_rate_imag_contact_abs
              << ", balance=" << right_summary.mean_phi_rate_imag_balance_abs << ")"
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-phi-verify-summary] MergedBody mean|grad(phi_re)-ref|="
              << merged_summary.mean_grad_real_error
              << ", mean|sigma*grad(phi_re)-ref|=" << merged_summary.mean_sigma_grad_real_error
              << ", mean|dphi_re/dtau|=" << merged_summary.mean_phi_rate_real_abs
              << " (inner=" << merged_summary.mean_phi_rate_real_inner_abs
              << ", contact=" << merged_summary.mean_phi_rate_real_contact_abs
              << ", balance=" << merged_summary.mean_phi_rate_real_balance_abs << ")"
              << ", mean|grad(phi_im)-ref|=" << merged_summary.mean_grad_imag_error
              << ", mean|dphi_im/dtau|=" << merged_summary.mean_phi_rate_imag_abs
              << " (inner=" << merged_summary.mean_phi_rate_imag_inner_abs
              << ", contact=" << merged_summary.mean_phi_rate_imag_contact_abs
              << ", balance=" << merged_summary.mean_phi_rate_imag_balance_abs << ")"
              << std::endl;
    std::cout << "[em-freq-phi-verify] summary file: "
              << io_environment.OutputFolder() + "/em_frequency_phi_only_linear_summary.csv"
              << std::endl;

    return 0;
}
