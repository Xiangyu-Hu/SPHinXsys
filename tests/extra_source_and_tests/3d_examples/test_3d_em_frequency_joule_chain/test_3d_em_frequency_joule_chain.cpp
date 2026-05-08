/**
 * @file test_3d_em_frequency_joule_chain.cpp
 * @brief Verify frequency-domain A/phi -> E/J/Joule post-processing chain.
 *
 * This case does NOT solve A-phi equations. It prescribes constant complex A and
 * zero grad(phi), then checks:
 *   E_re =  omega * A_im
 *   E_im = -omega * A_re
 *   J_re = sigma * E_re
 *   J_im = sigma * E_im
 *   Q    = 0.5 * (J_re dot E_re + J_im dot E_im)
 */

#include "sphinxsys.h"
#include "electromagnetic_team7_aphi_frequency_dynamics.hpp"
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

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

std::string get_env_string_local(const std::string &name,
                                 const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

const Real dp_0 = get_env_real_local("EM_FREQ_VERIFY_DP", 1.0);
const Real body_length = get_env_real_local("EM_FREQ_VERIFY_LENGTH", 20.0);
const Real body_height = get_env_real_local("EM_FREQ_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_FREQ_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_FREQ_VERIFY_BOUNDARY_WIDTH", 3.0 * dp_0);
const Real sigma = get_env_real_local("EM_FREQ_VERIFY_SIGMA", 3.0e6);
const Real rho_cp = get_env_real_local("EM_FREQ_VERIFY_RHO_CP", 1.0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Real frequency_hz = get_env_real_local("EM_FREQ_VERIFY_FREQUENCY_HZ", 50.0);
const Real omega = 2.0 * Pi * frequency_hz;

const Real a_real_y = get_env_real_local("EM_FREQ_VERIFY_A_REAL_Y", 2.0e-5);
const Real a_imag_y = get_env_real_local("EM_FREQ_VERIFY_A_IMAG_Y", 5.0e-5);

const Vec3d body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
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

class AssignAnalyticalFrequencyFields : public LocalDynamics
{
  public:
    explicit AssignAnalyticalFrequencyFields(SPHBody &sph_body,
                                             const Vec3d &a_real_value,
                                             const Vec3d &a_imag_value)
        : LocalDynamics(sph_body),
          a_real_value_(a_real_value),
          a_imag_value_(a_imag_value),
          a_real_(particles_->getVariableDataByName<Vecd>("VectorPotentialReal")),
          a_imag_(particles_->getVariableDataByName<Vecd>("VectorPotentialImag")),
          grad_phi_real_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientReal")),
          grad_phi_imag_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradientImag"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        a_real_[index_i] = a_real_value_;
        a_imag_[index_i] = a_imag_value_;
        grad_phi_real_[index_i] = ZeroData<Vecd>::value;
        grad_phi_imag_[index_i] = ZeroData<Vecd>::value;
    }

  protected:
    Vec3d a_real_value_;
    Vec3d a_imag_value_;
    Vecd *a_real_;
    Vecd *a_imag_;
    Vecd *grad_phi_real_;
    Vecd *grad_phi_imag_;
};

struct FrequencyChainSummary
{
    size_t total_particles = 0;
    Real expected_j_real_y = 0.0;
    Real expected_j_imag_y = 0.0;
    Real expected_joule = 0.0;
    Real mean_abs_error_j_real_y = 0.0;
    Real max_abs_error_j_real_y = 0.0;
    Real mean_abs_error_j_imag_y = 0.0;
    Real max_abs_error_j_imag_y = 0.0;
    Real mean_abs_error_joule = 0.0;
    Real max_abs_error_joule = 0.0;
};

FrequencyChainSummary evaluate_frequency_chain(BaseParticles &particles,
                                               Real sigma,
                                               Real omega,
                                               Real a_real_y,
                                               Real a_imag_y)
{
    FrequencyChainSummary summary;
    summary.expected_j_real_y = sigma * omega * a_imag_y;
    summary.expected_j_imag_y = -sigma * omega * a_real_y;
    summary.expected_joule = 0.5 * sigma * omega * omega *
                             (a_real_y * a_real_y + a_imag_y * a_imag_y);

    Vecd *current_density_real = particles.getVariableDataByName<Vecd>("CurrentDensityReal");
    Vecd *current_density_imag = particles.getVariableDataByName<Vecd>("CurrentDensityImag");
    Real *joule_heat = particles.getVariableDataByName<Real>("JouleHeatSource");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_abs_error_j_real_y = 0.0;
    Real sum_abs_error_j_imag_y = 0.0;
    Real sum_abs_error_joule = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        summary.total_particles++;
        Real error_j_real_y = fabs(current_density_real[i][1] - summary.expected_j_real_y);
        Real error_j_imag_y = fabs(current_density_imag[i][1] - summary.expected_j_imag_y);
        Real error_joule = fabs(joule_heat[i] - summary.expected_joule);
        sum_abs_error_j_real_y += error_j_real_y;
        sum_abs_error_j_imag_y += error_j_imag_y;
        sum_abs_error_joule += error_joule;
        summary.max_abs_error_j_real_y = SMAX(summary.max_abs_error_j_real_y, error_j_real_y);
        summary.max_abs_error_j_imag_y = SMAX(summary.max_abs_error_j_imag_y, error_j_imag_y);
        summary.max_abs_error_joule = SMAX(summary.max_abs_error_joule, error_joule);
    }

    summary.mean_abs_error_j_real_y =
        sum_abs_error_j_real_y / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_abs_error_j_imag_y =
        sum_abs_error_j_imag_y / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_abs_error_joule =
        sum_abs_error_joule / (static_cast<Real>(summary.total_particles) + TinyReal);
    return summary;
}

} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_FREQ_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-verify-config] dp=" << dp_0
              << ", sigma=" << sigma
              << ", rho_cp=" << rho_cp
              << ", frequency_hz=" << frequency_hz
              << ", omega=" << omega
              << ", A_real_y=" << a_real_y
              << ", A_imag_y=" << a_imag_y
              << std::fixed << std::setprecision(6)
              << std::endl;

    SolidBody conductor_body(sph_system, makeShared<ConductorShape>("Conductor"));
    conductor_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    conductor_body.defineMaterial<Solid>();
    conductor_body.defineBodyLevelSetShape();
    conductor_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation conductor_inner(conductor_body);

    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_frequency_em(conductor_body, sigma, rho_cp, magnetic_reluctivity);
    SimpleDynamics<AssignAnalyticalFrequencyFields> assign_frequency_fields(
        conductor_body, Vec3d(0.0, a_real_y, 0.0), Vec3d(0.0, a_imag_y, 0.0));
    InteractionDynamics<electromagnetics::FrequencyElectricFieldCurrentAndJouleHeatInner>
        update_joule_from_frequency_fields(conductor_inner, omega);

    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_frequency_em.exec();
    assign_frequency_fields.exec();
    update_joule_from_frequency_fields.exec();

    write_states.writeToFile(0.0);

    FrequencyChainSummary summary = evaluate_frequency_chain(
        conductor_body.getBaseParticles(), sigma, omega, a_real_y, a_imag_y);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_frequency_joule_chain_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(16);
    summary_file << "total_particles,"
                 << "expected_j_real_y,expected_j_imag_y,expected_joule,"
                 << "mean_abs_error_j_real_y,max_abs_error_j_real_y,"
                 << "mean_abs_error_j_imag_y,max_abs_error_j_imag_y,"
                 << "mean_abs_error_joule,max_abs_error_joule\n";
    summary_file << summary.total_particles << ","
                 << summary.expected_j_real_y << ","
                 << summary.expected_j_imag_y << ","
                 << summary.expected_joule << ","
                 << summary.mean_abs_error_j_real_y << ","
                 << summary.max_abs_error_j_real_y << ","
                 << summary.mean_abs_error_j_imag_y << ","
                 << summary.max_abs_error_j_imag_y << ","
                 << summary.mean_abs_error_joule << ","
                 << summary.max_abs_error_joule << "\n";
    summary_file.flush();

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-verify-summary] expected Jy_real=" << summary.expected_j_real_y
              << ", mean|error(Jy_real)|=" << summary.mean_abs_error_j_real_y
              << ", max|error(Jy_real)|=" << summary.max_abs_error_j_real_y
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-verify-summary] expected Jy_imag=" << summary.expected_j_imag_y
              << ", mean|error(Jy_imag)|=" << summary.mean_abs_error_j_imag_y
              << ", max|error(Jy_imag)|=" << summary.max_abs_error_j_imag_y
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-verify-summary] expected Joule=" << summary.expected_joule
              << ", mean|error(Joule)|=" << summary.mean_abs_error_joule
              << ", max|error(Joule)|=" << summary.max_abs_error_joule
              << std::endl;
    std::cout << "[em-freq-verify] summary file: "
              << io_environment.OutputFolder() + "/em_frequency_joule_chain_summary.csv"
              << std::endl;

    return 0;
}
