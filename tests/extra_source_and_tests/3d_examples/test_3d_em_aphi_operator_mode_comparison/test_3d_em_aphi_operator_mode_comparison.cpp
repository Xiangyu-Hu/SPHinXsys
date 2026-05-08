/**
 * @file test_3d_em_aphi_operator_mode_comparison.cpp
 * @brief Compare discrete Laplace and weak-inspired curl-curl operator energies
 * on smooth and checkerboard-like vector-potential modes.
 */

#include "sphinxsys.h"
#include "electromagnetic_aphi_operator_comparison_eigen.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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

const Real dp_0 = get_env_real_local("EM_APHI_OPERATOR_COMPARE_DP", 1.0);
const Real body_length = get_env_real_local("EM_APHI_OPERATOR_COMPARE_LENGTH", 10.0);
const Real body_height = get_env_real_local("EM_APHI_OPERATOR_COMPARE_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_APHI_OPERATOR_COMPARE_WIDTH", 6.0);
const Real boundary_width = get_env_real_local("EM_APHI_OPERATOR_COMPARE_BOUNDARY_WIDTH", 3.0 * dp_0);
const Real magnetic_reluctivity = get_env_real_local("EM_APHI_OPERATOR_COMPARE_NU", 1.0);
const Real gauge_penalty = get_env_real_local("EM_APHI_OPERATOR_COMPARE_GAUGE_PENALTY", 0.0);
const bool enable_gauge_penalty =
    get_env_bool_local("EM_APHI_OPERATOR_COMPARE_ENABLE_GAUGE_PENALTY", false);

const Vec3d body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_lower_bound = body_center - body_halfsize;
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

using Complex = electromagnetics::Complex;

Complex smooth_ax(const Vec3d &position)
{
    Real xi = Pi * (position[0] - body_lower_bound[0]) / body_length;
    Real yi = Pi * (position[1] - body_lower_bound[1]) / body_height;
    Real zi = Pi * (position[2] - body_lower_bound[2]) / body_width;
    return Complex(-2.0 * (Pi / body_height) *
                       std::sin(xi) * std::sin(xi) *
                       std::sin(yi) * std::cos(yi) *
                       std::sin(zi),
                   0.0);
}

Complex smooth_ay(const Vec3d &position)
{
    Real xi = Pi * (position[0] - body_lower_bound[0]) / body_length;
    Real yi = Pi * (position[1] - body_lower_bound[1]) / body_height;
    Real zi = Pi * (position[2] - body_lower_bound[2]) / body_width;
    return Complex(2.0 * (Pi / body_length) *
                       std::sin(xi) * std::cos(xi) *
                       std::sin(yi) * std::sin(yi) *
                       std::sin(zi),
                   0.0);
}

int parity_from_position(const Vec3d &position)
{
    int ix = static_cast<int>(std::llround((position[0] - 0.5 * dp_0) / dp_0));
    int iy = static_cast<int>(std::llround((position[1] - 0.5 * dp_0) / dp_0));
    int iz = static_cast<int>(std::llround((position[2] - 0.5 * dp_0) / dp_0));
    return ((ix + iy + iz) % 2 == 0) ? 1 : -1;
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_APHI_OPERATOR_COMPARE_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-aphi-operator-compare-config] dp=" << dp_0
              << ", nu=" << magnetic_reluctivity
              << ", gauge_penalty=" << gauge_penalty
              << ", enable_gauge_penalty=" << enable_gauge_penalty
              << std::endl;

    SolidBody conductor_body(sph_system, makeShared<ConductorShape>("Conductor"));
    conductor_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    conductor_body.defineMaterial<Solid>();
    conductor_body.defineBodyLevelSetShape();
    conductor_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(conductor_body);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = conductor_body.getBaseParticles();
    size_t total_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    StdVec<Real> nu(total_particles, magnetic_reluctivity);
    StdVec<Complex> smooth_ax_field(total_particles);
    StdVec<Complex> smooth_ay_field(total_particles);
    StdVec<Complex> smooth_az_field(total_particles, Complex(0.0, 0.0));
    StdVec<Complex> checker_ax_field(total_particles);
    StdVec<Complex> checker_ay_field(total_particles, Complex(0.0, 0.0));
    StdVec<Complex> checker_az_field(total_particles, Complex(0.0, 0.0));

    for (size_t i = 0; i != total_particles; ++i)
    {
        smooth_ax_field[i] = smooth_ax(positions[i]);
        smooth_ay_field[i] = smooth_ay(positions[i]);
        checker_ax_field[i] = Complex(static_cast<Real>(parity_from_position(positions[i])), 0.0);
    }

    electromagnetics::DiscreteOperatorComparisonParameters parameters;
    parameters.gauge_penalty = gauge_penalty;
    parameters.enable_gauge_penalty = enable_gauge_penalty;

    electromagnetics::DiscreteEMOperatorComparator comparator(inner_relation, parameters);
    comparator.buildOperators(nu);

    electromagnetics::DiscreteModeEnergySummary smooth_summary =
        comparator.evaluateMode(smooth_ax_field, smooth_ay_field, smooth_az_field);
    electromagnetics::DiscreteModeEnergySummary checker_summary =
        comparator.evaluateMode(checker_ax_field, checker_ay_field, checker_az_field);

    auto normalized_energy = [](const electromagnetics::DiscreteModeEnergySummary &summary,
                                Real energy) -> Real
    {
        return energy / (summary.max_abs_a * summary.max_abs_a + TinyReal);
    };

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_aphi_operator_mode_comparison.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "mode,laplace_energy,weak_curlcurl_energy,laplace_energy_norm,weak_curlcurl_energy_norm,"
                 << "divergence_l2,curl_l2,max_abs_a,weak_to_laplace_ratio\n";

    auto write_summary = [&](const std::string &mode_name,
                             const electromagnetics::DiscreteModeEnergySummary &summary)
    {
        Real laplace_energy_norm = normalized_energy(summary, summary.laplace_energy);
        Real weak_energy_norm = normalized_energy(summary, summary.weak_curlcurl_energy);
        Real ratio = weak_energy_norm / (laplace_energy_norm + TinyReal);
        summary_file << mode_name << ","
                     << summary.laplace_energy << ","
                     << summary.weak_curlcurl_energy << ","
                     << laplace_energy_norm << ","
                     << weak_energy_norm << ","
                     << summary.divergence_l2 << ","
                     << summary.curl_l2 << ","
                     << summary.max_abs_a << ","
                     << ratio << "\n";
    };

    write_summary("smooth_divfree", smooth_summary);
    write_summary("checkerboard_ax", checker_summary);

    std::cout << "[em-aphi-operator-compare] smooth ratio="
              << normalized_energy(smooth_summary, smooth_summary.weak_curlcurl_energy) /
                     (normalized_energy(smooth_summary, smooth_summary.laplace_energy) + TinyReal)
              << ", checkerboard ratio="
              << normalized_energy(checker_summary, checker_summary.weak_curlcurl_energy) /
                     (normalized_energy(checker_summary, checker_summary.laplace_energy) + TinyReal)
              << std::endl;

    return 0;
}
