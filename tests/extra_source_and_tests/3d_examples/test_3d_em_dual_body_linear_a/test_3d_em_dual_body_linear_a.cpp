/**
 * @file test_3d_em_dual_body_linear_a.cpp
 * @brief Minimal operator verification for a split-domain electromagnetic body.
 *
 * The purpose of this case is to verify the discrete interface operators before
 * returning to TEAM7. Two adjacent bodies share the same material properties and
 * are initialized with an analytical linear vector potential field A = 0.5 * B x r.
 *
 * For constant material reluctivity and constant target magnetic flux density B:
 *   1) curl(A) = B        (constant)
 *   2) curl(nu * curl(A)) = 0
 *
 * This case isolates the body-body interface contribution and reports errors in
 * both the interface shell and the body core.
 */

#include "sphinxsys.h"
#include "electromagnetic_team7_aphi_dynamics.hpp"
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

const Real dp_0 = get_env_real_local("EM_VERIFY_DP", 1.0);
const Real total_length = get_env_real_local("EM_VERIFY_TOTAL_LENGTH", 16.0);
const Real body_height = get_env_real_local("EM_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_VERIFY_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real interface_shell_thickness =
    get_env_real_local("EM_VERIFY_INTERFACE_SHELL_THICKNESS", 2.5 * dp_0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Vec3d target_b(get_env_real_local("EM_VERIFY_BX", 0.0),
                     get_env_real_local("EM_VERIFY_BY", 0.0),
                     get_env_real_local("EM_VERIFY_BZ", 1.0));
const bool use_inner_operator = get_env_bool_local("EM_VERIFY_USE_INNER", true);
const bool use_contact_operator = get_env_bool_local("EM_VERIFY_USE_CONTACT", true);

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

class AssignLinearVectorPotential : public LocalDynamics
{
  public:
    explicit AssignLinearVectorPotential(SPHBody &sph_body,
                                         const Vec3d &target_b,
                                         const Vec3d &gauge_origin)
        : LocalDynamics(sph_body),
          target_b_(target_b),
          gauge_origin_(gauge_origin),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          vector_potential_(particles_->getVariableDataByName<Vecd>("VectorPotential")),
          vector_potential_prev_(particles_->getVariableDataByName<Vecd>("VectorPotentialPrevious"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d relative_position = positions_[index_i] - gauge_origin_;
        Vec3d analytical_a = 0.5 * target_b_.cross(relative_position);
        vector_potential_[index_i] = analytical_a;
        vector_potential_prev_[index_i] = analytical_a;
    }

  protected:
    Vec3d target_b_;
    Vec3d gauge_origin_;
    Vecd *positions_;
    Vecd *vector_potential_;
    Vecd *vector_potential_prev_;
};

struct BodyOperatorSummary
{
    std::string body_name;
    size_t total_particles = 0;
    size_t interface_particles = 0;
    size_t core_particles = 0;
    Real mean_b_error = 0.0;
    Real max_b_error = 0.0;
    Real mean_b_error_interface = 0.0;
    Real max_b_error_interface = 0.0;
    Real mean_b_error_core = 0.0;
    Real max_b_error_core = 0.0;
    Real mean_curl_nu_b_norm = 0.0;
    Real max_curl_nu_b_norm = 0.0;
    Real mean_curl_nu_b_norm_interface = 0.0;
    Real max_curl_nu_b_norm_interface = 0.0;
    Real mean_curl_nu_b_norm_core = 0.0;
    Real max_curl_nu_b_norm_core = 0.0;
};

Vec3d angular_to_vec(const AngularVecd &value)
{
    return Vec3d(value[0], value[1], value[2]);
}

bool is_interface_particle(const Vec3d &position)
{
    return fabs(position[0]) <= interface_shell_thickness;
}

BodyOperatorSummary evaluate_body_summary(const std::string &body_name,
                                          BaseParticles &particles,
                                          const Vec3d &target_b)
{
    BodyOperatorSummary summary;
    summary.body_name = body_name;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *vector_potential_curl =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>("CurlNuB");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_b_error = 0.0;
    Real sum_b_error_interface = 0.0;
    Real sum_b_error_core = 0.0;
    Real sum_curl_nu_b_norm = 0.0;
    Real sum_curl_nu_b_norm_interface = 0.0;
    Real sum_curl_nu_b_norm_core = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d discrete_b = angular_to_vec(vector_potential_curl[i]);
        Real b_error = (discrete_b - target_b).norm();
        Real curl_nu_b_norm = curl_nu_b[i].norm();
        bool interface_particle = is_interface_particle(positions[i]);

        summary.total_particles++;
        sum_b_error += b_error;
        sum_curl_nu_b_norm += curl_nu_b_norm;
        summary.max_b_error = SMAX(summary.max_b_error, b_error);
        summary.max_curl_nu_b_norm = SMAX(summary.max_curl_nu_b_norm, curl_nu_b_norm);

        if (interface_particle)
        {
            summary.interface_particles++;
            sum_b_error_interface += b_error;
            sum_curl_nu_b_norm_interface += curl_nu_b_norm;
            summary.max_b_error_interface = SMAX(summary.max_b_error_interface, b_error);
            summary.max_curl_nu_b_norm_interface =
                SMAX(summary.max_curl_nu_b_norm_interface, curl_nu_b_norm);
        }
        else
        {
            summary.core_particles++;
            sum_b_error_core += b_error;
            sum_curl_nu_b_norm_core += curl_nu_b_norm;
            summary.max_b_error_core = SMAX(summary.max_b_error_core, b_error);
            summary.max_curl_nu_b_norm_core =
                SMAX(summary.max_curl_nu_b_norm_core, curl_nu_b_norm);
        }
    }

    summary.mean_b_error = sum_b_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_b_error_interface =
        sum_b_error_interface / (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_b_error_core =
        sum_b_error_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.mean_curl_nu_b_norm =
        sum_curl_nu_b_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_curl_nu_b_norm_interface =
        sum_curl_nu_b_norm_interface /
        (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_curl_nu_b_norm_core =
        sum_curl_nu_b_norm_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    return summary;
}

void write_particle_diagnostics(const std::string &file_path,
                                const std::string &body_name,
                                BaseParticles &particles,
                                const Vec3d &target_b)
{
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    file << std::setprecision(12);

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *vector_potential_curl =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>("CurlNuB");

    size_t total_real_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d discrete_b = angular_to_vec(vector_potential_curl[i]);
        Real b_error = (discrete_b - target_b).norm();
        Real curl_nu_b_norm = curl_nu_b[i].norm();
        file << body_name << ","
             << i << ","
             << positions[i][0] << ","
             << positions[i][1] << ","
             << positions[i][2] << ","
             << fabs(positions[i][0]) << ","
             << static_cast<int>(is_interface_particle(positions[i])) << ","
             << discrete_b[0] << ","
             << discrete_b[1] << ","
             << discrete_b[2] << ","
             << b_error << ","
             << curl_nu_b[i][0] << ","
             << curl_nu_b[i][1] << ","
             << curl_nu_b[i][2] << ","
             << curl_nu_b_norm << "\n";
    }
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-verify-config] dp=" << dp_0
              << ", total_length=" << total_length
              << ", height=" << body_height
              << ", width=" << body_width
              << ", interface_shell_thickness=" << interface_shell_thickness
              << ", use_inner=" << use_inner_operator
              << ", use_contact=" << use_contact_operator
              << ", output_tag=" << (output_tag.empty() ? std::string("default") : output_tag)
              << ", target_B=(" << target_b[0] << ","
              << target_b[1] << "," << target_b[2] << ")"
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

    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_left_em(left_body, 0.0, 1.0, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_right_em(right_body, 0.0, 1.0, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_merged_em(merged_body, 0.0, 1.0, magnetic_reluctivity);
    SimpleDynamics<AssignLinearVectorPotential>
        assign_left_linear_a(left_body, target_b, gauge_origin);
    SimpleDynamics<AssignLinearVectorPotential>
        assign_right_linear_a(right_body, target_b, gauge_origin);
    SimpleDynamics<AssignLinearVectorPotential>
        assign_merged_linear_a(merged_body, target_b, gauge_origin);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_left_inner(left_inner, "VectorPotential", "VectorPotentialCurl", 1.0);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_right_inner(right_inner, "VectorPotential", "VectorPotentialCurl", 1.0);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_left_contact(left_contact, "VectorPotential", "VectorPotentialCurl", 1.0);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_right_contact(right_contact, "VectorPotential", "VectorPotentialCurl", 1.0);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_merged_inner(merged_inner, "VectorPotential", "VectorPotentialCurl", 1.0);

    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_left_inner(left_inner, "VectorPotentialCurl", "CurlNuB", 1.0);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_right_inner(right_inner, "VectorPotentialCurl", "CurlNuB", 1.0);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_merged_inner(merged_inner, "VectorPotentialCurl", "CurlNuB", 1.0);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_left_contact(left_contact, "VectorPotentialCurl", "CurlNuB", 1.0);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_right_contact(right_contact, "VectorPotentialCurl", "CurlNuB", 1.0);

    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_left_em.exec();
    initialize_right_em.exec();
    initialize_merged_em.exec();
    assign_left_linear_a.exec();
    assign_right_linear_a.exec();
    assign_merged_linear_a.exec();

    if (use_inner_operator)
    {
        curl_a_left_inner.exec();
        curl_a_right_inner.exec();
        curl_a_merged_inner.exec();
        curl_nu_b_left_inner.exec();
        curl_nu_b_right_inner.exec();
        curl_nu_b_merged_inner.exec();
    }
    if (use_contact_operator)
    {
        curl_a_left_contact.exec();
        curl_a_right_contact.exec();
        curl_nu_b_left_contact.exec();
        curl_nu_b_right_contact.exec();
    }

    write_states.writeToFile(0.0);

    BodyOperatorSummary left_summary =
        evaluate_body_summary("LeftBody", left_body.getBaseParticles(), target_b);
    BodyOperatorSummary right_summary =
        evaluate_body_summary("RightBody", right_body.getBaseParticles(), target_b);
    BodyOperatorSummary merged_summary =
        evaluate_body_summary("MergedBody", merged_body.getBaseParticles(), target_b);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_dual_body_linear_a_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file
        << "body,total_particles,interface_particles,core_particles,"
        << "mean_b_error,max_b_error,mean_b_error_interface,max_b_error_interface,"
        << "mean_b_error_core,max_b_error_core,"
        << "mean_curl_nu_b_norm,max_curl_nu_b_norm,"
        << "mean_curl_nu_b_norm_interface,max_curl_nu_b_norm_interface,"
        << "mean_curl_nu_b_norm_core,max_curl_nu_b_norm_core\n";

    auto write_summary = [&](const BodyOperatorSummary &summary)
    {
        summary_file << summary.body_name << ","
                     << summary.total_particles << ","
                     << summary.interface_particles << ","
                     << summary.core_particles << ","
                     << summary.mean_b_error << ","
                     << summary.max_b_error << ","
                     << summary.mean_b_error_interface << ","
                     << summary.max_b_error_interface << ","
                     << summary.mean_b_error_core << ","
                     << summary.max_b_error_core << ","
                     << summary.mean_curl_nu_b_norm << ","
                     << summary.max_curl_nu_b_norm << ","
                     << summary.mean_curl_nu_b_norm_interface << ","
                     << summary.max_curl_nu_b_norm_interface << ","
                     << summary.mean_curl_nu_b_norm_core << ","
                     << summary.max_curl_nu_b_norm_core << "\n";
    };
    write_summary(left_summary);
    write_summary(right_summary);
    write_summary(merged_summary);
    summary_file.flush();

    std::ofstream particle_file(
        io_environment.OutputFolder() + "/em_dual_body_linear_a_particles.csv",
        std::ios::out | std::ios::trunc);
    particle_file << std::setprecision(12);
    particle_file
        << "body,particle_id,x,y,z,abs_interface_distance,is_interface,"
        << "curl_a_x,curl_a_y,curl_a_z,b_error,"
        << "curl_nu_b_x,curl_nu_b_y,curl_nu_b_z,curl_nu_b_norm\n";
    particle_file.close();

    write_particle_diagnostics(io_environment.OutputFolder() + "/em_dual_body_linear_a_particles.csv",
                               "LeftBody", left_body.getBaseParticles(), target_b);
    write_particle_diagnostics(io_environment.OutputFolder() + "/em_dual_body_linear_a_particles.csv",
                               "RightBody", right_body.getBaseParticles(), target_b);
    write_particle_diagnostics(io_environment.OutputFolder() + "/em_dual_body_linear_a_particles.csv",
                               "MergedBody", merged_body.getBaseParticles(), target_b);

    std::cout << std::scientific << std::setprecision(6)
              << "[em-verify-summary] LeftBody mean|curlA-B|=" << left_summary.mean_b_error
              << ", interface mean|curlA-B|=" << left_summary.mean_b_error_interface
              << ", mean|curlNuB|=" << left_summary.mean_curl_nu_b_norm
              << ", interface mean|curlNuB|="
              << left_summary.mean_curl_nu_b_norm_interface
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-verify-summary] RightBody mean|curlA-B|=" << right_summary.mean_b_error
              << ", interface mean|curlA-B|=" << right_summary.mean_b_error_interface
              << ", mean|curlNuB|=" << right_summary.mean_curl_nu_b_norm
              << ", interface mean|curlNuB|="
              << right_summary.mean_curl_nu_b_norm_interface
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-verify-summary] MergedBody mean|curlA-B|=" << merged_summary.mean_b_error
              << ", interface-plane mean|curlA-B|=" << merged_summary.mean_b_error_interface
              << ", mean|curlNuB|=" << merged_summary.mean_curl_nu_b_norm
              << ", interface-plane mean|curlNuB|="
              << merged_summary.mean_curl_nu_b_norm_interface
              << std::endl;
    std::cout << "[em-verify] summary file: "
              << io_environment.OutputFolder() + "/em_dual_body_linear_a_summary.csv"
              << std::endl;

    return 0;
}
