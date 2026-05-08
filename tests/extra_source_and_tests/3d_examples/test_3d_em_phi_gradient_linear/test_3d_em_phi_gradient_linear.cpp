/**
 * @file test_3d_em_phi_gradient_linear.cpp
 * @brief Verify scalar-potential operators with linear analytical phi.
 *
 * Analytical field:
 *   phi(x) = grad_phi_target dot (x - origin)
 *
 * Expected:
 *   1) grad(phi) = grad_phi_target (constant)
 *   2) div(sigma grad(phi)) = 0 for constant sigma
 *      -> ElectricPotentialChangeRate should be near zero when source term is zero.
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

const Real dp_0 = get_env_real_local("EM_PHI_VERIFY_DP", 1.0);
const Real total_length = get_env_real_local("EM_PHI_VERIFY_TOTAL_LENGTH", 16.0);
const Real body_height = get_env_real_local("EM_PHI_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_PHI_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_PHI_VERIFY_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real interface_shell_thickness =
    get_env_real_local("EM_PHI_VERIFY_INTERFACE_SHELL_THICKNESS", 2.5 * dp_0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Real conductivity = get_env_real_local("EM_PHI_VERIFY_SIGMA", 1.0);
const Vec3d target_grad_phi(get_env_real_local("EM_PHI_VERIFY_GRADPHI_X", 0.4),
                            get_env_real_local("EM_PHI_VERIFY_GRADPHI_Y", -0.2),
                            get_env_real_local("EM_PHI_VERIFY_GRADPHI_Z", 0.1));
const bool use_inner_operator = get_env_bool_local("EM_PHI_VERIFY_USE_INNER", true);
const bool use_contact_operator = get_env_bool_local("EM_PHI_VERIFY_USE_CONTACT", true);

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

class AssignLinearElectricPotential : public LocalDynamics
{
  public:
    explicit AssignLinearElectricPotential(SPHBody &sph_body,
                                           const Vec3d &target_grad_phi,
                                           const Vec3d &origin)
        : LocalDynamics(sph_body),
          target_grad_phi_(target_grad_phi),
          origin_(origin),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d relative_position = positions_[index_i] - origin_;
        electric_potential_[index_i] = target_grad_phi_.dot(relative_position);
    }

  protected:
    Vec3d target_grad_phi_;
    Vec3d origin_;
    Vecd *positions_;
    Real *electric_potential_;
};

class SetConstantConductivity : public LocalDynamics
{
  public:
    explicit SetConstantConductivity(SPHBody &sph_body, Real sigma)
        : LocalDynamics(sph_body),
          sigma_(sigma),
          electrical_conductivity_(
              particles_->getVariableDataByName<Real>("ElectricalConductivity"))
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

struct PhiOperatorSummary
{
    std::string body_name;
    size_t total_particles = 0;
    size_t interface_particles = 0;
    size_t core_particles = 0;
    Real mean_grad_error = 0.0;
    Real max_grad_error = 0.0;
    Real mean_grad_error_interface = 0.0;
    Real max_grad_error_interface = 0.0;
    Real mean_grad_error_core = 0.0;
    Real max_grad_error_core = 0.0;
    Real mean_phi_rate_abs = 0.0;
    Real max_phi_rate_abs = 0.0;
    Real mean_phi_rate_abs_interface = 0.0;
    Real max_phi_rate_abs_interface = 0.0;
    Real mean_phi_rate_abs_core = 0.0;
    Real max_phi_rate_abs_core = 0.0;
};

bool is_interface_particle(const Vec3d &position)
{
    return fabs(position[0]) <= interface_shell_thickness;
}

PhiOperatorSummary evaluate_phi_summary(const std::string &body_name,
                                        BaseParticles &particles,
                                        const Vec3d &target_grad_phi)
{
    PhiOperatorSummary summary;
    summary.body_name = body_name;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>("ElectricPotentialGradient");
    Real *phi_rate = particles.getVariableDataByName<Real>("ElectricPotentialChangeRate");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_grad_error = 0.0;
    Real sum_grad_error_interface = 0.0;
    Real sum_grad_error_core = 0.0;
    Real sum_phi_rate_abs = 0.0;
    Real sum_phi_rate_abs_interface = 0.0;
    Real sum_phi_rate_abs_core = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Real grad_error = (grad_phi[i] - target_grad_phi).norm();
        Real phi_rate_abs = fabs(phi_rate[i]);
        bool interface_particle = is_interface_particle(positions[i]);

        summary.total_particles++;
        sum_grad_error += grad_error;
        sum_phi_rate_abs += phi_rate_abs;
        summary.max_grad_error = SMAX(summary.max_grad_error, grad_error);
        summary.max_phi_rate_abs = SMAX(summary.max_phi_rate_abs, phi_rate_abs);

        if (interface_particle)
        {
            summary.interface_particles++;
            sum_grad_error_interface += grad_error;
            sum_phi_rate_abs_interface += phi_rate_abs;
            summary.max_grad_error_interface = SMAX(summary.max_grad_error_interface, grad_error);
            summary.max_phi_rate_abs_interface =
                SMAX(summary.max_phi_rate_abs_interface, phi_rate_abs);
        }
        else
        {
            summary.core_particles++;
            sum_grad_error_core += grad_error;
            sum_phi_rate_abs_core += phi_rate_abs;
            summary.max_grad_error_core = SMAX(summary.max_grad_error_core, grad_error);
            summary.max_phi_rate_abs_core =
                SMAX(summary.max_phi_rate_abs_core, phi_rate_abs);
        }
    }

    summary.mean_grad_error =
        sum_grad_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_grad_error_interface =
        sum_grad_error_interface / (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_grad_error_core =
        sum_grad_error_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.mean_phi_rate_abs =
        sum_phi_rate_abs / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_phi_rate_abs_interface =
        sum_phi_rate_abs_interface / (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_phi_rate_abs_core =
        sum_phi_rate_abs_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    return summary;
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_PHI_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-phi-verify-config] dp=" << dp_0
              << ", total_length=" << total_length
              << ", interface_shell_thickness=" << interface_shell_thickness
              << ", sigma=" << conductivity
              << ", use_inner=" << use_inner_operator
              << ", use_contact=" << use_contact_operator
              << ", target_grad_phi=(" << target_grad_phi[0] << ","
              << target_grad_phi[1] << "," << target_grad_phi[2] << ")"
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
        initialize_left_em(left_body, conductivity, 1.0, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_right_em(right_body, conductivity, 1.0, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_merged_em(merged_body, conductivity, 1.0, magnetic_reluctivity);
    SimpleDynamics<SetConstantConductivity>
        set_left_sigma(left_body, conductivity);
    SimpleDynamics<SetConstantConductivity>
        set_right_sigma(right_body, conductivity);
    SimpleDynamics<SetConstantConductivity>
        set_merged_sigma(merged_body, conductivity);
    SimpleDynamics<AssignLinearElectricPotential>
        assign_left_linear_phi(left_body, target_grad_phi, gauge_origin);
    SimpleDynamics<AssignLinearElectricPotential>
        assign_right_linear_phi(right_body, target_grad_phi, gauge_origin);
    SimpleDynamics<AssignLinearElectricPotential>
        assign_merged_linear_phi(merged_body, target_grad_phi, gauge_origin);

    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        grad_phi_left_inner(left_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        grad_phi_right_inner(right_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        grad_phi_merged_inner(merged_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientContact>
        grad_phi_left_contact(left_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientContact>
        grad_phi_right_contact(right_contact);

    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        phi_relax_left_inner(left_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        phi_relax_right_inner(right_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        phi_relax_merged_inner(merged_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationContact>
        phi_relax_left_contact(left_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationContact>
        phi_relax_right_contact(right_contact);

    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_left_em.exec();
    initialize_right_em.exec();
    initialize_merged_em.exec();
    set_left_sigma.exec();
    set_right_sigma.exec();
    set_merged_sigma.exec();
    assign_left_linear_phi.exec();
    assign_right_linear_phi.exec();
    assign_merged_linear_phi.exec();

    if (use_inner_operator)
    {
        grad_phi_left_inner.exec();
        grad_phi_right_inner.exec();
        grad_phi_merged_inner.exec();
        phi_relax_left_inner.exec();
        phi_relax_right_inner.exec();
        phi_relax_merged_inner.exec();
    }
    if (use_contact_operator)
    {
        grad_phi_left_contact.exec();
        grad_phi_right_contact.exec();
        phi_relax_left_contact.exec();
        phi_relax_right_contact.exec();
    }

    write_states.writeToFile(0.0);

    PhiOperatorSummary left_summary =
        evaluate_phi_summary("LeftBody", left_body.getBaseParticles(), target_grad_phi);
    PhiOperatorSummary right_summary =
        evaluate_phi_summary("RightBody", right_body.getBaseParticles(), target_grad_phi);
    PhiOperatorSummary merged_summary =
        evaluate_phi_summary("MergedBody", merged_body.getBaseParticles(), target_grad_phi);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_phi_gradient_linear_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "body,total_particles,interface_particles,core_particles,"
                 << "mean_grad_error,max_grad_error,mean_grad_error_interface,max_grad_error_interface,"
                 << "mean_grad_error_core,max_grad_error_core,"
                 << "mean_phi_rate_abs,max_phi_rate_abs,mean_phi_rate_abs_interface,max_phi_rate_abs_interface,"
                 << "mean_phi_rate_abs_core,max_phi_rate_abs_core\n";

    auto write_summary = [&](const PhiOperatorSummary &summary)
    {
        summary_file << summary.body_name << ","
                     << summary.total_particles << ","
                     << summary.interface_particles << ","
                     << summary.core_particles << ","
                     << summary.mean_grad_error << ","
                     << summary.max_grad_error << ","
                     << summary.mean_grad_error_interface << ","
                     << summary.max_grad_error_interface << ","
                     << summary.mean_grad_error_core << ","
                     << summary.max_grad_error_core << ","
                     << summary.mean_phi_rate_abs << ","
                     << summary.max_phi_rate_abs << ","
                     << summary.mean_phi_rate_abs_interface << ","
                     << summary.max_phi_rate_abs_interface << ","
                     << summary.mean_phi_rate_abs_core << ","
                     << summary.max_phi_rate_abs_core << "\n";
    };
    write_summary(left_summary);
    write_summary(right_summary);
    write_summary(merged_summary);
    summary_file.flush();

    std::cout << std::scientific << std::setprecision(6)
              << "[em-phi-verify-summary] LeftBody mean|gradphi-grad_ref|="
              << left_summary.mean_grad_error
              << ", mean|dphi_dtau|="
              << left_summary.mean_phi_rate_abs
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-phi-verify-summary] RightBody mean|gradphi-grad_ref|="
              << right_summary.mean_grad_error
              << ", mean|dphi_dtau|="
              << right_summary.mean_phi_rate_abs
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-phi-verify-summary] MergedBody mean|gradphi-grad_ref|="
              << merged_summary.mean_grad_error
              << ", mean|dphi_dtau|="
              << merged_summary.mean_phi_rate_abs
              << std::endl;
    std::cout << "[em-phi-verify] summary file: "
              << io_environment.OutputFolder() + "/em_phi_gradient_linear_summary.csv"
              << std::endl;

    return 0;
}
