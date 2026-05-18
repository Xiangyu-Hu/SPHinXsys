/**
 * @file test_3d_em_frequency_dual_body_linear_a.cpp
 * @brief Verify frequency-domain curl(A_real/A_imag) in split/merged bodies.
 *
 * Analytical field:
 *   A_real = 0.5 * (B_real x r), A_imag = 0.5 * (B_imag x r)
 *
 * Expected:
 *   curl(A_real) = B_real (constant)
 *   curl(A_imag) = B_imag (constant)
 *
 * This case focuses on first-order curl operator in frequency form and reports
 * split-vs-merged statistics with interface-shell decomposition.
 */

#include "sphinxsys.h"
#include "aphi_case_support/electromagnetic_team7_aphi_dynamics.hpp"
#include "aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.hpp"
#include "em_matched_stencil_body_relations.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
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

const Real dp_0 = get_env_real_local("EM_FREQ_DUAL_VERIFY_DP", 1.0);
const Real total_length = get_env_real_local("EM_FREQ_DUAL_VERIFY_TOTAL_LENGTH", 16.0);
const Real body_height = get_env_real_local("EM_FREQ_DUAL_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_FREQ_DUAL_VERIFY_WIDTH", 8.0);
const Real boundary_width =
    get_env_real_local("EM_FREQ_DUAL_VERIFY_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real interface_shell_thickness =
    get_env_real_local("EM_FREQ_DUAL_VERIFY_INTERFACE_SHELL_THICKNESS", 2.5 * dp_0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Real conductivity = get_env_real_local("EM_FREQ_DUAL_VERIFY_SIGMA", 1.0);
const Real rho_cp = get_env_real_local("EM_FREQ_DUAL_VERIFY_RHO_CP", 1.0);

const Vec3d target_b_real(get_env_real_local("EM_FREQ_DUAL_VERIFY_BREAL_X", 0.0),
                          get_env_real_local("EM_FREQ_DUAL_VERIFY_BREAL_Y", 0.0),
                          get_env_real_local("EM_FREQ_DUAL_VERIFY_BREAL_Z", 1.0));
const Vec3d target_b_imag(get_env_real_local("EM_FREQ_DUAL_VERIFY_BIMAG_X", 0.0),
                          get_env_real_local("EM_FREQ_DUAL_VERIFY_BIMAG_Y", 0.0),
                          get_env_real_local("EM_FREQ_DUAL_VERIFY_BIMAG_Z", 0.5));
const bool use_inner_operator = get_env_bool_local("EM_FREQ_DUAL_VERIFY_USE_INNER", true);
const bool use_contact_operator = get_env_bool_local("EM_FREQ_DUAL_VERIFY_USE_CONTACT", true);
const Real curl_scaling = get_env_real_local("EM_FREQ_DUAL_VERIFY_CURL_SCALING", 1.0);
const Real curl_nu_b_scaling = get_env_real_local("EM_FREQ_DUAL_VERIFY_CURL_NUB_SCALING", 1.0);
const bool curl_nu_b_symmetric_volume_weight =
    get_env_bool_local("EM_FREQ_DUAL_VERIFY_CURL_NUB_SYMMETRIC_VOL", false);
const std::string contact_correction_mode =
    get_env_string_local("EM_FREQ_DUAL_VERIFY_CONTACT_CORRECTION_MODE", "baseline");
const Real contact_damping_factor =
    get_env_real_local("EM_FREQ_DUAL_VERIFY_CONTACT_DAMPING_FACTOR", 0.5);
const Real contact_limit_factor =
    get_env_real_local("EM_FREQ_DUAL_VERIFY_CONTACT_LIMIT_FACTOR", 0.5);
/** 0=default. Bit1: contact stencil = inner (narrow). Bit2: inner stencil = contact (wide). Bit3=both. */
const int neighbor_stencil_mode =
    static_cast<int>(get_env_real_local("EM_FREQ_DUAL_VERIFY_NEIGHBOR_STENCIL", 0.0));
const bool position_consistency_diag = get_env_bool_local("EM_FREQ_DUAL_VERIFY_POSITION_DIAG", false);
const bool split_merged_field_diag = get_env_bool_local("EM_FREQ_DUAL_VERIFY_SPLIT_MERGED_FIELD_DIAG", false);
const bool minimal_io_mode = get_env_bool_local("EM_FREQ_DUAL_VERIFY_MINIMAL_IO", false);

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

class AssignLinearFrequencyVectorPotential : public LocalDynamics
{
  public:
    explicit AssignLinearFrequencyVectorPotential(SPHBody &sph_body,
                                                  const Vec3d &target_b_real,
                                                  const Vec3d &target_b_imag,
                                                  const Vec3d &gauge_origin)
        : LocalDynamics(sph_body),
          target_b_real_(target_b_real),
          target_b_imag_(target_b_imag),
          gauge_origin_(gauge_origin),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          vector_potential_real_(particles_->getVariableDataByName<Vecd>("VectorPotentialReal")),
          vector_potential_imag_(particles_->getVariableDataByName<Vecd>("VectorPotentialImag"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d relative_position = positions_[index_i] - gauge_origin_;
        vector_potential_real_[index_i] = 0.5 * target_b_real_.cross(relative_position);
        vector_potential_imag_[index_i] = 0.5 * target_b_imag_.cross(relative_position);
    }

  protected:
    Vec3d target_b_real_;
    Vec3d target_b_imag_;
    Vec3d gauge_origin_;
    Vecd *positions_;
    Vecd *vector_potential_real_;
    Vecd *vector_potential_imag_;
};

struct BodyCurlSummary
{
    std::string body_name;
    size_t total_particles = 0;
    size_t interface_particles = 0;
    size_t core_particles = 0;
    Real mean_b_real_error = 0.0;
    Real max_b_real_error = 0.0;
    Real mean_b_real_error_interface = 0.0;
    Real max_b_real_error_interface = 0.0;
    Real mean_b_real_error_core = 0.0;
    Real max_b_real_error_core = 0.0;
    Real mean_b_imag_error = 0.0;
    Real max_b_imag_error = 0.0;
    Real mean_b_imag_error_interface = 0.0;
    Real max_b_imag_error_interface = 0.0;
    Real mean_b_imag_error_core = 0.0;
    Real max_b_imag_error_core = 0.0;
    Real mean_curl_nu_b_real_norm = 0.0;
    Real max_curl_nu_b_real_norm = 0.0;
    Real mean_curl_nu_b_real_norm_interface = 0.0;
    Real max_curl_nu_b_real_norm_interface = 0.0;
    Real mean_curl_nu_b_real_norm_core = 0.0;
    Real max_curl_nu_b_real_norm_core = 0.0;
    Real mean_curl_nu_b_imag_norm = 0.0;
    Real max_curl_nu_b_imag_norm = 0.0;
    Real mean_curl_nu_b_imag_norm_interface = 0.0;
    Real max_curl_nu_b_imag_norm_interface = 0.0;
    Real mean_curl_nu_b_imag_norm_core = 0.0;
    Real max_curl_nu_b_imag_norm_core = 0.0;
};

Vec3d angular_to_vec(const AngularVecd &value)
{
    return Vec3d(value[0], value[1], value[2]);
}

bool is_interface_particle(const Vec3d &position)
{
    return fabs(position[0]) <= interface_shell_thickness;
}

BodyCurlSummary evaluate_body_curl_summary(const std::string &body_name,
                                           BaseParticles &particles,
                                           const Vec3d &target_b_real,
                                           const Vec3d &target_b_imag)
{
    BodyCurlSummary summary;
    summary.body_name = body_name;

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *curl_real =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
    AngularVecd *curl_imag =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
    Vecd *curl_nu_b_real = particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *curl_nu_b_imag = particles.getVariableDataByName<Vecd>("CurlNuBImag");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_b_real_error = 0.0;
    Real sum_b_real_error_interface = 0.0;
    Real sum_b_real_error_core = 0.0;
    Real sum_b_imag_error = 0.0;
    Real sum_b_imag_error_interface = 0.0;
    Real sum_b_imag_error_core = 0.0;
    Real sum_curl_nu_b_real_norm = 0.0;
    Real sum_curl_nu_b_real_norm_interface = 0.0;
    Real sum_curl_nu_b_real_norm_core = 0.0;
    Real sum_curl_nu_b_imag_norm = 0.0;
    Real sum_curl_nu_b_imag_norm_interface = 0.0;
    Real sum_curl_nu_b_imag_norm_core = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d discrete_b_real = angular_to_vec(curl_real[i]);
        Vec3d discrete_b_imag = angular_to_vec(curl_imag[i]);
        Real b_real_error = (discrete_b_real - target_b_real).norm();
        Real b_imag_error = (discrete_b_imag - target_b_imag).norm();
        Real curl_nu_b_real_norm = curl_nu_b_real[i].norm();
        Real curl_nu_b_imag_norm = curl_nu_b_imag[i].norm();
        bool interface_particle = is_interface_particle(positions[i]);

        summary.total_particles++;
        sum_b_real_error += b_real_error;
        sum_b_imag_error += b_imag_error;
        sum_curl_nu_b_real_norm += curl_nu_b_real_norm;
        sum_curl_nu_b_imag_norm += curl_nu_b_imag_norm;
        summary.max_b_real_error = SMAX(summary.max_b_real_error, b_real_error);
        summary.max_b_imag_error = SMAX(summary.max_b_imag_error, b_imag_error);
        summary.max_curl_nu_b_real_norm =
            SMAX(summary.max_curl_nu_b_real_norm, curl_nu_b_real_norm);
        summary.max_curl_nu_b_imag_norm =
            SMAX(summary.max_curl_nu_b_imag_norm, curl_nu_b_imag_norm);

        if (interface_particle)
        {
            summary.interface_particles++;
            sum_b_real_error_interface += b_real_error;
            sum_b_imag_error_interface += b_imag_error;
            sum_curl_nu_b_real_norm_interface += curl_nu_b_real_norm;
            sum_curl_nu_b_imag_norm_interface += curl_nu_b_imag_norm;
            summary.max_b_real_error_interface =
                SMAX(summary.max_b_real_error_interface, b_real_error);
            summary.max_b_imag_error_interface =
                SMAX(summary.max_b_imag_error_interface, b_imag_error);
            summary.max_curl_nu_b_real_norm_interface =
                SMAX(summary.max_curl_nu_b_real_norm_interface, curl_nu_b_real_norm);
            summary.max_curl_nu_b_imag_norm_interface =
                SMAX(summary.max_curl_nu_b_imag_norm_interface, curl_nu_b_imag_norm);
        }
        else
        {
            summary.core_particles++;
            sum_b_real_error_core += b_real_error;
            sum_b_imag_error_core += b_imag_error;
            sum_curl_nu_b_real_norm_core += curl_nu_b_real_norm;
            sum_curl_nu_b_imag_norm_core += curl_nu_b_imag_norm;
            summary.max_b_real_error_core = SMAX(summary.max_b_real_error_core, b_real_error);
            summary.max_b_imag_error_core = SMAX(summary.max_b_imag_error_core, b_imag_error);
            summary.max_curl_nu_b_real_norm_core =
                SMAX(summary.max_curl_nu_b_real_norm_core, curl_nu_b_real_norm);
            summary.max_curl_nu_b_imag_norm_core =
                SMAX(summary.max_curl_nu_b_imag_norm_core, curl_nu_b_imag_norm);
        }
    }

    summary.mean_b_real_error =
        sum_b_real_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_b_real_error_interface =
        sum_b_real_error_interface / (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_b_real_error_core =
        sum_b_real_error_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.mean_b_imag_error =
        sum_b_imag_error / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_b_imag_error_interface =
        sum_b_imag_error_interface / (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_b_imag_error_core =
        sum_b_imag_error_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.mean_curl_nu_b_real_norm =
        sum_curl_nu_b_real_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_curl_nu_b_real_norm_interface =
        sum_curl_nu_b_real_norm_interface /
        (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_curl_nu_b_real_norm_core =
        sum_curl_nu_b_real_norm_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.mean_curl_nu_b_imag_norm =
        sum_curl_nu_b_imag_norm / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.mean_curl_nu_b_imag_norm_interface =
        sum_curl_nu_b_imag_norm_interface /
        (static_cast<Real>(summary.interface_particles) + TinyReal);
    summary.mean_curl_nu_b_imag_norm_core =
        sum_curl_nu_b_imag_norm_core / (static_cast<Real>(summary.core_particles) + TinyReal);
    return summary;
}

std::vector<Vecd> snapshot_vec_field(BaseParticles &particles, const std::string &field_name)
{
    Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    size_t total_real_particles = particles.TotalRealParticles();
    std::vector<Vecd> snapshot(total_real_particles, ZeroData<Vecd>::value);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        snapshot[i] = field[i];
    }
    return snapshot;
}

/** Brute-force nearest merged particle distance for each split particle (N<=~2k: cheap). */
size_t nearest_merged_particle_index(Vecd *pos_merged, size_t nm, const Vecd &p)
{
    Real best = std::numeric_limits<Real>::max();
    size_t best_j = 0;
    for (size_t j = 0; j < nm; ++j)
    {
        Vec3d q(pos_merged[j][0], pos_merged[j][1], pos_merged[j][2]);
        Vec3d pi(p[0], p[1], p[2]);
        const Real d = (pi - q).norm();
        if (d < best)
        {
            best = d;
            best_j = j;
        }
    }
    return best_j;
}

void report_split_vs_merged_discrete_fields(const std::string &stage_label,
                                            BaseParticles &left_particles,
                                            BaseParticles &right_particles,
                                            BaseParticles &merged_particles,
                                            bool include_curl_nu_b)
{
    Vecd *pos_m = merged_particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *curl_r_m =
        merged_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
    AngularVecd *curl_i_m =
        merged_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
    Vecd *curl_nu_r_m = merged_particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *curl_nu_i_m = merged_particles.getVariableDataByName<Vecd>("CurlNuBImag");
    const size_t nm = merged_particles.TotalRealParticles();

    auto scan_body = [&](BaseParticles &split_p, const char *body_tag)
    {
        Vecd *pos = split_p.getVariableDataByName<Vecd>("Position");
        AngularVecd *curl_r = split_p.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
        AngularVecd *curl_i = split_p.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
        Vecd *curl_nu_r = split_p.getVariableDataByName<Vecd>("CurlNuBReal");
        Vecd *curl_nu_i = split_p.getVariableDataByName<Vecd>("CurlNuBImag");
        const size_t n = split_p.TotalRealParticles();
        Real sum_dr = 0, max_dr = 0, sum_dr_if = 0, max_dr_if = 0;
        Real sum_di = 0, max_di = 0, sum_di_if = 0, max_di_if = 0;
        Real sum_nr = 0, max_nr = 0, sum_nr_if = 0, max_nr_if = 0;
        Real sum_ni = 0, max_ni = 0, sum_ni_if = 0, max_ni_if = 0;
        size_t cnt_if = 0;
        for (size_t i = 0; i < n; ++i)
        {
            const size_t j = nearest_merged_particle_index(pos_m, nm, pos[i]);
            Vec3d pr_m = angular_to_vec(curl_r_m[j]);
            Vec3d pi_m = angular_to_vec(curl_i_m[j]);
            Vec3d pr = angular_to_vec(curl_r[i]);
            Vec3d pi = angular_to_vec(curl_i[i]);
            const Real dr = (pr - pr_m).norm();
            const Real di = (pi - pi_m).norm();
            sum_dr += dr;
            sum_di += di;
            max_dr = SMAX(max_dr, dr);
            max_di = SMAX(max_di, di);
            Real nr = 0, ni = 0;
            if (include_curl_nu_b)
            {
                nr = (curl_nu_r[i] - curl_nu_r_m[j]).norm();
                ni = (curl_nu_i[i] - curl_nu_i_m[j]).norm();
                sum_nr += nr;
                sum_ni += ni;
                max_nr = SMAX(max_nr, nr);
                max_ni = SMAX(max_ni, ni);
            }
            const Vec3d pos_i(pos[i][0], pos[i][1], pos[i][2]);
            if (is_interface_particle(pos_i))
            {
                ++cnt_if;
                sum_dr_if += dr;
                sum_di_if += di;
                max_dr_if = SMAX(max_dr_if, dr);
                max_di_if = SMAX(max_di_if, di);
                if (include_curl_nu_b)
                {
                    sum_nr_if += nr;
                    sum_ni_if += ni;
                    max_nr_if = SMAX(max_nr_if, nr);
                    max_ni_if = SMAX(max_ni_if, ni);
                }
            }
        }
        const Real mean_dr = sum_dr / (static_cast<Real>(n) + TinyReal);
        const Real mean_di = sum_di / (static_cast<Real>(n) + TinyReal);
        const Real mean_nr = sum_nr / (static_cast<Real>(n) + TinyReal);
        const Real mean_ni = sum_ni / (static_cast<Real>(n) + TinyReal);
        const Real mean_dr_if = sum_dr_if / (static_cast<Real>(cnt_if) + TinyReal);
        const Real mean_di_if = sum_di_if / (static_cast<Real>(cnt_if) + TinyReal);
        const Real mean_nr_if = sum_nr_if / (static_cast<Real>(cnt_if) + TinyReal);
        const Real mean_ni_if = sum_ni_if / (static_cast<Real>(cnt_if) + TinyReal);
        std::cout << std::scientific << std::setprecision(6) << "[em-freq-dual-verify-field-cross] stage=" << stage_label
                  << " body=" << body_tag << " n=" << n << " |B_real_split-B_real_merged| mean=" << mean_dr
                  << " max=" << max_dr << " if_mean=" << mean_dr_if << " if_max=" << max_dr_if << " if_cnt=" << cnt_if
                  << " |B_imag| mean=" << mean_di << " max=" << max_di << " if_mean=" << mean_di_if
                  << " if_max=" << max_di_if;
        if (include_curl_nu_b)
        {
            std::cout << " |CurlNuB_real| mean=" << mean_nr << " max=" << max_nr << " if_mean=" << mean_nr_if
                      << " if_max=" << max_nr_if << " |CurlNuB_imag| mean=" << mean_ni << " max=" << max_ni
                      << " if_mean=" << mean_ni_if << " if_max=" << max_ni_if;
        }
        else
        {
            std::cout << " |CurlNuB| skipped";
        }
        std::cout << std::endl;
    };

    scan_body(left_particles, "LeftBody");
    scan_body(right_particles, "RightBody");
}

void report_split_merged_position_consistency(BaseParticles &left_particles,
                                                BaseParticles &right_particles,
                                                BaseParticles &merged_particles)
{
    Vecd *pos_l = left_particles.getVariableDataByName<Vecd>("Position");
    Vecd *pos_r = right_particles.getVariableDataByName<Vecd>("Position");
    Vecd *pos_m = merged_particles.getVariableDataByName<Vecd>("Position");
    const size_t nl = left_particles.TotalRealParticles();
    const size_t nr = right_particles.TotalRealParticles();
    const size_t nm = merged_particles.TotalRealParticles();
    if (nl + nr != nm)
    {
        std::cout << "[em-freq-dual-verify-position] WARNING split total " << (nl + nr)
                  << " != merged " << nm << std::endl;
    }
    Real max_gap = 0.0;
    Real sum_gap = 0.0;
    size_t counted = 0;
    auto accumulate_gap = [&](const Vec3d &p)
    {
        Real best = std::numeric_limits<Real>::max();
        for (size_t j = 0; j < nm; ++j)
        {
            Vec3d q(pos_m[j][0], pos_m[j][1], pos_m[j][2]);
            best = SMIN(best, (p - q).norm());
        }
        max_gap = SMAX(max_gap, best);
        sum_gap += best;
        ++counted;
    };
    for (size_t i = 0; i < nl; ++i)
    {
        accumulate_gap(Vec3d(pos_l[i][0], pos_l[i][1], pos_l[i][2]));
    }
    for (size_t i = 0; i < nr; ++i)
    {
        accumulate_gap(Vec3d(pos_r[i][0], pos_r[i][1], pos_r[i][2]));
    }
    const Real mean_gap = sum_gap / (static_cast<Real>(counted) + TinyReal);
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-dual-verify-position] split->merged nearest max_gap=" << max_gap
              << ", mean_gap=" << mean_gap << ", split_particles=" << counted << ", merged_particles=" << nm
              << std::endl;
}

void apply_vec_field(BaseParticles &particles,
                     const std::string &field_name,
                     const std::vector<Vecd> &values)
{
    Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    size_t total_real_particles = particles.TotalRealParticles();
    size_t n = SMIN(total_real_particles, values.size());
    for (size_t i = 0; i != n; ++i)
    {
        field[i] = values[i];
    }
}

void apply_contact_correction_to_curl_nub(BaseParticles &particles,
                                          const std::string &field_name,
                                          const std::vector<Vecd> &inner_snapshot,
                                          const std::string &mode,
                                          Real damping_factor,
                                          Real limit_factor)
{
    if (mode == "baseline")
    {
        return;
    }
    Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    size_t total_real_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vecd inner_value = i < inner_snapshot.size() ? inner_snapshot[i] : ZeroData<Vecd>::value;
        Vecd contact_delta = field[i] - inner_value;
        Vecd corrected_delta = contact_delta;
        if (mode == "symmetric_damped")
        {
            corrected_delta = damping_factor * contact_delta;
        }
        else if (mode == "limited_contact")
        {
            Real inner_norm = inner_value.norm();
            Real max_contact_norm = limit_factor * inner_norm;
            Real contact_norm = contact_delta.norm();
            if (contact_norm > TinyReal && contact_norm > max_contact_norm)
            {
                corrected_delta *= (max_contact_norm / contact_norm);
            }
        }
        field[i] = inner_value + corrected_delta;
    }
}

void write_particle_diagnostics(const std::string &file_path,
                                const std::string &body_name,
                                BaseParticles &particles,
                                const Vec3d &target_b_real,
                                const Vec3d &target_b_imag)
{
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    file << std::setprecision(12);

    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *curl_real =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
    AngularVecd *curl_imag =
        particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
    Vecd *curl_nu_b_real = particles.getVariableDataByName<Vecd>("CurlNuBReal");
    Vecd *curl_nu_b_imag = particles.getVariableDataByName<Vecd>("CurlNuBImag");

    size_t total_real_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d discrete_b_real = angular_to_vec(curl_real[i]);
        Vec3d discrete_b_imag = angular_to_vec(curl_imag[i]);
        Real b_real_error = (discrete_b_real - target_b_real).norm();
        Real b_imag_error = (discrete_b_imag - target_b_imag).norm();
        Real curl_nu_b_real_norm = curl_nu_b_real[i].norm();
        Real curl_nu_b_imag_norm = curl_nu_b_imag[i].norm();
        file << body_name << ","
             << i << ","
             << positions[i][0] << ","
             << positions[i][1] << ","
             << positions[i][2] << ","
             << fabs(positions[i][0]) << ","
             << static_cast<int>(is_interface_particle(positions[i])) << ","
             << discrete_b_real[0] << ","
             << discrete_b_real[1] << ","
             << discrete_b_real[2] << ","
             << b_real_error << ","
             << discrete_b_imag[0] << ","
             << discrete_b_imag[1] << ","
             << discrete_b_imag[2] << ","
             << b_imag_error << ","
             << curl_nu_b_real[i][0] << ","
             << curl_nu_b_real[i][1] << ","
             << curl_nu_b_real[i][2] << ","
             << curl_nu_b_real_norm << ","
             << curl_nu_b_imag[i][0] << ","
             << curl_nu_b_imag[i][1] << ","
             << curl_nu_b_imag[i][2] << ","
             << curl_nu_b_imag_norm << "\n";
    }
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_FREQ_DUAL_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-dual-verify-config] dp=" << dp_0
              << ", total_length=" << total_length
              << ", sigma=" << conductivity
              << ", rho_cp=" << rho_cp
              << ", use_inner=" << use_inner_operator
              << ", use_contact=" << use_contact_operator
              << ", curl_scaling=" << curl_scaling
              << ", curl_nu_b_scaling=" << curl_nu_b_scaling
              << ", curl_nu_b_symmetric_volume_weight=" << curl_nu_b_symmetric_volume_weight
              << ", neighbor_stencil_mode=" << neighbor_stencil_mode
              << ", position_consistency_diag=" << position_consistency_diag
              << ", split_merged_field_diag=" << split_merged_field_diag
              << ", minimal_io_mode=" << minimal_io_mode
              << ", contact_correction_mode=" << contact_correction_mode
              << ", contact_damping_factor=" << contact_damping_factor
              << ", contact_limit_factor=" << contact_limit_factor
              << ", B_real=(" << target_b_real[0] << ","
              << target_b_real[1] << "," << target_b_real[2] << ")"
              << ", B_imag=(" << target_b_imag[0] << ","
              << target_b_imag[1] << "," << target_b_imag[2] << ")"
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

    if (position_consistency_diag)
    {
        report_split_merged_position_consistency(left_body.getBaseParticles(), right_body.getBaseParticles(),
                                                 merged_body.getBaseParticles());
    }

    const bool use_inner_contact_stencil = (neighbor_stencil_mode & 2) != 0;
    const bool use_contact_inner_stencil = (neighbor_stencil_mode & 1) != 0;

    std::unique_ptr<InnerRelation> left_inner_holder;
    std::unique_ptr<InnerRelation> right_inner_holder;
    std::unique_ptr<InnerRelation> merged_inner_holder;
    if (use_inner_contact_stencil)
    {
        left_inner_holder =
            std::make_unique<electromagnetics::InnerRelationContactSearchDepth>(left_body);
        right_inner_holder =
            std::make_unique<electromagnetics::InnerRelationContactSearchDepth>(right_body);
        merged_inner_holder =
            std::make_unique<electromagnetics::InnerRelationContactSearchDepth>(merged_body);
    }
    else
    {
        left_inner_holder = std::make_unique<InnerRelation>(left_body);
        right_inner_holder = std::make_unique<InnerRelation>(right_body);
        merged_inner_holder = std::make_unique<InnerRelation>(merged_body);
    }
    InnerRelation &left_inner = *left_inner_holder;
    InnerRelation &right_inner = *right_inner_holder;
    InnerRelation &merged_inner = *merged_inner_holder;

    std::unique_ptr<ContactRelation> left_contact_holder;
    std::unique_ptr<ContactRelation> right_contact_holder;
    if (use_contact_inner_stencil)
    {
        left_contact_holder =
            std::make_unique<electromagnetics::ContactRelationInnerSearchDepth>(left_body, RealBodyVector{&right_body});
        right_contact_holder =
            std::make_unique<electromagnetics::ContactRelationInnerSearchDepth>(right_body, RealBodyVector{&left_body});
    }
    else
    {
        left_contact_holder = std::make_unique<ContactRelation>(left_body, RealBodyVector{&right_body});
        right_contact_holder = std::make_unique<ContactRelation>(right_body, RealBodyVector{&left_body});
    }
    ContactRelation &left_contact = *left_contact_holder;
    ContactRelation &right_contact = *right_contact_holder;

    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_left_em(left_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_right_em(right_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_merged_em(merged_body, conductivity, rho_cp, magnetic_reluctivity);
    SimpleDynamics<AssignLinearFrequencyVectorPotential>
        assign_left_linear_a(left_body, target_b_real, target_b_imag, gauge_origin);
    SimpleDynamics<AssignLinearFrequencyVectorPotential>
        assign_right_linear_a(right_body, target_b_real, target_b_imag, gauge_origin);
    SimpleDynamics<AssignLinearFrequencyVectorPotential>
        assign_merged_linear_a(merged_body, target_b_real, target_b_imag, gauge_origin);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_real_left_inner(left_inner, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_real_right_inner(right_inner, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_real_merged_inner(merged_inner, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_real_left_contact(left_contact, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_real_right_contact(right_contact, "VectorPotentialReal", "VectorPotentialCurlReal", curl_scaling);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_imag_left_inner(left_inner, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_imag_right_inner(right_inner, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a_imag_merged_inner(merged_inner, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_imag_left_contact(left_contact, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        curl_a_imag_right_contact(right_contact, "VectorPotentialImag", "VectorPotentialCurlImag", curl_scaling);

    /** Curl(nu B) must use the final B = curl(A) after contact on split bodies. */
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_left_inner(left_inner, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling,
                                  curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_right_inner(right_inner, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling,
                                   curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_merged_inner(merged_inner, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling,
                                    curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_left_contact(left_contact, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling,
                                    curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_right_contact(right_contact, "VectorPotentialCurlReal", "CurlNuBReal", curl_nu_b_scaling,
                                     curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_left_inner(left_inner, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling,
                                  curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_right_inner(right_inner, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling,
                                   curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_merged_inner(merged_inner, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling,
                                    curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_left_contact(left_contact, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling,
                                    curl_nu_b_symmetric_volume_weight);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_right_contact(right_contact, "VectorPotentialCurlImag", "CurlNuBImag", curl_nu_b_scaling,
                                     curl_nu_b_symmetric_volume_weight);

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
        curl_a_real_left_inner.exec();
        curl_a_real_right_inner.exec();
        curl_a_real_merged_inner.exec();
        curl_a_imag_left_inner.exec();
        curl_a_imag_right_inner.exec();
        curl_a_imag_merged_inner.exec();
        if (split_merged_field_diag)
        {
            report_split_vs_merged_discrete_fields("curlA_inner_only_split_vs_merged_inner", left_body.getBaseParticles(),
                                                   right_body.getBaseParticles(), merged_body.getBaseParticles(), false);
        }
        if (use_contact_operator)
        {
            curl_a_real_left_contact.exec();
            curl_a_real_right_contact.exec();
            curl_a_imag_left_contact.exec();
            curl_a_imag_right_contact.exec();
            if (split_merged_field_diag)
            {
                report_split_vs_merged_discrete_fields(
                    "curlA_split_inner_plus_contact_vs_merged_inner", left_body.getBaseParticles(),
                    right_body.getBaseParticles(), merged_body.getBaseParticles(), false);
            }
        }
        curl_nu_b_real_left_inner.exec();
        curl_nu_b_real_right_inner.exec();
        curl_nu_b_real_merged_inner.exec();
        curl_nu_b_imag_left_inner.exec();
        curl_nu_b_imag_right_inner.exec();
        curl_nu_b_imag_merged_inner.exec();
        if (split_merged_field_diag)
        {
            report_split_vs_merged_discrete_fields("CurlNuB_inner_using_final_B", left_body.getBaseParticles(),
                                                   right_body.getBaseParticles(), merged_body.getBaseParticles(), true);
        }
    }
    std::vector<Vecd> left_curl_nub_real_inner =
        snapshot_vec_field(left_body.getBaseParticles(), "CurlNuBReal");
    std::vector<Vecd> left_curl_nub_imag_inner =
        snapshot_vec_field(left_body.getBaseParticles(), "CurlNuBImag");
    std::vector<Vecd> right_curl_nub_real_inner =
        snapshot_vec_field(right_body.getBaseParticles(), "CurlNuBReal");
    std::vector<Vecd> right_curl_nub_imag_inner =
        snapshot_vec_field(right_body.getBaseParticles(), "CurlNuBImag");
    if (use_inner_operator && use_contact_operator)
    {
        curl_nu_b_real_left_contact.exec();
        curl_nu_b_real_right_contact.exec();
        curl_nu_b_imag_left_contact.exec();
        curl_nu_b_imag_right_contact.exec();
        apply_contact_correction_to_curl_nub(left_body.getBaseParticles(), "CurlNuBReal",
                                             left_curl_nub_real_inner, contact_correction_mode,
                                             contact_damping_factor, contact_limit_factor);
        apply_contact_correction_to_curl_nub(left_body.getBaseParticles(), "CurlNuBImag",
                                             left_curl_nub_imag_inner, contact_correction_mode,
                                             contact_damping_factor, contact_limit_factor);
        apply_contact_correction_to_curl_nub(right_body.getBaseParticles(), "CurlNuBReal",
                                             right_curl_nub_real_inner, contact_correction_mode,
                                             contact_damping_factor, contact_limit_factor);
        apply_contact_correction_to_curl_nub(right_body.getBaseParticles(), "CurlNuBImag",
                                             right_curl_nub_imag_inner, contact_correction_mode,
                                             contact_damping_factor, contact_limit_factor);
        if (split_merged_field_diag)
        {
            report_split_vs_merged_discrete_fields("after_CurlNuB_contact_and_optional_correction",
                                                   left_body.getBaseParticles(), right_body.getBaseParticles(),
                                                   merged_body.getBaseParticles(), true);
        }
    }

    if (!minimal_io_mode)
    {
        write_states.writeToFile(0.0);
    }

    BodyCurlSummary left_summary = evaluate_body_curl_summary(
        "LeftBody", left_body.getBaseParticles(), target_b_real, target_b_imag);
    BodyCurlSummary right_summary = evaluate_body_curl_summary(
        "RightBody", right_body.getBaseParticles(), target_b_real, target_b_imag);
    BodyCurlSummary merged_summary = evaluate_body_curl_summary(
        "MergedBody", merged_body.getBaseParticles(), target_b_real, target_b_imag);

    if (!minimal_io_mode)
    {
        std::ofstream summary_file(
            io_environment.OutputFolder() + "/em_frequency_dual_body_linear_a_summary.csv",
            std::ios::out | std::ios::trunc);
        summary_file << std::setprecision(12);
        summary_file
            << "body,total_particles,interface_particles,core_particles,"
            << "mean_b_real_error,max_b_real_error,mean_b_real_error_interface,max_b_real_error_interface,"
            << "mean_b_real_error_core,max_b_real_error_core,"
            << "mean_b_imag_error,max_b_imag_error,mean_b_imag_error_interface,max_b_imag_error_interface,"
            << "mean_b_imag_error_core,max_b_imag_error_core,"
            << "mean_curl_nu_b_real_norm,max_curl_nu_b_real_norm,mean_curl_nu_b_real_norm_interface,max_curl_nu_b_real_norm_interface,"
            << "mean_curl_nu_b_real_norm_core,max_curl_nu_b_real_norm_core,"
            << "mean_curl_nu_b_imag_norm,max_curl_nu_b_imag_norm,mean_curl_nu_b_imag_norm_interface,max_curl_nu_b_imag_norm_interface,"
            << "mean_curl_nu_b_imag_norm_core,max_curl_nu_b_imag_norm_core\n";

        auto write_summary = [&](const BodyCurlSummary &summary)
        {
            summary_file << summary.body_name << ","
                         << summary.total_particles << ","
                         << summary.interface_particles << ","
                         << summary.core_particles << ","
                         << summary.mean_b_real_error << ","
                         << summary.max_b_real_error << ","
                         << summary.mean_b_real_error_interface << ","
                         << summary.max_b_real_error_interface << ","
                         << summary.mean_b_real_error_core << ","
                         << summary.max_b_real_error_core << ","
                         << summary.mean_b_imag_error << ","
                         << summary.max_b_imag_error << ","
                         << summary.mean_b_imag_error_interface << ","
                         << summary.max_b_imag_error_interface << ","
                         << summary.mean_b_imag_error_core << ","
                         << summary.max_b_imag_error_core << ","
                         << summary.mean_curl_nu_b_real_norm << ","
                         << summary.max_curl_nu_b_real_norm << ","
                         << summary.mean_curl_nu_b_real_norm_interface << ","
                         << summary.max_curl_nu_b_real_norm_interface << ","
                         << summary.mean_curl_nu_b_real_norm_core << ","
                         << summary.max_curl_nu_b_real_norm_core << ","
                         << summary.mean_curl_nu_b_imag_norm << ","
                         << summary.max_curl_nu_b_imag_norm << ","
                         << summary.mean_curl_nu_b_imag_norm_interface << ","
                         << summary.max_curl_nu_b_imag_norm_interface << ","
                         << summary.mean_curl_nu_b_imag_norm_core << ","
                         << summary.max_curl_nu_b_imag_norm_core << "\n";
        };
        write_summary(left_summary);
        write_summary(right_summary);
        write_summary(merged_summary);
        summary_file.flush();

        std::ofstream particle_file(
            io_environment.OutputFolder() + "/em_frequency_dual_body_linear_a_particles.csv",
            std::ios::out | std::ios::trunc);
        particle_file << std::setprecision(12);
        particle_file
            << "body,particle_id,x,y,z,abs_interface_distance,is_interface,"
            << "curl_a_real_x,curl_a_real_y,curl_a_real_z,b_real_error,"
            << "curl_a_imag_x,curl_a_imag_y,curl_a_imag_z,b_imag_error,"
            << "curl_nu_b_real_x,curl_nu_b_real_y,curl_nu_b_real_z,curl_nu_b_real_norm,"
            << "curl_nu_b_imag_x,curl_nu_b_imag_y,curl_nu_b_imag_z,curl_nu_b_imag_norm\n";
        particle_file.close();

        write_particle_diagnostics(
            io_environment.OutputFolder() + "/em_frequency_dual_body_linear_a_particles.csv",
            "LeftBody", left_body.getBaseParticles(), target_b_real, target_b_imag);
        write_particle_diagnostics(
            io_environment.OutputFolder() + "/em_frequency_dual_body_linear_a_particles.csv",
            "RightBody", right_body.getBaseParticles(), target_b_real, target_b_imag);
        write_particle_diagnostics(
            io_environment.OutputFolder() + "/em_frequency_dual_body_linear_a_particles.csv",
            "MergedBody", merged_body.getBaseParticles(), target_b_real, target_b_imag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-dual-verify-summary] LeftBody mean|curlA_real-B_real|="
              << left_summary.mean_b_real_error
              << ", mean|curlA_imag-B_imag|=" << left_summary.mean_b_imag_error
              << ", mean|curlNuB_real|=" << left_summary.mean_curl_nu_b_real_norm
              << ", mean|curlNuB_imag|=" << left_summary.mean_curl_nu_b_imag_norm
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-dual-verify-summary] RightBody mean|curlA_real-B_real|="
              << right_summary.mean_b_real_error
              << ", mean|curlA_imag-B_imag|=" << right_summary.mean_b_imag_error
              << ", mean|curlNuB_real|=" << right_summary.mean_curl_nu_b_real_norm
              << ", mean|curlNuB_imag|=" << right_summary.mean_curl_nu_b_imag_norm
              << std::endl;
    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-dual-verify-summary] MergedBody mean|curlA_real-B_real|="
              << merged_summary.mean_b_real_error
              << ", mean|curlA_imag-B_imag|=" << merged_summary.mean_b_imag_error
              << ", mean|curlNuB_real|=" << merged_summary.mean_curl_nu_b_real_norm
              << ", mean|curlNuB_imag|=" << merged_summary.mean_curl_nu_b_imag_norm
              << std::endl;
    std::cout << "[em-freq-dual-verify] summary file: "
              << io_environment.OutputFolder() +
                     "/em_frequency_dual_body_linear_a_summary.csv"
              << std::endl;

    return 0;
}
