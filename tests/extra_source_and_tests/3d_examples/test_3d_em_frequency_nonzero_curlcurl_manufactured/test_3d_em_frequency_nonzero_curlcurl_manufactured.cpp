/**
 * @file test_3d_em_frequency_nonzero_curlcurl_manufactured.cpp
 * @brief Manufactured non-zero curl(nu curl A) verification for the SPH A-phi chain.
 *
 * Analytical field:
 *   A = (0, 0, A0 sin(kx x) sin(ky y))
 *   B = curl(A) = (dA_z/dy, -dA_z/dx, 0)
 *   curl(nu B) = nu * (kx^2 + ky^2) * A_z * e_z, for constant nu.
 *
 * This complements the existing linear-A tests, where curl(nu curl A)=0.
 */

#include "sphinxsys.h"
#include "electromagnetic_team7_aphi_dynamics.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace SPH;

namespace
{
Real get_env_real_local(const std::string &name, Real default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr) return default_value;
    char *end_ptr = nullptr;
    Real parsed = static_cast<Real>(std::strtod(value, &end_ptr));
    if (end_ptr == value || !std::isfinite(parsed)) return default_value;
    return parsed;
}

std::string get_env_string_local(const std::string &name,
                                 const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

const Real dp_0 = get_env_real_local("EM_NONZERO_CURLCURL_DP", 0.5);
const Real lx = get_env_real_local("EM_NONZERO_CURLCURL_LX", 8.0);
const Real ly = get_env_real_local("EM_NONZERO_CURLCURL_LY", 8.0);
const Real lz = get_env_real_local("EM_NONZERO_CURLCURL_LZ", 4.0);
const Real boundary_width = get_env_real_local("EM_NONZERO_CURLCURL_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real boundary_shell = get_env_real_local("EM_NONZERO_CURLCURL_BOUNDARY_SHELL", 2.5 * dp_0);
const Real amplitude = get_env_real_local("EM_NONZERO_CURLCURL_A0", 1.0);
const Real mode_x = get_env_real_local("EM_NONZERO_CURLCURL_MODE_X", 1.0);
const Real mode_y = get_env_real_local("EM_NONZERO_CURLCURL_MODE_Y", 1.0);
const Real magnetic_reluctivity = get_env_real_local(
    "EM_NONZERO_CURLCURL_NU", 1.0 / (4.0 * Pi * 1.0e-7));

const Real kx = mode_x * Pi / lx;
const Real ky = mode_y * Pi / ly;
const Vec3d body_center(0.5 * lx, 0.5 * ly, 0.5 * lz);
const Vec3d body_halfsize(0.5 * lx, 0.5 * ly, 0.5 * lz);

BoundingBoxd system_domain_bounds(
    Vec3d(-boundary_width, -boundary_width, -boundary_width),
    Vec3d(lx + boundary_width, ly + boundary_width, lz + boundary_width));

class BoxShape : public ComplexShape
{
  public:
    explicit BoxShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(body_center), body_halfsize);
    }
};

Real az_exact(const Vec3d &p)
{
    return amplitude * std::sin(kx * p[0]) * std::sin(ky * p[1]);
}

Vec3d a_exact(const Vec3d &p)
{
    return Vec3d(0.0, 0.0, az_exact(p));
}

Vec3d b_exact(const Vec3d &p)
{
    Real bx = amplitude * ky * std::sin(kx * p[0]) * std::cos(ky * p[1]);
    Real by = -amplitude * kx * std::cos(kx * p[0]) * std::sin(ky * p[1]);
    return Vec3d(bx, by, 0.0);
}

Vec3d curl_nu_b_exact(const Vec3d &p)
{
    return Vec3d(0.0, 0.0,
                 magnetic_reluctivity * (kx * kx + ky * ky) * az_exact(p));
}

bool is_boundary_shell_particle(const Vec3d &p)
{
    Real distance = std::min({p[0], lx - p[0], p[1], ly - p[1], p[2], lz - p[2]});
    return distance <= boundary_shell;
}

Vec3d angular_to_vec(const AngularVecd &v)
{
    return Vec3d(v[0], v[1], v[2]);
}

class AssignSinusoidalVectorPotential : public LocalDynamics
{
  public:
    explicit AssignSinusoidalVectorPotential(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          a_(particles_->getVariableDataByName<Vecd>("VectorPotential")),
          a_prev_(particles_->getVariableDataByName<Vecd>("VectorPotentialPrevious")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d value = a_exact(pos_[index_i]);
        a_[index_i] = value;
        a_prev_[index_i] = value;
    }

  protected:
    Vecd *pos_, *a_, *a_prev_;
};

struct Summary
{
    size_t total = 0, core = 0, boundary = 0;
    Real mean_b_error = 0.0, max_b_error = 0.0;
    Real mean_b_error_core = 0.0, max_b_error_core = 0.0;
    Real mean_b_error_boundary = 0.0, max_b_error_boundary = 0.0;
    Real mean_curl_nu_b_error = 0.0, max_curl_nu_b_error = 0.0;
    Real mean_curl_nu_b_error_core = 0.0, max_curl_nu_b_error_core = 0.0;
    Real mean_curl_nu_b_error_boundary = 0.0, max_curl_nu_b_error_boundary = 0.0;
};

Summary evaluate_summary(BaseParticles &particles)
{
    Summary s;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *curl_a = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>("CurlNuB");

    Real sum_b = 0.0, sum_b_core = 0.0, sum_b_boundary = 0.0;
    Real sum_cnb = 0.0, sum_cnb_core = 0.0, sum_cnb_boundary = 0.0;
    size_t n_real = particles.TotalRealParticles();
    for (size_t i = 0; i != n_real; ++i)
    {
        Vec3d p = pos[i];
        Real b_err = (angular_to_vec(curl_a[i]) - b_exact(p)).norm();
        Real cnb_err = (curl_nu_b[i] - curl_nu_b_exact(p)).norm();
        bool boundary = is_boundary_shell_particle(p);

        s.total++;
        sum_b += b_err;
        sum_cnb += cnb_err;
        s.max_b_error = SMAX(s.max_b_error, b_err);
        s.max_curl_nu_b_error = SMAX(s.max_curl_nu_b_error, cnb_err);

        if (boundary)
        {
            s.boundary++;
            sum_b_boundary += b_err;
            sum_cnb_boundary += cnb_err;
            s.max_b_error_boundary = SMAX(s.max_b_error_boundary, b_err);
            s.max_curl_nu_b_error_boundary = SMAX(s.max_curl_nu_b_error_boundary, cnb_err);
        }
        else
        {
            s.core++;
            sum_b_core += b_err;
            sum_cnb_core += cnb_err;
            s.max_b_error_core = SMAX(s.max_b_error_core, b_err);
            s.max_curl_nu_b_error_core = SMAX(s.max_curl_nu_b_error_core, cnb_err);
        }
    }
    s.mean_b_error = sum_b / (static_cast<Real>(s.total) + TinyReal);
    s.mean_curl_nu_b_error = sum_cnb / (static_cast<Real>(s.total) + TinyReal);
    s.mean_b_error_core = sum_b_core / (static_cast<Real>(s.core) + TinyReal);
    s.mean_curl_nu_b_error_core = sum_cnb_core / (static_cast<Real>(s.core) + TinyReal);
    s.mean_b_error_boundary = sum_b_boundary / (static_cast<Real>(s.boundary) + TinyReal);
    s.mean_curl_nu_b_error_boundary = sum_cnb_boundary / (static_cast<Real>(s.boundary) + TinyReal);
    return s;
}

void write_particles(const std::string &path, BaseParticles &particles)
{
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    f << std::setprecision(12);
    f << "particle_id,x,y,z,is_boundary,"
      << "b_x,b_y,b_z,b_exact_x,b_exact_y,b_exact_z,b_error,"
      << "curl_nu_b_x,curl_nu_b_y,curl_nu_b_z,"
      << "curl_nu_b_exact_x,curl_nu_b_exact_y,curl_nu_b_exact_z,curl_nu_b_error\n";

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    AngularVecd *curl_a = particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>("CurlNuB");
    size_t n_real = particles.TotalRealParticles();
    for (size_t i = 0; i != n_real; ++i)
    {
        Vec3d p = pos[i];
        Vec3d b = angular_to_vec(curl_a[i]);
        Vec3d be = b_exact(p);
        Vec3d cnbe = curl_nu_b_exact(p);
        f << i << "," << p[0] << "," << p[1] << "," << p[2] << ","
          << static_cast<int>(is_boundary_shell_particle(p)) << ","
          << b[0] << "," << b[1] << "," << b[2] << ","
          << be[0] << "," << be[1] << "," << be[2] << ","
          << (b - be).norm() << ","
          << curl_nu_b[i][0] << "," << curl_nu_b[i][1] << "," << curl_nu_b[i][2] << ","
          << cnbe[0] << "," << cnbe[1] << "," << cnbe[2] << ","
          << (curl_nu_b[i] - cnbe).norm() << "\n";
    }
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_NONZERO_CURLCURL_OUTPUT_TAG", "");
    if (!output_tag.empty()) io_environment.appendOutputFolder(output_tag);

    SolidBody body(sph_system, makeShared<BoxShape>("NonzeroCurlCurlBody"));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner(body);

    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_em(body, 0.0, 1.0, magnetic_reluctivity);
    SimpleDynamics<AssignSinusoidalVectorPotential> assign_a(body);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        curl_a(inner, "VectorPotential", "VectorPotentialCurl", 1.0);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b(inner, "VectorPotentialCurl", "CurlNuB", 1.0);
    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_em.exec();
    assign_a.exec();
    curl_a.exec();
    curl_nu_b.exec();
    write_states.writeToFile(0.0);

    Summary s = evaluate_summary(body.getBaseParticles());
    std::ofstream summary(io_environment.OutputFolder() +
                              "/em_frequency_nonzero_curlcurl_manufactured_summary.csv",
                          std::ios::out | std::ios::trunc);
    summary << std::setprecision(12);
    summary << "total_particles,core_particles,boundary_particles,"
            << "mean_b_error,max_b_error,mean_b_error_core,max_b_error_core,"
            << "mean_b_error_boundary,max_b_error_boundary,"
            << "mean_curl_nu_b_error,max_curl_nu_b_error,"
            << "mean_curl_nu_b_error_core,max_curl_nu_b_error_core,"
            << "mean_curl_nu_b_error_boundary,max_curl_nu_b_error_boundary\n";
    summary << s.total << "," << s.core << "," << s.boundary << ","
            << s.mean_b_error << "," << s.max_b_error << ","
            << s.mean_b_error_core << "," << s.max_b_error_core << ","
            << s.mean_b_error_boundary << "," << s.max_b_error_boundary << ","
            << s.mean_curl_nu_b_error << "," << s.max_curl_nu_b_error << ","
            << s.mean_curl_nu_b_error_core << "," << s.max_curl_nu_b_error_core << ","
            << s.mean_curl_nu_b_error_boundary << "," << s.max_curl_nu_b_error_boundary << "\n";
    summary.flush();

    write_particles(io_environment.OutputFolder() +
                        "/em_frequency_nonzero_curlcurl_manufactured_particles.csv",
                    body.getBaseParticles());

    std::cout << std::scientific << std::setprecision(6)
              << "[nonzero-curlcurl] mean_b_error=" << s.mean_b_error
              << ", mean_b_error_core=" << s.mean_b_error_core
              << ", mean_curl_nu_b_error=" << s.mean_curl_nu_b_error
              << ", mean_curl_nu_b_error_core=" << s.mean_curl_nu_b_error_core
              << std::endl;
    return 0;
}
