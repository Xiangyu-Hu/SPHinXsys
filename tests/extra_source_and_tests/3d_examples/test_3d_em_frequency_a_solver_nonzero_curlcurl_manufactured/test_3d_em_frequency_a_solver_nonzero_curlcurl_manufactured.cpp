/**
 * @file test_3d_em_frequency_a_solver_nonzero_curlcurl_manufactured.cpp
 * @brief Manufactured verification for the frequency-domain A solver with nonzero curl(nu curl A).
 *
 * The analytical field is chosen as a z-only sinusoidal vector potential:
 *
 *   A_z(x, y) = amplitude * sin(kx * x) * sin(ky * y)
 *
 * so that div(A) = 0 and
 *
 *   curl(nu curl(A)) = nu * (kx^2 + ky^2) * A_z * e_z
 *
 * is nonzero. This exercises the main magnetic operator directly, unlike the
 * earlier linear-A manufactured case where curl(nu curl(A)) = 0.
 */

#include "sphinxsys.h"
#include "electromagnetic_aphi_global_implicit_solver.hpp"
#include "electromagnetic_component_hessian_ck.hpp"
#include "electromagnetic_team7_aphi_dynamics.hpp"
#include "electromagnetic_team7_aphi_frequency_dynamics.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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

std::string to_lower_copy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

bool solver_mode_contains(const std::string &solver_mode, const std::string &token)
{
    return solver_mode == token || solver_mode == "both" || solver_mode == "all";
}

const Real dp_0 = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_DP", 1.0);
const Real body_length = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_LENGTH", 20.0);
const Real body_height = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_WIDTH", 8.0);
const Real boundary_width = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_BOUNDARY_WIDTH", 3.0 * dp_0);
const Real boundary_shell_thickness =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_BOUNDARY_SHELL_THICKNESS", 2.5 * dp_0);
const Real conductivity = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_SIGMA", 3.0e6);
const Real rho_cp = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_RHO_CP", 1.0);
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);
const Real frequency_hz = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_FREQUENCY_HZ", 50.0);
const Real omega = 2.0 * Pi * frequency_hz;
const Real dt_pseudo = get_env_real_local("EM_FREQ_CURLCURL_VERIFY_DT", 1.0e-6);
const Real reference_dt_pseudo =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_REFERENCE_DT", dt_pseudo);
const int solve_iterations =
    static_cast<int>(get_env_real_local("EM_FREQ_CURLCURL_VERIFY_ITERATIONS", 2000.0));
const Real sigma_relaxation_scaling =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_SIGMA_RELAX_SCALE", 1.0);
const Real sigma_relaxation_floor =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_SIGMA_RELAX_FLOOR", TinyReal);
const Real magnetic_diagonal_scaling =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_MAG_DIAG_SCALE", 1.0);
const Real relaxation_scaling =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_RELAXATION_SCALING", 1.0);
const Real max_change_rate =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_MAX_CHANGE_RATE", 1.0e6);
const Real curl_scaling =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_CURL_SCALING", 1.0);
const Real curl_nu_b_scaling =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_CURL_NUB_SCALING", 1.0);
const Real initial_guess_scale_real =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_INITIAL_GUESS_SCALE_REAL", 0.0);
const Real initial_guess_scale_imag =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_INITIAL_GUESS_SCALE_IMAG", 0.0);
const bool write_vtp = get_env_bool_local("EM_FREQ_CURLCURL_VERIFY_WRITE_VTP", true);
const bool write_particles = get_env_bool_local("EM_FREQ_CURLCURL_VERIFY_WRITE_PARTICLES", true);
const bool use_boundary_constraint =
    get_env_bool_local("EM_FREQ_CURLCURL_VERIFY_USE_BOUNDARY_CONSTRAINT", true);
const std::string solver_mode_raw =
    get_env_string_local("EM_FREQ_CURLCURL_VERIFY_SOLVER_MODE", "coupled");
const std::string operator_mode_raw =
    get_env_string_local("EM_FREQ_CURLCURL_VERIFY_OPERATOR_MODE", "direct");
const std::string global_alpha_formula_raw =
    get_env_string_local("EM_FREQ_CURLCURL_VERIFY_GLOBAL_ALPHA_FORMULA", "preconditioned_sd");
const Real global_alpha_max =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GLOBAL_ALPHA_MAX", 1.0);
const int history_interval =
    static_cast<int>(get_env_real_local("EM_FREQ_CURLCURL_VERIFY_HISTORY_INTERVAL", 50.0));

const Real amplitude_real =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_AMPLITUDE_REAL", 1.0);
const Real amplitude_imag =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_AMPLITUDE_IMAG", 0.5);
const Real kx_multiplier =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_KX_MULTIPLIER", 1.0);
const Real ky_multiplier =
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_KY_MULTIPLIER", 1.0);

const Vec3d target_grad_phi_real(
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_REAL_X", 0.0),
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_REAL_Y", 0.0),
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_REAL_Z", 0.0));
const Vec3d target_grad_phi_imag(
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_IMAG_X", 0.0),
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_IMAG_Y", 0.0),
    get_env_real_local("EM_FREQ_CURLCURL_VERIFY_GRADPHI_IMAG_Z", 0.0));

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

struct SinusoidalAzField
{
    Real amplitude;
    Real kx;
    Real ky;
    Real magnetic_reluctivity;

    Vec3d evaluateA(const Vec3d &position) const
    {
        Real shape = amplitude * sin(kx * position[0]) * sin(ky * position[1]);
        return Vec3d(0.0, 0.0, shape);
    }

    Vec3d evaluateB(const Vec3d &position) const
    {
        Real sin_x = sin(kx * position[0]);
        Real cos_x = cos(kx * position[0]);
        Real sin_y = sin(ky * position[1]);
        Real cos_y = cos(ky * position[1]);
        return Vec3d(amplitude * ky * sin_x * cos_y,
                     -amplitude * kx * cos_x * sin_y,
                     0.0);
    }

    Vec3d evaluateCurlNuB(const Vec3d &position) const
    {
        Vec3d a = evaluateA(position);
        return magnetic_reluctivity * (kx * kx + ky * ky) * a;
    }
};

const Real kx = kx_multiplier * Pi / body_length;
const Real ky = ky_multiplier * Pi / body_height;
const SinusoidalAzField a_real_model{amplitude_real, kx, ky, magnetic_reluctivity};
const SinusoidalAzField a_imag_model{amplitude_imag, kx, ky, magnetic_reluctivity};
const electromagnetics::ComponentVariableNames<> a_real_component_names = {
    "VerifyAxReal", "VerifyAyReal", "VerifyAzReal"};
const electromagnetics::ComponentVariableNames<> a_imag_component_names = {
    "VerifyAxImag", "VerifyAyImag", "VerifyAzImag"};

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
        const Vec3d a_real_exact = a_real_model.evaluateA(position);
        const Vec3d a_imag_exact = a_imag_model.evaluateA(position);
        const Vec3d curl_nu_b_real_exact = a_real_model.evaluateCurlNuB(position);
        const Vec3d curl_nu_b_imag_exact = a_imag_model.evaluateCurlNuB(position);
        const Real sigma_i = conductivity_[index_i];

        a_real_[index_i] = a_real_exact;
        a_imag_[index_i] = a_imag_exact;
        grad_phi_real_[index_i] = target_grad_phi_real;
        grad_phi_imag_[index_i] = target_grad_phi_imag;
        source_current_density_real_[index_i] =
            curl_nu_b_real_exact + sigma_i * target_grad_phi_real - sigma_i * omega * a_imag_exact;
        source_current_density_imag_[index_i] =
            curl_nu_b_imag_exact + sigma_i * target_grad_phi_imag + sigma_i * omega * a_real_exact;
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
                                         const SinusoidalAzField &model,
                                         Real scale)
        : LocalDynamics(sph_body),
          model_(model),
          scale_(scale),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          vector_field_(particles_->getVariableDataByName<Vecd>(field_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        vector_field_[index_i] = scale_ * model_.evaluateA(positions_[index_i]);
    }

  protected:
    SinusoidalAzField model_;
    Real scale_;
    Vecd *positions_;
    Vecd *vector_field_;
};

class ConstrainManufacturedVectorPotential : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainManufacturedVectorPotential(BodyPartByParticle &body_part,
                                                  const std::string &field_name,
                                                  const SinusoidalAzField &model)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          model_(model),
          positions_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          vector_field_(this->particles_->template getVariableDataByName<Vecd>(field_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        vector_field_[index_i] = model_.evaluateA(positions_[index_i]);
    }

  protected:
    SinusoidalAzField model_;
    Vecd *positions_;
    Vecd *vector_field_;
};

class SelectCurlNuBHybridByBoundaryShell : public LocalDynamics
{
  public:
    explicit SelectCurlNuBHybridByBoundaryShell(SPHBody &sph_body,
                                                const std::string &target_name,
                                                const std::string &core_name,
                                                const Vec3d &box_center,
                                                const Vec3d &box_halfsize,
                                                Real shell_thickness)
        : LocalDynamics(sph_body),
          box_center_(box_center),
          box_halfsize_(box_halfsize),
          shell_thickness_(shell_thickness),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          target_(particles_->getVariableDataByName<Vecd>(target_name)),
          core_(particles_->getVariableDataByName<Vecd>(core_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        Vec3d offset = positions_[index_i] - box_center_;
        Vec3d margin = box_halfsize_ - offset.cwiseAbs();
        Real distance_to_boundary = margin.minCoeff();
        bool in_shell = distance_to_boundary <= shell_thickness_ + TinyReal;
        if (!in_shell)
        {
            target_[index_i] = core_[index_i];
        }
    }

  protected:
    Vec3d box_center_;
    Vec3d box_halfsize_;
    Real shell_thickness_;
    Vecd *positions_;
    Vecd *target_;
    Vecd *core_;
};

class ExtendCurlNuBFromCoreNeighbors : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ExtendCurlNuBFromCoreNeighbors(BaseInnerRelation &inner_relation,
                                            const std::string &target_name,
                                            const std::string &core_name,
                                            Real shell_thickness)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          shell_thickness_(shell_thickness),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          target_(particles_->getVariableDataByName<Vecd>(target_name)),
          core_(particles_->getVariableDataByName<Vecd>(core_name))
    {
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        if (!is_shell_position(positions_[index_i]))
        {
            target_[index_i] = core_[index_i];
            return;
        }

        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vecd weighted_core_value = ZeroData<Vecd>::value;
        Real weight_sum = 0.0;
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (is_shell_position(positions_[index_j]))
            {
                continue;
            }
            Real weight = inner_neighborhood.W_ij_[n] * Vol_[index_j];
            if (weight <= TinyReal)
            {
                continue;
            }
            weighted_core_value += weight * core_[index_j];
            weight_sum += weight;
        }

        if (weight_sum > TinyReal)
        {
            target_[index_i] = weighted_core_value / weight_sum;
        }
    }

  protected:
    Real distance_to_boundary(const Vecd &position) const
    {
        Real dx = SMIN(position[0], body_length - position[0]);
        Real dy = SMIN(position[1], body_height - position[1]);
        Real dz = SMIN(position[2], body_width - position[2]);
        return SMAX(static_cast<Real>(0.0), SMIN(dx, SMIN(dy, dz)));
    }

    bool is_shell_position(const Vecd &position) const
    {
        return distance_to_boundary(position) <= shell_thickness_ + TinyReal;
    }

    Real shell_thickness_;
    Real *Vol_;
    Vecd *positions_;
    Vecd *target_;
    Vecd *core_;
};

class ExtendCurlNuBFromCoreNeighborsInwardBiased : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ExtendCurlNuBFromCoreNeighborsInwardBiased(BaseInnerRelation &inner_relation,
                                                        const std::string &target_name,
                                                        const std::string &core_name,
                                                        Real shell_thickness,
                                                        const Vec3d &box_center,
                                                        const Vec3d &box_halfsize)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          shell_thickness_(shell_thickness),
          box_center_(box_center),
          box_halfsize_(box_halfsize),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          target_(particles_->getVariableDataByName<Vecd>(target_name)),
          core_(particles_->getVariableDataByName<Vecd>(core_name))
    {
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        if (!is_shell_position(positions_[index_i]))
        {
            target_[index_i] = core_[index_i];
            return;
        }

        Vecd inward_normal = compute_inward_normal(positions_[index_i]);
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vecd weighted_inward_value = ZeroData<Vecd>::value;
        Real inward_weight_sum = 0.0;
        Vecd weighted_core_value = ZeroData<Vecd>::value;
        Real core_weight_sum = 0.0;

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (is_shell_position(positions_[index_j]))
            {
                continue;
            }
            Real base_weight = inner_neighborhood.W_ij_[n] * Vol_[index_j];
            if (base_weight <= TinyReal)
            {
                continue;
            }

            weighted_core_value += base_weight * core_[index_j];
            core_weight_sum += base_weight;

            Vecd inward_displacement = positions_[index_j] - positions_[index_i];
            Real inward_alignment = inward_displacement.dot(inward_normal) /
                                     (inward_displacement.norm() + TinyReal);
            if (inward_alignment <= TinyReal)
            {
                continue;
            }

            Real inward_weight = base_weight * inward_alignment;
            weighted_inward_value += inward_weight * core_[index_j];
            inward_weight_sum += inward_weight;
        }

        if (inward_weight_sum > TinyReal)
        {
            target_[index_i] = weighted_inward_value / inward_weight_sum;
        }
        else if (core_weight_sum > TinyReal)
        {
            target_[index_i] = weighted_core_value / core_weight_sum;
        }
    }

  protected:
    Real distance_to_boundary(const Vecd &position) const
    {
        Vec3d offset = position - box_center_;
        Vec3d margin = box_halfsize_ - offset.cwiseAbs();
        return SMAX(static_cast<Real>(0.0), margin.minCoeff());
    }

    bool is_shell_position(const Vecd &position) const
    {
        return distance_to_boundary(position) <= shell_thickness_ + TinyReal;
    }

    Vecd compute_inward_normal(const Vecd &position) const
    {
        Vec3d offset = position - box_center_;
        Vec3d margin = box_halfsize_ - offset.cwiseAbs();
        int axis = 0;
        if (margin[1] < margin[axis]) axis = 1;
        if constexpr (Vecd::RowsAtCompileTime > 2)
        {
            if (margin[2] < margin[axis]) axis = 2;
        }
        Vecd normal = ZeroData<Vecd>::value;
        normal[axis] = (offset[axis] >= 0.0) ? -1.0 : 1.0;
        return normal;
    }

    Real shell_thickness_;
    Vec3d box_center_;
    Vec3d box_halfsize_;
    Real *Vol_;
    Vecd *positions_;
    Vecd *target_;
    Vecd *core_;
};

class ApplyDirichletBoundaryClosureFromExactA : public LocalDynamics
{
  public:
    explicit ApplyDirichletBoundaryClosureFromExactA(SPHBody &sph_body,
                                                     const std::string &target_name,
                                                     const std::string &core_name,
                                                     const std::string &vector_potential_name,
                                                     const SinusoidalAzField &model,
                                                     Real shell_thickness,
                                                     const Vec3d &box_center,
                                                     const Vec3d &box_halfsize,
                                                     Real magnetic_reluctivity)
        : LocalDynamics(sph_body),
          model_(model),
          shell_thickness_(shell_thickness),
          box_center_(box_center),
          box_halfsize_(box_halfsize),
          magnetic_reluctivity_(magnetic_reluctivity),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          target_(particles_->getVariableDataByName<Vecd>(target_name)),
          core_(particles_->getVariableDataByName<Vecd>(core_name)),
          vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        const Vec3d position = positions_[index_i];
        if (!is_shell_position(position))
        {
            target_[index_i] = core_[index_i];
            return;
        }

        Vec3d corrected = core_[index_i];
        const Vec3d a_i = vector_potential_[index_i];
        const Vec3d lower = box_center_ - box_halfsize_;
        const Vec3d upper = box_center_ + box_halfsize_;

        for (int axis = 0; axis != 3; ++axis)
        {
            Real dist_lower = position[axis] - lower[axis];
            if (dist_lower <= shell_thickness_ + TinyReal)
            {
                Vec3d boundary_point = position;
                boundary_point[axis] = lower[axis];
                Vec3d a_boundary = model_.evaluateA(boundary_point);
                Real denom = dist_lower * dist_lower + static_cast<Real>(0.01) * dp_0 * dp_0;
                corrected += static_cast<Real>(2.0) * magnetic_reluctivity_ * (a_i - a_boundary) / denom;
            }

            Real dist_upper = upper[axis] - position[axis];
            if (dist_upper <= shell_thickness_ + TinyReal)
            {
                Vec3d boundary_point = position;
                boundary_point[axis] = upper[axis];
                Vec3d a_boundary = model_.evaluateA(boundary_point);
                Real denom = dist_upper * dist_upper + static_cast<Real>(0.01) * dp_0 * dp_0;
                corrected += static_cast<Real>(2.0) * magnetic_reluctivity_ * (a_i - a_boundary) / denom;
            }
        }

        target_[index_i] = corrected;
    }

  protected:
    Real distance_to_boundary(const Vecd &position) const
    {
        Vec3d offset = position - box_center_;
        Vec3d margin = box_halfsize_ - offset.cwiseAbs();
        return SMAX(static_cast<Real>(0.0), margin.minCoeff());
    }

    bool is_shell_position(const Vecd &position) const
    {
        return distance_to_boundary(position) <= shell_thickness_ + TinyReal;
    }

    SinusoidalAzField model_;
    Real shell_thickness_;
    Vec3d box_center_;
    Vec3d box_halfsize_;
    Real magnetic_reluctivity_;
    Vecd *positions_;
    Vecd *target_;
    Vecd *core_;
    Vecd *vector_potential_;
};

class ApplyBoundaryFluxClosureFromExactA : public LocalDynamics
{
  public:
    explicit ApplyBoundaryFluxClosureFromExactA(SPHBody &sph_body,
                                                const std::string &target_name,
                                                const std::string &core_name,
                                                const std::string &vector_potential_name,
                                                const SinusoidalAzField &model,
                                                Real shell_thickness,
                                                const Vec3d &box_center,
                                                const Vec3d &box_halfsize,
                                                Real magnetic_reluctivity,
                                                Real smoothing_length)
        : LocalDynamics(sph_body),
          model_(model),
          shell_thickness_(shell_thickness),
          box_center_(box_center),
          box_halfsize_(box_halfsize),
          magnetic_reluctivity_(magnetic_reluctivity),
          smoothing_length_(smoothing_length),
          positions_(particles_->getVariableDataByName<Vecd>("Position")),
          target_(particles_->getVariableDataByName<Vecd>(target_name)),
          core_(particles_->getVariableDataByName<Vecd>(core_name)),
          vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        const Vec3d position = positions_[index_i];
        if (!is_shell_position(position))
        {
            target_[index_i] = core_[index_i];
            return;
        }

        Vec3d corrected = core_[index_i];
        const Vec3d a_i = vector_potential_[index_i];
        const Vec3d lower = box_center_ - box_halfsize_;
        const Vec3d upper = box_center_ + box_halfsize_;
        const Real h = smoothing_length_;
        const Real eps = static_cast<Real>(0.01) * h * h;

        for (int axis = 0; axis != 3; ++axis)
        {
            Real dist_lower = position[axis] - lower[axis];
            if (dist_lower <= shell_thickness_ + TinyReal)
            {
                Vec3d boundary_point = position;
                boundary_point[axis] = lower[axis];
                Vec3d a_boundary = model_.evaluateA(boundary_point);
                Real taper = SMAX(static_cast<Real>(0.0), static_cast<Real>(1.0) - dist_lower / (shell_thickness_ + TinyReal));
                Real flux_prefactor = static_cast<Real>(2.0) * taper * dist_lower /
                                      (h * (dist_lower * dist_lower + eps));
                corrected += magnetic_reluctivity_ * flux_prefactor * (a_i - a_boundary);
            }

            Real dist_upper = upper[axis] - position[axis];
            if (dist_upper <= shell_thickness_ + TinyReal)
            {
                Vec3d boundary_point = position;
                boundary_point[axis] = upper[axis];
                Vec3d a_boundary = model_.evaluateA(boundary_point);
                Real taper = SMAX(static_cast<Real>(0.0), static_cast<Real>(1.0) - dist_upper / (shell_thickness_ + TinyReal));
                Real flux_prefactor = static_cast<Real>(2.0) * taper * dist_upper /
                                      (h * (dist_upper * dist_upper + eps));
                corrected += magnetic_reluctivity_ * flux_prefactor * (a_i - a_boundary);
            }
        }

        target_[index_i] = corrected;
    }

  protected:
    Real distance_to_boundary(const Vecd &position) const
    {
        Vec3d offset = position - box_center_;
        Vec3d margin = box_halfsize_ - offset.cwiseAbs();
        return SMAX(static_cast<Real>(0.0), margin.minCoeff());
    }

    bool is_shell_position(const Vecd &position) const
    {
        return distance_to_boundary(position) <= shell_thickness_ + TinyReal;
    }

    SinusoidalAzField model_;
    Real shell_thickness_;
    Vec3d box_center_;
    Vec3d box_halfsize_;
    Real magnetic_reluctivity_;
    Real smoothing_length_;
    Vecd *positions_;
    Vecd *target_;
    Vecd *core_;
    Vecd *vector_potential_;
};

Vec3d angular_to_vec(const AngularVecd &value)
{
    return Vec3d(value[0], value[1], value[2]);
}

Real compute_relative_error(const Vec3d &discrete, const Vec3d &exact)
{
    return (discrete - exact).norm() / (exact.norm() + TinyReal);
}

Real distance_to_body_boundary(const Vec3d &position)
{
    Real dx = SMIN(position[0], body_length - position[0]);
    Real dy = SMIN(position[1], body_height - position[1]);
    Real dz = SMIN(position[2], body_width - position[2]);
    return SMAX(static_cast<Real>(0.0), SMIN(dx, SMIN(dy, dz)));
}

bool is_boundary_shell_particle(const Vec3d &position)
{
    return distance_to_body_boundary(position) <= boundary_shell_thickness;
}

struct ComponentFieldConfig
{
    std::string component_name;
    std::string vector_potential_name;
    std::string coupled_vector_potential_name;
    std::string grad_phi_name;
    std::string source_name;
    std::string curl_name;
    std::string curl_nu_b_name;
    std::string change_rate_name;
    std::string relative_residual_name;
    const SinusoidalAzField *model;
    Real omega_coupling_sign;
};

struct ComponentSummary
{
    std::string component_name;
    std::string solver_mode;
    std::string operator_mode;
    std::string solve_phase;
    size_t total_particles = 0;
    int iterations = 0;
    Real initial_guess_scale = 0.0;
    Real mean_a_error = 0.0;
    Real max_a_error = 0.0;
    Real mean_b_error = 0.0;
    Real max_b_error = 0.0;
    Real mean_curl_nu_b_error = 0.0;
    Real max_curl_nu_b_error = 0.0;
    Real mean_relative_curl_nu_b_error = 0.0;
    Real max_relative_curl_nu_b_error = 0.0;
    Real mean_residual_norm = 0.0;
    Real max_residual_norm = 0.0;
    Real mean_relative_residual = 0.0;
    Real max_relative_residual = 0.0;
    Real mean_change_rate_norm = 0.0;
    Real max_change_rate_norm = 0.0;
    size_t max_curl_nu_b_error_particle_id = 0;
    Vec3d max_curl_nu_b_error_position = ZeroData<Vec3d>::value;
    Real max_curl_nu_b_error_distance_to_boundary = 0.0;
    int max_curl_nu_b_error_in_shell = 0;
    Real max_curl_nu_b_error_exact_norm = 0.0;
};

struct ComponentIterationMetrics
{
    Real mean_b_error = 0.0;
    Real max_b_error = 0.0;
    Real mean_curl_nu_b_error = 0.0;
    Real max_curl_nu_b_error = 0.0;
    Real mean_relative_curl_nu_b_error = 0.0;
    Real max_relative_curl_nu_b_error = 0.0;
    Real mean_residual_norm = 0.0;
    Real max_residual_norm = 0.0;
    Real mean_relative_residual = 0.0;
    Real max_relative_residual = 0.0;
    Real mean_change_rate_norm = 0.0;
    Real max_change_rate_norm = 0.0;
};

ComponentIterationMetrics evaluate_iteration_metrics(BaseParticles &particles,
                                                     const ComponentFieldConfig &config)
{
    ComponentIterationMetrics metrics;
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vector_potential = particles.getVariableDataByName<Vecd>(config.vector_potential_name);
    Vecd *coupled_vector_potential = particles.getVariableDataByName<Vecd>(config.coupled_vector_potential_name);
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>(config.grad_phi_name);
    Vecd *source = particles.getVariableDataByName<Vecd>(config.source_name);
    AngularVecd *curl_field = particles.getVariableDataByName<AngularVecd>(config.curl_name);
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>(config.curl_nu_b_name);
    Vecd *change_rate = particles.getVariableDataByName<Vecd>(config.change_rate_name);
    Real *relative_residual = particles.getVariableDataByName<Real>(config.relative_residual_name);
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_b_error = 0.0;
    Real sum_curl_nu_b_error = 0.0;
    Real sum_relative_curl_nu_b_error = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d b_exact = config.model->evaluateB(positions[i]);
        Vec3d curl_nu_b_exact = config.model->evaluateCurlNuB(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_field[i]);
        Vec3d residual = source[i] - curl_nu_b[i] - sigma[i] * grad_phi[i] +
                         config.omega_coupling_sign * sigma[i] * omega * coupled_vector_potential[i];
        Real b_error = (b_discrete - b_exact).norm();
        Real curl_nu_b_error = (curl_nu_b[i] - curl_nu_b_exact).norm();
        Real relative_curl_nu_b_error = compute_relative_error(curl_nu_b[i], curl_nu_b_exact);
        Real residual_norm = residual.norm();
        Real relative_residual_norm = relative_residual[i];
        Real change_rate_norm = change_rate[i].norm();

        sum_b_error += b_error;
        sum_curl_nu_b_error += curl_nu_b_error;
        sum_relative_curl_nu_b_error += relative_curl_nu_b_error;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual_norm;
        sum_change_rate_norm += change_rate_norm;
        metrics.max_b_error = SMAX(metrics.max_b_error, b_error);
        metrics.max_curl_nu_b_error = SMAX(metrics.max_curl_nu_b_error, curl_nu_b_error);
        metrics.max_relative_curl_nu_b_error = SMAX(metrics.max_relative_curl_nu_b_error, relative_curl_nu_b_error);
        metrics.max_residual_norm = SMAX(metrics.max_residual_norm, residual_norm);
        metrics.max_relative_residual = SMAX(metrics.max_relative_residual, relative_residual_norm);
        metrics.max_change_rate_norm = SMAX(metrics.max_change_rate_norm, change_rate_norm);
    }

    Real count = static_cast<Real>(total_real_particles) + TinyReal;
    metrics.mean_b_error = sum_b_error / count;
    metrics.mean_curl_nu_b_error = sum_curl_nu_b_error / count;
    metrics.mean_relative_curl_nu_b_error = sum_relative_curl_nu_b_error / count;
    metrics.mean_residual_norm = sum_residual_norm / count;
    metrics.mean_relative_residual = sum_relative_residual / count;
    metrics.mean_change_rate_norm = sum_change_rate_norm / count;
    return metrics;
}

ComponentSummary evaluate_component_summary(BaseParticles &particles,
                                           const ComponentFieldConfig &config)
{
    ComponentSummary summary;
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vector_potential = particles.getVariableDataByName<Vecd>(config.vector_potential_name);
    Vecd *coupled_vector_potential = particles.getVariableDataByName<Vecd>(config.coupled_vector_potential_name);
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>(config.grad_phi_name);
    Vecd *source = particles.getVariableDataByName<Vecd>(config.source_name);
    AngularVecd *curl_field = particles.getVariableDataByName<AngularVecd>(config.curl_name);
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>(config.curl_nu_b_name);
    Vecd *change_rate = particles.getVariableDataByName<Vecd>(config.change_rate_name);
    Real *relative_residual = particles.getVariableDataByName<Real>(config.relative_residual_name);
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");

    size_t total_real_particles = particles.TotalRealParticles();
    Real sum_a_error = 0.0;
    Real sum_b_error = 0.0;
    Real sum_curl_nu_b_error = 0.0;
    Real sum_relative_curl_nu_b_error = 0.0;
    Real sum_residual_norm = 0.0;
    Real sum_relative_residual = 0.0;
    Real sum_change_rate_norm = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = config.model->evaluateA(positions[i]);
        Vec3d b_exact = config.model->evaluateB(positions[i]);
        Vec3d curl_nu_b_exact = config.model->evaluateCurlNuB(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_field[i]);
        Vec3d residual = source[i] - curl_nu_b[i] - sigma[i] * grad_phi[i] +
                         config.omega_coupling_sign * sigma[i] * omega * coupled_vector_potential[i];
        Real a_error = (vector_potential[i] - a_exact).norm();
        Real b_error = (b_discrete - b_exact).norm();
        Real curl_nu_b_error = (curl_nu_b[i] - curl_nu_b_exact).norm();
        Real relative_curl_nu_b_error = compute_relative_error(curl_nu_b[i], curl_nu_b_exact);
        Real residual_norm = residual.norm();
        Real relative_residual_norm = relative_residual[i];
        Real change_rate_norm = change_rate[i].norm();

        summary.total_particles++;
        sum_a_error += a_error;
        sum_b_error += b_error;
        sum_curl_nu_b_error += curl_nu_b_error;
        sum_relative_curl_nu_b_error += relative_curl_nu_b_error;
        sum_residual_norm += residual_norm;
        sum_relative_residual += relative_residual_norm;
        sum_change_rate_norm += change_rate_norm;
        summary.max_a_error = SMAX(summary.max_a_error, a_error);
        summary.max_b_error = SMAX(summary.max_b_error, b_error);
        if (curl_nu_b_error > summary.max_curl_nu_b_error)
        {
            summary.max_curl_nu_b_error = curl_nu_b_error;
            summary.max_curl_nu_b_error_particle_id = i;
            summary.max_curl_nu_b_error_position = positions[i];
            summary.max_curl_nu_b_error_distance_to_boundary = distance_to_body_boundary(positions[i]);
            summary.max_curl_nu_b_error_in_shell = is_boundary_shell_particle(positions[i]) ? 1 : 0;
            summary.max_curl_nu_b_error_exact_norm = curl_nu_b_exact.norm();
        }
        summary.max_relative_curl_nu_b_error = SMAX(summary.max_relative_curl_nu_b_error, relative_curl_nu_b_error);
        summary.max_residual_norm = SMAX(summary.max_residual_norm, residual_norm);
        summary.max_relative_residual = SMAX(summary.max_relative_residual, relative_residual_norm);
        summary.max_change_rate_norm = SMAX(summary.max_change_rate_norm, change_rate_norm);
    }

    Real count = static_cast<Real>(summary.total_particles) + TinyReal;
    summary.mean_a_error = sum_a_error / count;
    summary.mean_b_error = sum_b_error / count;
    summary.mean_curl_nu_b_error = sum_curl_nu_b_error / count;
    summary.mean_relative_curl_nu_b_error = sum_relative_curl_nu_b_error / count;
    summary.mean_residual_norm = sum_residual_norm / count;
    summary.mean_relative_residual = sum_relative_residual / count;
    summary.mean_change_rate_norm = sum_change_rate_norm / count;
    return summary;
}

void write_particle_header(const std::string &file_path)
{
    std::ofstream file(file_path, std::ios::out | std::ios::trunc);
    file << std::setprecision(12);
    file << "component,particle_id,x,y,z,distance_to_boundary,is_shell,"
         << "a_x,a_y,a_z,a_exact_x,a_exact_y,a_exact_z,a_error,"
         << "curl_a_x,curl_a_y,curl_a_z,curl_a_exact_x,curl_a_exact_y,curl_a_exact_z,b_error,"
         << "curl_nu_b_x,curl_nu_b_y,curl_nu_b_z,curl_nu_b_exact_x,curl_nu_b_exact_y,curl_nu_b_exact_z,curl_nu_b_exact_norm,curl_nu_b_error,curl_nu_b_relative_error,"
         << "residual_x,residual_y,residual_z,residual_norm,relative_residual,"
         << "change_rate_x,change_rate_y,change_rate_z,change_rate_norm\n";
}

void write_particle_diagnostics(const std::string &file_path,
                                BaseParticles &particles,
                                const ComponentFieldConfig &config)
{
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    file << std::setprecision(12);
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vector_potential = particles.getVariableDataByName<Vecd>(config.vector_potential_name);
    Vecd *coupled_vector_potential = particles.getVariableDataByName<Vecd>(config.coupled_vector_potential_name);
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>(config.grad_phi_name);
    Vecd *source = particles.getVariableDataByName<Vecd>(config.source_name);
    AngularVecd *curl_field = particles.getVariableDataByName<AngularVecd>(config.curl_name);
    Vecd *curl_nu_b = particles.getVariableDataByName<Vecd>(config.curl_nu_b_name);
    Vecd *change_rate = particles.getVariableDataByName<Vecd>(config.change_rate_name);
    Real *relative_residual = particles.getVariableDataByName<Real>(config.relative_residual_name);
    Real *sigma = particles.getVariableDataByName<Real>("ElectricalConductivity");
    size_t total_real_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d a_exact = config.model->evaluateA(positions[i]);
        Vec3d b_exact = config.model->evaluateB(positions[i]);
        Vec3d curl_nu_b_exact = config.model->evaluateCurlNuB(positions[i]);
        Vec3d b_discrete = angular_to_vec(curl_field[i]);
        Vec3d residual = source[i] - curl_nu_b[i] - sigma[i] * grad_phi[i] +
                         config.omega_coupling_sign * sigma[i] * omega * coupled_vector_potential[i];
        Real boundary_distance = distance_to_body_boundary(positions[i]);
        Real curl_nu_b_relative_error = compute_relative_error(curl_nu_b[i], curl_nu_b_exact);
        file << config.component_name << ","
             << i << ","
             << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << ","
             << boundary_distance << ","
             << (is_boundary_shell_particle(positions[i]) ? 1 : 0) << ","
             << vector_potential[i][0] << "," << vector_potential[i][1] << "," << vector_potential[i][2] << ","
             << a_exact[0] << "," << a_exact[1] << "," << a_exact[2] << ","
             << (vector_potential[i] - a_exact).norm() << ","
             << b_discrete[0] << "," << b_discrete[1] << "," << b_discrete[2] << ","
             << b_exact[0] << "," << b_exact[1] << "," << b_exact[2] << ","
             << (b_discrete - b_exact).norm() << ","
             << curl_nu_b[i][0] << "," << curl_nu_b[i][1] << "," << curl_nu_b[i][2] << ","
             << curl_nu_b_exact[0] << "," << curl_nu_b_exact[1] << "," << curl_nu_b_exact[2] << ","
             << curl_nu_b_exact.norm() << ","
             << (curl_nu_b[i] - curl_nu_b_exact).norm() << ","
             << curl_nu_b_relative_error << ","
             << residual[0] << "," << residual[1] << "," << residual[2] << ","
             << residual.norm() << ","
             << relative_residual[i] << ","
             << change_rate[i][0] << "," << change_rate[i][1] << "," << change_rate[i][2] << ","
             << change_rate[i].norm() << "\n";
    }
}

} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_FREQ_CURLCURL_VERIFY_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::string solver_mode = to_lower_copy(solver_mode_raw);
    if (solver_mode != "single" && solver_mode != "coupled" && solver_mode != "both" &&
        solver_mode != "global" && solver_mode != "all")
    {
        std::cout << "[em-freq-curlcurl-verify-warning] unsupported solver mode '"
                  << solver_mode_raw << "', fallback to 'coupled'." << std::endl;
        solver_mode = "coupled";
    }
    const bool run_single = solver_mode_contains(solver_mode, "single");
    const bool run_coupled = solver_mode_contains(solver_mode, "coupled");
    const bool run_global = solver_mode_contains(solver_mode, "global");
    std::string operator_mode = to_lower_copy(operator_mode_raw);
    if (operator_mode != "direct" && operator_mode != "component_hessian" &&
        operator_mode != "hybrid_core_hessian_shell_direct" &&
        operator_mode != "hybrid_core_hessian_shell_coreavg" &&
        operator_mode != "hybrid_core_hessian_shell_inwardavg" &&
        operator_mode != "hybrid_core_hessian_shell_dirichletbc" &&
        operator_mode != "hybrid_core_hessian_shell_boundaryflux")
    {
        std::cout << "[em-freq-curlcurl-verify-warning] unsupported operator mode '"
                  << operator_mode_raw << "', fallback to 'direct'." << std::endl;
        operator_mode = "direct";
    }
    const bool use_component_hessian_operator = (operator_mode == "component_hessian");
    const bool use_hybrid_operator = (operator_mode == "hybrid_core_hessian_shell_direct");
    const bool use_hybrid_shell_coreavg_operator = (operator_mode == "hybrid_core_hessian_shell_coreavg");
    const bool use_hybrid_shell_inwardavg_operator = (operator_mode == "hybrid_core_hessian_shell_inwardavg");
    const bool use_hybrid_shell_dirichletbc_operator = (operator_mode == "hybrid_core_hessian_shell_dirichletbc");
    const bool use_hybrid_shell_boundaryflux_operator = (operator_mode == "hybrid_core_hessian_shell_boundaryflux");
    std::string global_operator_mode = operator_mode;
    if (run_global && global_operator_mode != "direct")
    {
        std::cout << "[em-freq-curlcurl-verify-warning] global solver currently uses the direct magnetic operator action; "
                  << "fallback from operator_mode='" << operator_mode << "' to 'direct'." << std::endl;
        global_operator_mode = "direct";
    }
    std::string global_alpha_formula = to_lower_copy(global_alpha_formula_raw);
    if (global_alpha_formula != "preconditioned_sd" && global_alpha_formula != "residual_min")
    {
        std::cout << "[em-freq-curlcurl-verify-warning] unsupported global alpha formula '"
                  << global_alpha_formula_raw << "', fallback to 'preconditioned_sd'." << std::endl;
        global_alpha_formula = "preconditioned_sd";
    }
    const int effective_history_interval = SMAX(1, history_interval);

    std::cout << std::scientific << std::setprecision(6)
              << "[em-freq-curlcurl-verify-config] dp=" << dp_0
              << ", sigma=" << conductivity
              << ", frequency_hz=" << frequency_hz
              << ", omega=" << omega
              << ", kx=" << kx
              << ", ky=" << ky
              << ", amplitude_real=" << amplitude_real
              << ", amplitude_imag=" << amplitude_imag
              << ", iterations=" << solve_iterations
              << ", dt=" << dt_pseudo
              << ", reference_dt=" << reference_dt_pseudo
              << ", solver_mode=" << solver_mode
              << ", operator_mode=" << operator_mode
              << ", global_alpha_formula=" << global_alpha_formula
              << ", global_alpha_max=" << global_alpha_max
              << std::endl;

    SolidBody conductor_body(sph_system, makeShared<ConductorShape>("Conductor"));
    conductor_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    conductor_body.defineMaterial<Solid>();
    conductor_body.defineBodyLevelSetShape();
    conductor_body.generateParticles<BaseParticles, Lattice>();
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_real_component_names[0], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_real_component_names[1], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_real_component_names[2], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_imag_component_names[0], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_imag_component_names[1], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>(a_imag_component_names[2], Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<AngularVecd>("GlobalSearchCurlReal", ZeroData<AngularVecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<AngularVecd>("GlobalSearchCurlImag", ZeroData<AngularVecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("GlobalSearchCurlNuBReal", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("GlobalSearchCurlNuBImag", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("GlobalOperatorActionReal", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("GlobalOperatorActionImag", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("CurlNuBRealHybridCore", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("CurlNuBImagHybridCore", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("AEquationResidualVectorReal", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Vecd>("AEquationResidualVectorImag", ZeroData<Vecd>::value);
    conductor_body.getBaseParticles().registerStateVariableData<Real>("AEquationRelativeResidualReal", Real(0));
    conductor_body.getBaseParticles().registerStateVariableData<Real>("AEquationRelativeResidualImag", Real(0));

    BodyRegionByParticle boundary_region(
        conductor_body,
        makeShared<OuterBoundaryShellShape>("ConductorBoundaryShell"));
    InnerRelation conductor_inner(conductor_body);
    using MainExecutionPolicy = execution::ParallelPolicy;
    Inner<> conductor_inner_ck(conductor_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_conductor_cell_linked_list_ck(conductor_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_conductor_inner_ck(conductor_inner_ck);

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
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        conductor_linear_correction_ck(DynamicsArgs(conductor_inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>>>
        conductor_displacement_matrix_gradient_ck(conductor_inner_ck);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>>>
        conductor_hessian_correction_ck(DynamicsArgs(conductor_inner_ck, 0.0));
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_a_real_components(conductor_body, "VectorPotentialReal", a_real_component_names);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_real_gradient_x_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_real_hessian_x_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_real_gradient_y_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_real_hessian_y_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_real_gradient_z_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_real_hessian_z_ck(DynamicsArgs(conductor_inner_ck, a_real_component_names[2]));
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_curl_nu_b_real_from_components(
            conductor_body, a_real_component_names, curl_nu_b_scaling, "CurlNuBReal");
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_curl_nu_b_real_hybrid_core_from_components(
            conductor_body, a_real_component_names, curl_nu_b_scaling, "CurlNuBRealHybridCore");
    SimpleDynamics<SelectCurlNuBHybridByBoundaryShell>
        select_curl_nu_b_real_hybrid(conductor_body, "CurlNuBReal", "CurlNuBRealHybridCore",
                                     body_center, body_halfsize, boundary_shell_thickness);
    InteractionDynamics<ExtendCurlNuBFromCoreNeighbors>
        extend_curl_nu_b_real_from_core_neighbors(conductor_inner, "CurlNuBReal", "CurlNuBRealHybridCore",
                                                  boundary_shell_thickness);
    InteractionDynamics<ExtendCurlNuBFromCoreNeighborsInwardBiased>
        extend_curl_nu_b_real_from_core_neighbors_inward(conductor_inner, "CurlNuBReal", "CurlNuBRealHybridCore",
                                                         boundary_shell_thickness, body_center, body_halfsize);
    SimpleDynamics<ApplyDirichletBoundaryClosureFromExactA>
        apply_dirichlet_boundary_closure_real(conductor_body, "CurlNuBReal", "CurlNuBRealHybridCore",
                                              "VectorPotentialReal", a_real_model, boundary_shell_thickness,
                                              body_center, body_halfsize, magnetic_reluctivity);
    SimpleDynamics<ApplyBoundaryFluxClosureFromExactA>
        apply_boundary_flux_closure_real(conductor_body, "CurlNuBReal", "CurlNuBRealHybridCore",
                                         "VectorPotentialReal", a_real_model, boundary_shell_thickness,
                                         body_center, body_halfsize, magnetic_reluctivity,
                                         conductor_body.getSPHAdaptation().ReferenceSmoothingLength());
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_a_imag_components(conductor_body, "VectorPotentialImag", a_imag_component_names);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_imag_gradient_x_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_imag_hessian_x_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_imag_gradient_y_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_imag_hessian_y_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        a_imag_gradient_z_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        a_imag_hessian_z_ck(DynamicsArgs(conductor_inner_ck, a_imag_component_names[2]));
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_curl_nu_b_imag_from_components(
            conductor_body, a_imag_component_names, curl_nu_b_scaling, "CurlNuBImag");
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_curl_nu_b_imag_hybrid_core_from_components(
            conductor_body, a_imag_component_names, curl_nu_b_scaling, "CurlNuBImagHybridCore");
    SimpleDynamics<SelectCurlNuBHybridByBoundaryShell>
        select_curl_nu_b_imag_hybrid(conductor_body, "CurlNuBImag", "CurlNuBImagHybridCore",
                                     body_center, body_halfsize, boundary_shell_thickness);
    InteractionDynamics<ExtendCurlNuBFromCoreNeighbors>
        extend_curl_nu_b_imag_from_core_neighbors(conductor_inner, "CurlNuBImag", "CurlNuBImagHybridCore",
                                                  boundary_shell_thickness);
    InteractionDynamics<ExtendCurlNuBFromCoreNeighborsInwardBiased>
        extend_curl_nu_b_imag_from_core_neighbors_inward(conductor_inner, "CurlNuBImag", "CurlNuBImagHybridCore",
                                                         boundary_shell_thickness, body_center, body_halfsize);
    SimpleDynamics<ApplyDirichletBoundaryClosureFromExactA>
        apply_dirichlet_boundary_closure_imag(conductor_body, "CurlNuBImag", "CurlNuBImagHybridCore",
                                              "VectorPotentialImag", a_imag_model, boundary_shell_thickness,
                                              body_center, body_halfsize, magnetic_reluctivity);
    SimpleDynamics<ApplyBoundaryFluxClosureFromExactA>
        apply_boundary_flux_closure_imag(conductor_body, "CurlNuBImag", "CurlNuBImagHybridCore",
                                         "VectorPotentialImag", a_imag_model, boundary_shell_thickness,
                                         body_center, body_halfsize, magnetic_reluctivity,
                                         conductor_body.getSPHAdaptation().ReferenceSmoothingLength());

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
    InteractionDynamics<electromagnetics::VectorPotentialFrequencyCoupledPreconditionerInner>
        global_preconditioned_search(conductor_inner, omega,
                                     "AEquationResidualVectorReal", "AEquationResidualVectorImag",
                                     "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
                                     sigma_relaxation_scaling, sigma_relaxation_floor,
                                     magnetic_diagonal_scaling, max_change_rate);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        global_search_curl_real_inner(conductor_inner,
                                      "VectorPotentialChangeRateReal", "GlobalSearchCurlReal", curl_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        global_search_curl_imag_inner(conductor_inner,
                                      "VectorPotentialChangeRateImag", "GlobalSearchCurlImag", curl_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        global_search_curl_nu_b_real_inner(conductor_inner,
                                           "GlobalSearchCurlReal", "GlobalSearchCurlNuBReal", curl_nu_b_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        global_search_curl_nu_b_imag_inner(conductor_inner,
                                           "GlobalSearchCurlImag", "GlobalSearchCurlNuBImag", curl_nu_b_scaling);
    SimpleDynamics<electromagnetics::FrequencyVectorPotentialLinearOperatorComplex>
        global_search_operator(conductor_body, omega,
                               "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
                               "GlobalSearchCurlNuBReal", "GlobalSearchCurlNuBImag",
                               "GlobalOperatorActionReal", "GlobalOperatorActionImag");
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
    write_states.addToWrite<Vecd>(conductor_body, "AEquationResidualVectorReal");
    write_states.addToWrite<Vecd>(conductor_body, "AEquationResidualVectorImag");
    write_states.addToWrite<Real>(conductor_body, "AEquationRelativeResidualReal");
    write_states.addToWrite<Real>(conductor_body, "AEquationRelativeResidualImag");

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    initialize_frequency_em.exec();
    assign_manufactured_auxiliary_fields.exec();

    const ComponentFieldConfig real_config{
        "A_real", "VectorPotentialReal", "VectorPotentialImag",
        "ElectricPotentialGradientReal", "SourceCurrentDensityReal",
        "VectorPotentialCurlReal", "CurlNuBReal", "VectorPotentialChangeRateReal",
        "AEquationRelativeResidualReal",
        &a_real_model, 1.0};
    const ComponentFieldConfig imag_config{
        "A_imag", "VectorPotentialImag", "VectorPotentialReal",
        "ElectricPotentialGradientImag", "SourceCurrentDensityImag",
        "VectorPotentialCurlImag", "CurlNuBImag", "VectorPotentialChangeRateImag",
        "AEquationRelativeResidualImag",
        &a_imag_model, -1.0};
    auto refresh_equation_residual_diagnostics = [&]()
    {
        diagnose_a_real.exec();
        diagnose_a_imag.exec();
    };

    const std::string history_path =
        io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_history.csv";
    std::ofstream history_file(history_path, std::ios::out | std::ios::trunc);
    history_file << std::setprecision(12);
    history_file << "solver_mode,operator_mode,solve_phase,iteration,component,"
                 << "mean_b_error,max_b_error,"
                 << "mean_curl_nu_b_error,max_curl_nu_b_error,"
                 << "mean_relative_curl_nu_b_error,max_relative_curl_nu_b_error,"
                 << "mean_residual_norm,max_residual_norm,"
                 << "mean_relative_residual,max_relative_residual,"
                 << "mean_change_rate_norm,max_change_rate_norm\n";
    history_file.close();

    auto write_history_row = [&](const std::string &mode,
                                 const std::string &operator_label,
                                 const std::string &phase,
                                 int iteration,
                                 const ComponentFieldConfig &config)
    {
        ComponentIterationMetrics metrics =
            evaluate_iteration_metrics(conductor_body.getBaseParticles(), config);
        std::ofstream append_file(history_path, std::ios::out | std::ios::app);
        append_file << std::setprecision(12);
        append_file << mode << ","
                    << operator_label << ","
                    << phase << ","
                    << iteration << ","
                    << config.component_name << ","
                    << metrics.mean_b_error << ","
                    << metrics.max_b_error << ","
                    << metrics.mean_curl_nu_b_error << ","
                    << metrics.max_curl_nu_b_error << ","
                    << metrics.mean_relative_curl_nu_b_error << ","
                    << metrics.max_relative_curl_nu_b_error << ","
                    << metrics.mean_residual_norm << ","
                    << metrics.max_residual_norm << ","
                    << metrics.mean_relative_residual << ","
                    << metrics.max_relative_residual << ","
                    << metrics.mean_change_rate_norm << ","
                    << metrics.max_change_rate_norm << "\n";
    };

    auto set_summary_meta = [&](ComponentSummary &summary,
                                const ComponentFieldConfig &config,
                                const std::string &mode,
                                const std::string &operator_label,
                                const std::string &phase,
                                Real initial_scale)
    {
        summary.component_name = config.component_name;
        summary.solver_mode = mode;
        summary.operator_mode = operator_label;
        summary.solve_phase = phase;
        summary.iterations = solve_iterations;
        summary.initial_guess_scale = initial_scale;
    };

    bool component_hessian_ck_prepared = false;
    auto ensure_component_hessian_ck_prepared = [&]()
    {
        if (component_hessian_ck_prepared)
        {
            return;
        }
        update_conductor_cell_linked_list_ck.exec();
        update_conductor_inner_ck.exec();
        conductor_linear_correction_ck.exec();
        conductor_displacement_matrix_gradient_ck.exec();
        conductor_hessian_correction_ck.exec();
        component_hessian_ck_prepared = true;
    };
    auto refresh_real_magnetic_operator = [&]()
    {
        curl_a_real_inner.exec();
        if (use_component_hessian_operator)
        {
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_from_components.exec();
        }
        else if (use_hybrid_operator)
        {
            curl_nu_b_real_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_hybrid_core_from_components.exec();
            select_curl_nu_b_real_hybrid.exec();
        }
        else if (use_hybrid_shell_coreavg_operator)
        {
            curl_nu_b_real_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_hybrid_core_from_components.exec();
            extend_curl_nu_b_real_from_core_neighbors.exec();
        }
        else if (use_hybrid_shell_inwardavg_operator)
        {
            curl_nu_b_real_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_hybrid_core_from_components.exec();
            extend_curl_nu_b_real_from_core_neighbors_inward.exec();
        }
        else if (use_hybrid_shell_dirichletbc_operator)
        {
            curl_nu_b_real_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_hybrid_core_from_components.exec();
            apply_dirichlet_boundary_closure_real.exec();
        }
        else if (use_hybrid_shell_boundaryflux_operator)
        {
            curl_nu_b_real_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_real_components.exec();
            a_real_gradient_x_ck.exec();
            a_real_hessian_x_ck.exec();
            a_real_gradient_y_ck.exec();
            a_real_hessian_y_ck.exec();
            a_real_gradient_z_ck.exec();
            a_real_hessian_z_ck.exec();
            reconstruct_curl_nu_b_real_hybrid_core_from_components.exec();
            apply_boundary_flux_closure_real.exec();
        }
        else
        {
            curl_nu_b_real_inner.exec();
        }
    };
    auto refresh_imag_magnetic_operator = [&]()
    {
        curl_a_imag_inner.exec();
        if (use_component_hessian_operator)
        {
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_from_components.exec();
        }
        else if (use_hybrid_operator)
        {
            curl_nu_b_imag_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_hybrid_core_from_components.exec();
            select_curl_nu_b_imag_hybrid.exec();
        }
        else if (use_hybrid_shell_coreavg_operator)
        {
            curl_nu_b_imag_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_hybrid_core_from_components.exec();
            extend_curl_nu_b_imag_from_core_neighbors.exec();
        }
        else if (use_hybrid_shell_inwardavg_operator)
        {
            curl_nu_b_imag_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_hybrid_core_from_components.exec();
            extend_curl_nu_b_imag_from_core_neighbors_inward.exec();
        }
        else if (use_hybrid_shell_dirichletbc_operator)
        {
            curl_nu_b_imag_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_hybrid_core_from_components.exec();
            apply_dirichlet_boundary_closure_imag.exec();
        }
        else if (use_hybrid_shell_boundaryflux_operator)
        {
            curl_nu_b_imag_inner.exec();
            ensure_component_hessian_ck_prepared();
            copy_a_imag_components.exec();
            a_imag_gradient_x_ck.exec();
            a_imag_hessian_x_ck.exec();
            a_imag_gradient_y_ck.exec();
            a_imag_hessian_y_ck.exec();
            a_imag_gradient_z_ck.exec();
            a_imag_hessian_z_ck.exec();
            reconstruct_curl_nu_b_imag_hybrid_core_from_components.exec();
            apply_boundary_flux_closure_imag.exec();
        }
        else
        {
            curl_nu_b_imag_inner.exec();
        }
    };
    auto refresh_global_search_operator_direct = [&]()
    {
        global_search_curl_real_inner.exec();
        global_search_curl_imag_inner.exec();
        global_search_curl_nu_b_real_inner.exec();
        global_search_curl_nu_b_imag_inner.exec();
        global_search_operator.exec();
    };
    auto compute_global_alpha = [&]()
    {
        BaseParticles &particles = conductor_body.getBaseParticles();
        Vecd *residual_real = particles.getVariableDataByName<Vecd>("AEquationResidualVectorReal");
        Vecd *residual_imag = particles.getVariableDataByName<Vecd>("AEquationResidualVectorImag");
        Vecd *search_real = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
        Vecd *search_imag = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
        Vecd *operator_action_real = particles.getVariableDataByName<Vecd>("GlobalOperatorActionReal");
        Vecd *operator_action_imag = particles.getVariableDataByName<Vecd>("GlobalOperatorActionImag");
        size_t total_real_particles = particles.TotalRealParticles();
        Real r_dot_z = 0.0;
        Real z_dot_az = 0.0;
        Real r_dot_az = 0.0;
        Real az_dot_az = 0.0;
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            r_dot_z += residual_real[i].dot(search_real[i]) +
                       residual_imag[i].dot(search_imag[i]);
            z_dot_az += search_real[i].dot(operator_action_real[i]) +
                        search_imag[i].dot(operator_action_imag[i]);
            r_dot_az += residual_real[i].dot(operator_action_real[i]) +
                        residual_imag[i].dot(operator_action_imag[i]);
            az_dot_az += operator_action_real[i].squaredNorm() +
                         operator_action_imag[i].squaredNorm();
        }

        Real alpha = 0.0;
        if (global_alpha_formula == "preconditioned_sd")
        {
            if (std::isfinite(z_dot_az) && fabs(z_dot_az) > TinyReal)
            {
                alpha = r_dot_z / z_dot_az;
            }
            else if (std::isfinite(az_dot_az) && az_dot_az > TinyReal)
            {
                alpha = r_dot_az / az_dot_az;
            }
        }
        else
        {
            if (std::isfinite(az_dot_az) && az_dot_az > TinyReal)
            {
                alpha = r_dot_az / az_dot_az;
            }
        }

        if (!std::isfinite(alpha))
        {
            alpha = 0.0;
        }
        Real alpha_cap = SMAX(static_cast<Real>(0.0), global_alpha_max);
        if (alpha_cap > TinyReal)
        {
            alpha = SMAX(-alpha_cap, SMIN(alpha_cap, alpha));
        }
        return alpha;
    };
    auto apply_global_search_update = [&](Real alpha)
    {
        BaseParticles &particles = conductor_body.getBaseParticles();
        Vecd *a_real = particles.getVariableDataByName<Vecd>("VectorPotentialReal");
        Vecd *a_imag = particles.getVariableDataByName<Vecd>("VectorPotentialImag");
        Vecd *search_real = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
        Vecd *search_imag = particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
        Real pseudo_step_scale = relaxation_scaling *
            electromagnetics::ComputeNormalizedPseudoTimeStep(dt_pseudo, reference_dt_pseudo);
        size_t total_real_particles = particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            a_real[i] += pseudo_step_scale * alpha * search_real[i];
            a_imag[i] += pseudo_step_scale * alpha * search_imag[i];
            if (!std::isfinite(a_real[i].squaredNorm()))
            {
                a_real[i] = ZeroData<Vecd>::value;
            }
            if (!std::isfinite(a_imag[i].squaredNorm()))
            {
                a_imag[i] = ZeroData<Vecd>::value;
            }
        }
    };

    std::vector<ComponentSummary> summaries;

    if (run_single)
    {
        set_initial_guess_real.exec();
        if (use_boundary_constraint)
        {
            constrain_a_real_boundary.exec();
        }
        refresh_real_magnetic_operator();
        refresh_equation_residual_diagnostics();
        write_history_row("single", operator_mode, "phase1_real", 0, real_config);
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            refresh_real_magnetic_operator();
            solve_a_real_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_real_boundary.exec();
            }
            if ((iter + 1) % effective_history_interval == 0 || (iter + 1 == solve_iterations))
            {
                refresh_real_magnetic_operator();
                refresh_equation_residual_diagnostics();
                write_history_row("single", operator_mode, "phase1_real", iter + 1, real_config);
            }
        }
        refresh_real_magnetic_operator();
        refresh_equation_residual_diagnostics();
        ComponentSummary real_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), real_config);
        set_summary_meta(real_summary, real_config, "single", operator_mode, "phase1_real", initial_guess_scale_real);
        summaries.push_back(real_summary);
        if (write_particles)
        {
            const std::string particles_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_particles_phase1_real.csv";
            write_particle_header(particles_path);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), real_config);
        }

        assign_manufactured_auxiliary_fields.exec();
        set_initial_guess_imag.exec();
        if (use_boundary_constraint)
        {
            constrain_a_imag_boundary.exec();
        }
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();
        write_history_row("single", operator_mode, "phase2_imag", 0, imag_config);
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            refresh_imag_magnetic_operator();
            solve_a_imag_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_imag_boundary.exec();
            }
            if ((iter + 1) % effective_history_interval == 0 || (iter + 1 == solve_iterations))
            {
                refresh_imag_magnetic_operator();
                refresh_equation_residual_diagnostics();
                write_history_row("single", operator_mode, "phase2_imag", iter + 1, imag_config);
            }
        }
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();
        ComponentSummary imag_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), imag_config);
        set_summary_meta(imag_summary, imag_config, "single", operator_mode, "phase2_imag", initial_guess_scale_imag);
        summaries.push_back(imag_summary);
        if (write_particles)
        {
            const std::string particles_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_particles_phase2_imag.csv";
            write_particle_header(particles_path);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), imag_config);
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
        refresh_real_magnetic_operator();
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();
        write_history_row("coupled", operator_mode, "coupled", 0, real_config);
        write_history_row("coupled", operator_mode, "coupled", 0, imag_config);
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            refresh_real_magnetic_operator();
            refresh_imag_magnetic_operator();
            solve_a_coupled_inner.exec(dt_pseudo);
            if (use_boundary_constraint)
            {
                constrain_a_real_boundary.exec();
                constrain_a_imag_boundary.exec();
            }
            if ((iter + 1) % effective_history_interval == 0 || (iter + 1 == solve_iterations))
            {
                refresh_real_magnetic_operator();
                refresh_imag_magnetic_operator();
                refresh_equation_residual_diagnostics();
                write_history_row("coupled", operator_mode, "coupled", iter + 1, real_config);
                write_history_row("coupled", operator_mode, "coupled", iter + 1, imag_config);
            }
        }
        refresh_real_magnetic_operator();
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();

        ComponentSummary real_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), real_config);
        set_summary_meta(real_summary, real_config, "coupled", operator_mode, "coupled", initial_guess_scale_real);
        summaries.push_back(real_summary);

        ComponentSummary imag_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), imag_config);
        set_summary_meta(imag_summary, imag_config, "coupled", operator_mode, "coupled", initial_guess_scale_imag);
        summaries.push_back(imag_summary);

        if (write_particles)
        {
            const std::string particles_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_particles_coupled.csv";
            write_particle_header(particles_path);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), real_config);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), imag_config);
        }
    }

    if (run_global)
    {
        assign_manufactured_auxiliary_fields.exec();
        set_initial_guess_real.exec();
        set_initial_guess_imag.exec();
        if (use_boundary_constraint)
        {
            constrain_a_real_boundary.exec();
            constrain_a_imag_boundary.exec();
        }
        refresh_real_magnetic_operator();
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();
        write_history_row("global", global_operator_mode, "global", 0, real_config);
        write_history_row("global", global_operator_mode, "global", 0, imag_config);
        for (int iter = 0; iter != solve_iterations; ++iter)
        {
            refresh_real_magnetic_operator();
            refresh_imag_magnetic_operator();
            refresh_equation_residual_diagnostics();
            global_preconditioned_search.exec();
            refresh_global_search_operator_direct();
            Real alpha = compute_global_alpha();
            apply_global_search_update(alpha);
            if (use_boundary_constraint)
            {
                constrain_a_real_boundary.exec();
                constrain_a_imag_boundary.exec();
            }
            if ((iter + 1) % effective_history_interval == 0 || (iter + 1 == solve_iterations))
            {
                refresh_real_magnetic_operator();
                refresh_imag_magnetic_operator();
                refresh_equation_residual_diagnostics();
                write_history_row("global", global_operator_mode, "global", iter + 1, real_config);
                write_history_row("global", global_operator_mode, "global", iter + 1, imag_config);
                std::cout << "[em-freq-curlcurl-global] iter=" << (iter + 1)
                          << ", alpha=" << alpha << std::endl;
            }
        }
        refresh_real_magnetic_operator();
        refresh_imag_magnetic_operator();
        refresh_equation_residual_diagnostics();

        ComponentSummary real_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), real_config);
        set_summary_meta(real_summary, real_config, "global", global_operator_mode, "global", initial_guess_scale_real);
        summaries.push_back(real_summary);

        ComponentSummary imag_summary =
            evaluate_component_summary(conductor_body.getBaseParticles(), imag_config);
        set_summary_meta(imag_summary, imag_config, "global", global_operator_mode, "global", initial_guess_scale_imag);
        summaries.push_back(imag_summary);

        if (write_particles)
        {
            const std::string particles_path =
                io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_particles_global.csv";
            write_particle_header(particles_path);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), real_config);
            write_particle_diagnostics(particles_path, conductor_body.getBaseParticles(), imag_config);
        }
    }

    if (write_vtp)
    {
        write_states.writeToFile(0.0);
    }

    const std::string summary_path =
        io_environment.OutputFolder() + "/em_frequency_a_solver_nonzero_curlcurl_summary.csv";
    std::ofstream summary_file(summary_path, std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "component,solver_mode,operator_mode,solve_phase,total_particles,iterations,initial_guess_scale,"
                 << "mean_a_error,max_a_error,"
                 << "mean_b_error,max_b_error,"
                 << "mean_curl_nu_b_error,max_curl_nu_b_error,"
                 << "mean_relative_curl_nu_b_error,max_relative_curl_nu_b_error,"
                 << "mean_residual_norm,max_residual_norm,"
                 << "mean_relative_residual,max_relative_residual,"
                 << "mean_change_rate_norm,max_change_rate_norm,"
                 << "max_curl_nu_b_error_particle_id,max_curl_nu_b_error_x,max_curl_nu_b_error_y,max_curl_nu_b_error_z,"
                 << "max_curl_nu_b_error_distance_to_boundary,max_curl_nu_b_error_in_shell,max_curl_nu_b_error_exact_norm\n";
    for (const ComponentSummary &summary : summaries)
    {
        summary_file << summary.component_name << ","
                     << summary.solver_mode << ","
                     << summary.operator_mode << ","
                     << summary.solve_phase << ","
                     << summary.total_particles << ","
                     << summary.iterations << ","
                     << summary.initial_guess_scale << ","
                     << summary.mean_a_error << ","
                     << summary.max_a_error << ","
                     << summary.mean_b_error << ","
                     << summary.max_b_error << ","
                     << summary.mean_curl_nu_b_error << ","
                     << summary.max_curl_nu_b_error << ","
                     << summary.mean_relative_curl_nu_b_error << ","
                     << summary.max_relative_curl_nu_b_error << ","
                     << summary.mean_residual_norm << ","
                     << summary.max_residual_norm << ","
                     << summary.mean_relative_residual << ","
                     << summary.max_relative_residual << ","
                     << summary.mean_change_rate_norm << ","
                     << summary.max_change_rate_norm << ","
                     << summary.max_curl_nu_b_error_particle_id << ","
                     << summary.max_curl_nu_b_error_position[0] << ","
                     << summary.max_curl_nu_b_error_position[1] << ","
                     << summary.max_curl_nu_b_error_position[2] << ","
                     << summary.max_curl_nu_b_error_distance_to_boundary << ","
                     << summary.max_curl_nu_b_error_in_shell << ","
                     << summary.max_curl_nu_b_error_exact_norm << "\n";
    }
    summary_file.flush();

    for (const ComponentSummary &summary : summaries)
    {
        std::cout << std::scientific << std::setprecision(6)
                  << "[em-freq-curlcurl-verify-summary] mode=" << summary.solver_mode
                  << ", operator_mode=" << summary.operator_mode
                  << ", phase=" << summary.solve_phase
                  << ", component=" << summary.component_name
                  << ", mean|A-A_exact|=" << summary.mean_a_error
                  << ", mean|curlA-B_exact|=" << summary.mean_b_error
                  << ", mean|curlNuB-curlNuB_exact|=" << summary.mean_curl_nu_b_error
                  << ", mean(relative curlNuB error)=" << summary.mean_relative_curl_nu_b_error
                  << ", mean|residual|=" << summary.mean_residual_norm
                  << ", mean(relative residual)=" << summary.mean_relative_residual
                  << ", mean|dA/dtau|=" << summary.mean_change_rate_norm
                  << std::endl;
    }

    std::cout << "[em-freq-curlcurl-verify] summary file: " << summary_path << std::endl;
    std::cout << "[em-freq-curlcurl-verify] history file: " << history_path << std::endl;
    return 0;
}
