#include "sphinxsys.h"
#include <unsupported/Eigen/Splines>
using namespace SPH;

constexpr Real cm_to_m = 1e-2;
constexpr Real g_to_kg = 1e-3;
constexpr Real poise = 0.1;
constexpr Real dyne = 1e-5;

// geometry
const Real height = 1.61 * cm_to_m;
const Real length = 8 * cm_to_m;
const Real shell_length = 0.7 * cm_to_m;
const Real shell_thickness = 0.0212 * cm_to_m;
const Real shell_pos_x = 2 * cm_to_m;

const Real dp_fluid = height / 40.0;
const Real dp_solid = dp_fluid;
const Real dp_shell = dp_fluid;

const Real buffer_length = dp_fluid * 3.0;
const Real wall_thickness = dp_fluid * 4.0;

// material
// fluid
const Real rho0_f = 100 * g_to_kg / std::pow(cm_to_m, 3);  /**< Density. */
const Real a = 5.0 / cm_to_m;                              // parameter for parabolic inflow profile
const Real shift = 1.1;                                    // parameter for parabolic inflow profile
const Real U_f = a * (1 + shift) * 0.25 * height * height; /**< Characteristic velocity. */
const Real U_max = 2 * U_f;                                // 2 is a safety factor for the maximum velocity
const Real c_f = 10.0 * U_max;                             /**< Speed of sound. */
const Real mu_f = 10 * poise;                              /**< Dynamics viscosity. */

// solid
const Real rho0_s = rho0_f;                                      /**< Reference density.*/
const Real poisson = 0.4;                                        /**< Poisson ratio.*/
const Real youngs_modulus = 5.6e7 * dyne / std::pow(cm_to_m, 2); /**< Youngs modulus.*/

// Cycle
const Real time_flow_init = 0;
const Real time_cycle = 1.0;
const size_t num_cycles = 2;

namespace SPH
{
class Shell;
template <>
class ParticleGenerator<SurfaceParticles, Shell> : public ParticleGenerator<SurfaceParticles>
{
    const StdVec<Vec2d> &positions_;
    const StdVec<Vec2d> &normals_;
    Real dp_;
    Real thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               const StdVec<Vec2d> &positions,
                               const StdVec<Vec2d> &normals,
                               Real dp,
                               Real thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          positions_(positions), normals_(normals), dp_(dp), thickness_(thickness)
    {
        if (positions_.size() != normals_.size())
        {
            std::cout << "Error: In ParticleGenerator<Shell>, positions size is not equal to normals size!" << std::endl;
            exit(1);
        }
    };
    void prepareGeometricData() override
    {
        const auto particle_number = positions_.size();
        // generate particles for the elastic gate
        for (size_t i = 0; i < particle_number; i++)
        {
            addPositionAndVolumetricMeasure(positions_[i], dp_);
            addSurfaceProperties(normals_[i], thickness_);
        }
    }
};
} // namespace SPH

class ShellFluidMixtureMass : public LocalDynamics
{
  private:
    Real rho_f0_; // assume the density of fluid is constant for now
    Real dp_;     // initial particle spacing
    Real *thickness_;
    Real *mass_;
    Real *Vol_;

  public:
    ShellFluidMixtureMass(SPHBody &shell_body, Real rho_f0)
        : LocalDynamics(shell_body),
          rho_f0_(rho_f0),
          dp_(shell_body.getSPHAdaptation().ReferenceSpacing()),
          thickness_(shell_body.getBaseParticles().getVariableDataByName<Real>("Thickness")),
          mass_(shell_body.getBaseParticles().getVariableDataByName<Real>("Mass")),
          Vol_(shell_body.getBaseParticles().getVariableDataByName<Real>("VolumetricMeasure"))
    {
    }

    void update(size_t index_i, Real)
    {
        Real dp_m_t = dp_ - thickness_[index_i];
        if (dp_m_t < 0)
            throw std::runtime_error("Error: In ShellFluidMixtureMass, dp - thickness < 0!");
        Real V_f = Vol_[index_i] * dp_m_t; // fluid volume
        Real mass_f = V_f * rho_f0_;
        mass_[index_i] += mass_f;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    explicit LeftInflowPressure(BoundaryConditionType &) {}

    Real operator()(Real p, Real)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    explicit RightInflowPressure(BoundaryConditionType &) {}

    Real operator()(Real p, Real)
    {
        return p;
    }
};

struct InflowVelocity
{
    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &) {}

    Vecd operator()(Vecd &pos, Vecd &, Real time)
    {
        Real y = pos.y();
        Real u = a * (sin(2.0 * Pi * time / time_cycle) + shift) * (0.25 * height * height - y * y);
        return u * Vec2d::UnitX();
    }
};

struct OutflowVelocity
{
    template <class BoundaryConditionType>
    explicit OutflowVelocity(BoundaryConditionType &) {}

    Vecd operator()(Vecd &, Vecd &vel, Real)
    {
        Vecd target_velocity = vel;
        target_velocity.y() = 0.0;
        return target_velocity;
    }
};

// Shell wrappers
inline Real get_physical_viscosity()
{
    return 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * shell_thickness * shell_thickness;
}

struct ShellAlgorithms
{
    InnerRelation inner_;
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration_;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half_;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half_;
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_computing_time_step_size_;
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> shell_position_damping_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> shell_rotation_damping_;

    explicit ShellAlgorithms(RealBody &body)
        : inner_(body),
          corrected_configuration_(inner_),
          stress_relaxation_first_half_(inner_, 3, true),
          stress_relaxation_second_half_(inner_),
          shell_computing_time_step_size_(body),
          update_normal_(body),
          shell_position_damping_(0.2, inner_, "Velocity", get_physical_viscosity()),
          shell_rotation_damping_(0.2, inner_, "AngularVelocity", get_physical_viscosity())
    {
    }
};

struct ShellObject
{
    SolidBody body_;
    std::unique_ptr<ShellAlgorithms> algs_;

    ShellObject(SPHSystem &sph_system, const std::string &name, const StdVec<Vec2d> &positions, const StdVec<Vec2d> &normals)
        : body_(sph_system, makeShared<DefaultShape>(name))
    {
        body_.defineAdaptation<SPHAdaptation>(1.15, dp_fluid / dp_shell);
        body_.defineMatterMaterial<NeoHookeanSolid>(rho0_s, youngs_modulus, poisson);
        body_.generateParticles<SurfaceParticles, Shell>(positions, normals, dp_shell, shell_thickness);
        SimpleDynamics<ShellFluidMixtureMass> reset_shell_mass(body_, rho0_f);
        reset_shell_mass.exec();
        algs_ = std::make_unique<ShellAlgorithms>(body_);
    }
};

template <class FluidIntegration2ndHalfType>
struct ShellFluidAlgorithms
{
    ContactRelationSFI2 contact_relation_;
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_;
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid_;
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<FluidIntegration2ndHalfType>> pressure_force_from_fluid_;

    ShellFluidAlgorithms(RealBody &shell_body, const RealBodyVector &fluid_bodies)
        : contact_relation_(shell_body, fluid_bodies),
          average_velocity_and_acceleration_(shell_body),
          viscous_force_from_fluid_(contact_relation_),
          pressure_force_from_fluid_(contact_relation_) {}
};

// Real reference data
// Colume 1: time, column 2: displacement
inline std::vector<std::vector<Real>> read_csv_data(const std::string &file_name, bool has_header = true)
{
    std::cout << "read_csv_data started" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_csv_data: file doesn't exist");

    std::vector<std::vector<Real>> data;

    std::string line;
    if (has_header)
        std::getline(my_file, line); // skip header
    while (std::getline(my_file, line))
    {
        std::vector<Real> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ','))
            row.emplace_back(std::stod(cell));

        data.push_back(row);
    }
    my_file.close();
    // convert to column-wise data
    std::vector<std::vector<Real>> column_data;
    if (!data.empty())
    {
        size_t num_columns = data[0].size();
        column_data.resize(num_columns);
        for (const auto &row : data)
        {
            for (size_t i = 0; i < num_columns; ++i)
            {
                column_data[i].push_back(row[i]);
            }
        }
    }
    std::cout << "read_csv_data finished" << std::endl;
    return column_data;
}

class LinearInterpolator
{
  private:
    std::vector<Real> x_;
    std::vector<Real> y_;

  public:
    LinearInterpolator(std::vector<Real> x,
                       std::vector<Real> y)
        : x_(std::move(x)), y_(std::move(y))
    {
        if (x_.size() != y_.size())
            throw std::runtime_error("Error: In LinearInterpolator, x and y must have the same size!");
    }

    Real operator()(Real xq) const
    {
        if (xq <= x_.front())
            return y_.front();

        if (xq >= x_.back())
            return y_.back();

        auto it = std::lower_bound(x_.begin(), x_.end(), xq);
        size_t i = std::distance(x_.begin(), it) - 1;

        Real t = (xq - x_[i]) / (x_[i + 1] - x_[i]);

        return y_[i] + t * (y_[i + 1] - y_[i]);
    }
};

struct return_data
{
    std::vector<Real> time_vec;
    std::vector<Real> disp_x_vec;
    std::vector<Real> disp_y_vec;
    std::vector<Real> disp_x_ref_vec;
    std::vector<Real> disp_y_ref_vec;
};

inline void write_data_to_csv(const std::string &file_name, const return_data &data)
{
    std::ofstream file(file_name);
    if (!file.is_open())
        throw std::runtime_error("write_to_csv: file doesn't exist");

    // Write header
    file << "time [s],disp_x [cm],disp_y [cm],disp_x_ref [cm],disp_y_ref [cm]\n";

    // Write data
    for (size_t i = 0; i < data.time_vec.size(); ++i)
    {
        file << data.time_vec[i] << ","
             << data.disp_x_vec[i] << ","
             << data.disp_y_vec[i] << ","
             << data.disp_x_ref_vec[i] << ","
             << data.disp_y_ref_vec[i] << "\n";
    }
    file.close();
}

inline auto get_average_error(const return_data &data, Real start_time, Real end_time)
{
    Real error_sum_x = 0.0;
    Real error_sum_y = 0.0;
    size_t n = data.time_vec.size();
    size_t count = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (data.time_vec[i] < start_time || data.time_vec[i] > end_time)
            continue;
        Real dx = data.disp_x_vec[i] - data.disp_x_ref_vec[i];
        Real dy = data.disp_y_vec[i] - data.disp_y_ref_vec[i];
        error_sum_x += std::abs(dx);
        error_sum_y += std::abs(dy);
        ++count;
    }
    Real max_disp_x_ref = *std::max_element(data.disp_x_ref_vec.begin(), data.disp_x_ref_vec.end());
    Real max_disp_y_ref = *std::max_element(data.disp_y_ref_vec.begin(), data.disp_y_ref_vec.end());
    Real error_x_ave = count > 0 ? error_sum_x / Real(count) / max_disp_x_ref : 0.0;
    Real error_y_ave = count > 0 ? error_sum_y / Real(count) / max_disp_y_ref : 0.0;
    return std::make_pair(error_x_ave, error_y_ave);
}