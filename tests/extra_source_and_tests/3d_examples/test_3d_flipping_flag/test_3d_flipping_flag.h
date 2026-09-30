/**
 * @file 	test_3d_flipping_flag.h
 * @brief 	3D flipping flag example
 * @details This is the one test case for the 3D shell immersed in fluid.
 * Reference: Tianrun Gao, Lin Fu, https://doi.org/10.1016/j.cma.2024.117179.
 * A three-dimensional meshless fluid–shell interaction framework based on smoothed particle hydrodynamics coupled with semi-meshless thin shell
 * HUANG W-X, SUNG HJ. doi:10.1017/S0022112010000248
 * Three-dimensional simulation of a flapping flag in a uniform flow. Journal of Fluid Mechanics.
 * @author 	Weiyi Kong (Virtonomy GmbH) and Xiangyu Hu
 */
#include "fluid_shell_interaction.h"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Parameters
//----------------------------------------------------------------------
// geometry
const Real L = 1;
const Real inclination_angle = Pi / 10.0;
const Real plate_thickness = 0.01 * L;
// The domain is shrank to save computational cost.
const Real L_x_before = L;
const Real L_x_after = 4 * L;
const Real L_y = 2 * L;
const Real L_z = 3 * L;
const auto beam_observation_location = StdVec<Vec3d>{Vec3d(L * cos(inclination_angle), -0.5 * L, L *sin(inclination_angle))};

// material
const Real rho0_f = 1.0;                 /**< Reference density of fluid.*/
const Real U_f = 1.0;                    /**< Reference velocity of fluid. */
const Real Re = 200;                     /**< Reynolds number. */
const Real mu_f = rho0_f * U_f * L / Re; /**< Dynamic viscosity of fluid. */
const Real c_f = 25.0;                   /**< Reference sound speed of fluid. */

// solid
const Real rho0_s = 100 * rho0_f;                                                                                            /**< Reference density.*/
const Real poisson = 0.4;                                                                                                    /**< Poisson ratio.*/
const Real youngs_modulus = 0.0001 * 12.0 * (1 - poisson * poisson) * rho0_f * U_f * U_f * std::pow(L / plate_thickness, 3); /**< Youngs modulus.*/

// Cycle
const Real flow_init_time = 3.0;
const Real fsi_start_time = flow_init_time + 0;
const Real end_time = fsi_start_time + 7.0;

// Reference solution, Huang et al.
const Real A_dL_ref = 0.78; // amplitude/L
const Real St_ref = 0.260;  // Strouhal number
//----------------------------------------------------------------------
//	Boundary conditions
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
    Real u_ref_ = U_f;
    Real t_ref_ = flow_init_time;

    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &) {}

    Vecd operator()(Vecd &, Vecd &, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        target_velocity[0] = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	Shell particle generator
//----------------------------------------------------------------------
namespace SPH
{
class Shell;
template <>
class ParticleGenerator<SurfaceParticles, Shell> : public ParticleGenerator<SurfaceParticles>
{
    const StdVec<Vec3d> &positions_;
    const StdVec<Vec3d> &normals_;
    Real dp_;
    Real thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               const StdVec<Vec3d> &positions,
                               const StdVec<Vec3d> &normals,
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
            addPositionAndVolumetricMeasure(positions_[i], dp_ * dp_);
            addSurfaceProperties(normals_[i], thickness_);
        }
    }
};
} // namespace SPH

//----------------------------------------------------------------------
//	Shell algorithms
//----------------------------------------------------------------------
inline Real get_physical_viscosity()
{
    return 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * plate_thickness * plate_thickness;
}

struct ShellAlgorithms
{
    InnerRelation inner_;
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration_;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half_;
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half_;
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_computing_time_step_size_;
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> shell_position_damping_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> shell_rotation_damping_;

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
        body.updateCellLinkedList();
        inner_.updateConfiguration();
        corrected_configuration_.exec();
    }
};

struct shell_inputs
{
    std::string name;
    StdVec<Vec3d> positions;
    StdVec<Vec3d> normals;
    Real dp;
    Real thickness;
};

struct ShellObject
{
    SolidBody body_;
    std::unique_ptr<ShellAlgorithms> algs_;

    ShellObject(SPHSystem &sph_system, const shell_inputs &inputs)
        : body_(sph_system, makeShared<DefaultShape>(inputs.name))
    {
        body_.defineAdaptation<SPHAdaptation>(1.15, sph_system.GlobalResolution() / inputs.dp);
        body_.defineMatterMaterial<LinearElasticSolid>(rho0_s, youngs_modulus, poisson);
        body_.generateParticles<SurfaceParticles, Shell>(inputs.positions, inputs.normals, inputs.dp, inputs.thickness);
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
          pressure_force_from_fluid_(contact_relation_)
    {
        SimpleDynamics<ShellFluidMixtureMass> reset_shell_mass(shell_body, rho0_f);
        reset_shell_mass.exec();
    }
};
//----------------------------------------------------------------------
//	Fluid algorithms
//----------------------------------------------------------------------
namespace SPH
{
class FreeSlipWall;

namespace fluid_dynamics
{
inline Vecd get_slip_vel(const Vec3d &v_i, const Vec3d &v_k, const Vec3d &n_k)
{
    Vecd v_n = v_i.dot(n_k) * n_k;
    Vecd v_t = v_i - v_n;
    Vecd v_n_wall = 2 * v_k.dot(n_k) * n_k - v_n;
    return v_t + v_n_wall;
}

template <class RiemannSolverType>
class Integration2ndHalf<Contact<FreeSlipWall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  private:
    RiemannSolverType riemann_solver_;

  public:
    explicit Integration2ndHalf(BaseContactRelation &wall_contact_relation)
        : BaseIntegrationWithWall(wall_contact_relation),
          riemann_solver_(this->fluid_, this->fluid_) {};
    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real density_change_rate = 0.0;
        Vecd p_dissipation = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const Vecd *vel_ave_k = this->wall_vel_ave_[k];
            const Vecd *n_k = this->wall_n_[k];
            const Real *wall_Vol_k = this->wall_Vol_[k];
            Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                size_t index_j = wall_neighborhood.j_[n];
                Vecd &e_ij = wall_neighborhood.e_ij_[n];
                Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];

                Vecd face_to_fluid_n = SGN(e_ij.dot(n_k[index_j])) * n_k[index_j];

                Vecd vel_j_in_wall = get_slip_vel(this->vel_[index_i], vel_ave_k[index_j], n_k[index_j]);
                density_change_rate += (this->vel_[index_i] - vel_j_in_wall).dot(e_ij) * dW_ijV_j;
                Real u_jump = 2.0 * (this->vel_[index_i] - vel_ave_k[index_j]).dot(face_to_fluid_n);
                p_dissipation += this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * face_to_fluid_n;
            }
        }
        this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
        this->force_[index_i] += p_dissipation * this->Vol_[index_i];
    }
};

template <typename ViscosityType, class KernelCorrectionType>
class ViscousForce<Contact<FreeSlipWall>, ViscosityType, KernelCorrectionType>
    : public BaseViscousForceWithWall
{
  private:
    ViscosityType mu_;
    KernelCorrectionType kernel_correction_;

  public:
    explicit ViscousForce(BaseContactRelation &wall_contact_relation)
        : BaseViscousForceWithWall(wall_contact_relation),
          mu_(particles_), kernel_correction_(particles_) {}
    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const Vecd *vel_ave_k = wall_vel_ave_[k];
            const Real *wall_Vol_k = wall_Vol_[k];
            const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real r_ij = contact_neighborhood.r_ij_[n];
                const Vecd &e_ij = contact_neighborhood.e_ij_[n];

                Vec3d vel_wall = get_slip_vel(vel_[index_i], vel_ave_k[index_j], wall_n_[k][index_j]);
                Vecd vel_derivative = (vel_[index_i] - vel_wall) /
                                      (r_ij + 0.01 * smoothing_length_);
                force += 2.0 * e_ij.dot(kernel_correction_(index_i) * e_ij) * mu_(index_i, index_i) *
                         vel_derivative * contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            }
        }

        viscous_force_[index_i] += force * Vol_[index_i];
    }
};
} // namespace fluid_dynamics
} // namespace SPH

struct FlagObject : public ShellObject
{
    auto get_inputs(Real dp) const
    {
        Real cos_theta = cos(inclination_angle);
        Real sin_theta = sin(inclination_angle);
        StdVec<Vec3d> positions;
        StdVec<Vec3d> normals;
        {
            Real y = 0.5 * L - 0.5 * dp;
            while (y > -0.5 * L)
            {
                Real r = 0.5 * dp;
                while (r < L)
                {
                    Real x = r * cos_theta;
                    Real z = r * sin_theta;
                    positions.emplace_back(x, y, z);
                    normals.emplace_back(sin_theta, 0.0, -cos_theta);
                    r += dp;
                }
                y -= dp;
            }
        };
        return shell_inputs{"flag", positions, normals, dp, plate_thickness};
    }
    FlagObject(SPHSystem &sph_system, Real dp_s)
        : ShellObject(sph_system, get_inputs(dp_s)) {}
};

struct FlagPinCondition
{
    BodyPartByParticle part;
    SimpleDynamics<FixBodyPartConstraint> constraint;

    explicit FlagPinCondition(SPHBody &sph_body)
        : part(sph_body), constraint(part)
    {
        Real dp = sph_body.getSPHBodyResolutionRef();
        auto clamp_id = [&, dp2 = dp * dp]()
        {
            IndexVector ids;
            const auto *pos = sph_body.getBaseParticles().getVariableDataByName<Vec3d>("Position");
            for (size_t i = 0; i < sph_body.getBaseParticles().TotalRealParticles(); ++i)
            {
                Real r2 = pos[i].x() * pos[i].x() + pos[i].z() * pos[i].z();
                if (r2 < dp2)
                    ids.push_back(i);
            }
            return ids;
        }();
        part.body_part_particles_ = clamp_id;
    }

    void exec()
    {
        constraint.exec();
    }
};