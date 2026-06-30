#include "sphinxsys.h"
using namespace SPH;

// geometry
const Real L = 1;
const Vec3d domain_lower_bound(3 * L, 2 * L, 0 * L);
const Vec3d domain_upper_bound(10 * L, 5.6 * L, 5 * L);
const auto plate_tip_pos = Vec3d(7 * L, 3.3 * L, 0 * L);
const Real plate_width = L;
const Real plate_height = L;
const Real plate_y = L / 4.0;
const Real plate_z = L / 2.0;
const Real plate_thickness = 0.01 * L;

// material
// fluid
const Real rho0_f = 1.0;
const Real U_f = 5.48;
const Real Re = 200;
const Real mu_f = rho0_f * U_f * L / Re;
const Real c_f = 10.0 * U_f;

// solid
const Real rho0_s = 1200;          /**< Reference density.*/
const Real youngs_modulus = 3.5e9; /**< Youngs modulus.*/
const Real poisson = 0.32;         /**< Poisson ratio.*/

// Cycle
const Real flow_init_time = 0.5;
const Real fsi_start_time = 2.0;
const Real end_time = fsi_start_time + 9.0;

// Shell structure
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
    }
};

struct ShellObject
{
    SolidBody body_;
    std::unique_ptr<ShellAlgorithms> algs_;

    ShellObject(SPHSystem &sph_system, const std::string &name, const StdVec<Vec3d> &positions, const StdVec<Vec3d> &normals, Real dp, Real thickness)
        : body_(sph_system, makeShared<DefaultShape>(name))
    {
        body_.defineAdaptation<SPHAdaptation>(1.15, sph_system.GlobalResolution() / dp);
        body_.defineMatterMaterial<LinearElasticSolid>(rho0_s, youngs_modulus, poisson);
        body_.generateParticles<SurfaceParticles, Shell>(positions, normals, dp, thickness);
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

inline void relax_solid(RealBody &body, BaseInnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    ReloadParticleIO write_particle_reload_files(body);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    body.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
        }
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    write_particle_reload_files.writeToFile(0);
}

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