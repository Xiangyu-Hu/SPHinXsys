/**
 * @file 	case.h
 * @brief 	This is the case study the curvature of shell.
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real radius = 1;
const int particle_number_mid_surface = 100;
const Real resolution_ref = 2 * M_PI * radius / Real(particle_number_mid_surface);
const Real shell_thickness = resolution_ref;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-radius, -radius), Vec2d(radius, radius));
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** Particle generator and constraint boundary for shell baffle. */
// x=R*cos(theta), y=R*sin(theta)
class ShellParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real theta = 2 * M_PI / Real(particle_number_mid_surface) * (i + 0.5);
            Real x = radius * cos(theta);
            Real y = radius * sin(theta);
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_ref);
            Vec2d normal_direction = Vec2d(x / radius, y / radius);
            initializeSurfaceProperties(normal_direction, shell_thickness);
        }
    }
};

inline Real get_mean_curvature(const Matd &dn)
{
    return 0.5 * dn.trace();
}

inline Real get_Gaussian_curvature(Real H, const Matd &dn)
{
    Real sum = 0;
    for (int i = 0; i < Dimensions; i++)
        for (int j = 0; j < Dimensions; j++)
            sum += dn(i, j) * dn(i, j);
    return 0.5 * (4 * H * H - sum);
}

inline Real get_Gaussian_curvature(const Matd &dn)
{
    return get_Gaussian_curvature(dn.trace(), dn);
}

/**
 * @class ShellDynamicsInitialCondition
 * @brief  set initial condition for shell particles
 * This is a abstract class to be override for case specific initial conditions.
 */
class ShellInitialCurvature : public LocalDynamics, public thin_structure_dynamics::ShellDataInner
{
  public:
    explicit ShellInitialCurvature(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), thin_structure_dynamics::ShellDataInner(inner_relation),
          n0_(particles_->n0_), B_(particles_->B_), transformation_matrix_(particles_->transformation_matrix_)
    {
        particles_->registerVariable(dn_0_, "InitialNormalGradient");
        particles_->registerVariable(H_, "MeanCurvature");
        particles_->registerVariable(K_, "GaussianCurvature");
    };
    virtual ~ShellInitialCurvature(){};

    inline void update(size_t index_i, Real dt)
    {
        Matd dn_0_i = Matd::Zero();
        // transform initial local B_ to global B_
        const Matd B_global_i = transformation_matrix_[index_i].transpose() * B_[index_i] * transformation_matrix_[index_i];
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            const size_t index_j = inner_neighborhood.j_[n];
            const Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            dn_0_i -= (n0_[index_i] - n0_[index_j]) * gradW_ijV_j.transpose();
        }
        dn_0_[index_i] = dn_0_i * B_global_i;
        H_[index_i] = get_mean_curvature(dn_0_[index_i]);
        K_[index_i] = get_Gaussian_curvature(H_[index_i], dn_0_[index_i]);
    }

  protected:
    StdLargeVec<Vecd> &n0_;
    StdLargeVec<Matd> &B_;
    StdLargeVec<Matd> &transformation_matrix_;
    StdLargeVec<Matd> dn_0_;
    StdLargeVec<Real> H_;
    StdLargeVec<Real> K_;
};