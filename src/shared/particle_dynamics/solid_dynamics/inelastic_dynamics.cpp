#include "inelastic_dynamics.h"

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
DecomposedPlasticIntegration1stHalf::
    DecomposedPlasticIntegration1stHalf(BaseInnerRelation &inner_relation)
    : DecomposedIntegration1stHalf(inner_relation),
      plastic_solid_(DynamicCast<PlasticSolid>(this, elastic_solid_)),
      scaling_matrix_(particles_->registerStateVariable<Matd>("ScalingMatrix")),
      inverse_F_(particles_->registerStateVariable<Matd>("InverseDeformation")) {}
//=================================================================================================//
void DecomposedPlasticIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Real one_over_J = 1.0 / J;
    rho_[index_i] = rho0_ * one_over_J;

    Matd normalized_be = plastic_solid_.ElasticLeftCauchy(F_[index_i], index_i, dt);
    inverse_F_[index_i] = F_[index_i].inverse();
    Matd inverse_F_T = inverse_F_[index_i].transpose();
    scaling_matrix_[index_i] = normalized_be * inverse_F_T;
    Real isotropic_stress = plastic_solid_.ShearModulus() * normalized_be.trace() * OneOverDimensions;
    // Note that as we use small numerical damping here, the time step size (CFL number) may need to be decreased.
    stress_on_particle_[index_i] =
        inverse_F_T * (plastic_solid_.VolumetricKirchhoff(J) - isotropic_stress) +
        0.125 * plastic_solid_.NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i) * inverse_F_T;
}
//=================================================================================================//
} // namespace solid_dynamics
  //=====================================================================================================//
} // namespace SPH
