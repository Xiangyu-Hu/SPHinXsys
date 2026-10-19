// SPDX-License-Identifier: Apache-2.0
#include "hjc_dynamics.h"
#include "adaptation.h"

namespace SPH
{
namespace solid_dynamics
{
HJCIntegration1stHalf::HJCIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation), hjc_solid_(DynamicCast<HJCSolid>(this, elastic_solid_)) {}

void HJCIntegration1stHalf::initialization(size_t i, Real dt)
{
    pos_[i] += vel_[i] * dt * 0.5;
    F_[i] += dF_dt_[i] * dt * 0.5;
    Matd sigma = hjc_solid_.UpdateStress(F_[i], i, dt);
    Real J = F_[i].determinant();
    rho_[i] = rho0_ / J;
    stress_PK1_B_[i] = J * sigma * F_[i].inverse().transpose() * B_[i].transpose();
}

HJCAcousticTimeStep::HJCAcousticTimeStep(SPHBody &body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMin<Real>>(body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      modulus_(particles_->getVariableDataByName<Real>("HJCAcousticModulus")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      smoothing_length_(getSPHAdaptation().ReferenceSmoothingLength()), acousticCFL_(acousticCFL) {}

Real HJCAcousticTimeStep::reduce(size_t i, Real dt)
{
    Real acceleration = (force_[i] + force_prior_[i]).norm() / mass_[i];
    return acousticCFL_ * std::min(std::sqrt(smoothing_length_ / (acceleration + TinyReal)),
        smoothing_length_ / (std::sqrt(modulus_[i] / rho_[i]) + vel_[i].norm()));
}
} // namespace solid_dynamics
} // namespace SPH
