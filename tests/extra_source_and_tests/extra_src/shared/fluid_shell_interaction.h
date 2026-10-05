/**
 * @file fluid_shell_interaction.h
 * @brief This file defines the immersed fluid-shell interaction.
 */

#ifndef FLUID_SHELL_INTERACTION_H
#define FLUID_SHELL_INTERACTION_H

#include "base_local_dynamics.h"

namespace SPH
{
/**
 * @class ShellFluidMixtureMass
 * @brief Reset the mass of shell particles considering the fluid volume.
 */
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
          dp_(shell_body.getSPHBodyResolutionRef()),
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
} // namespace SPH

#endif // FLUID_SHELL_INTERACTION_H
