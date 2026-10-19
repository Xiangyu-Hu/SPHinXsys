// SPDX-License-Identifier: Apache-2.0
#ifndef HJC_SOLID_H
#define HJC_SOLID_H

#include "elastic_solid.h"

namespace SPH
{
/** HJC parameters in consistent density, stress and time units.
 * mu_lock is permanent volumetric compaction, not total compression at P_lock.
 * See the HJC impact example for equations, integration conventions and references.
 */
struct HJCParameters
{
    Real shear_modulus, A, B, C, N, compressive_strength, tensile_strength;
    Real reference_strain_rate, minimum_fracture_strain, maximum_strength;
    Real crush_pressure, crush_strain, lock_pressure, lock_strain;
    Real D1, D2, K1, K2, K3;
};

struct HJCState
{
    Mat3d stress = Mat3d::Zero();
    Real damage = 0, plastic_strain = 0, plastic_volume = 0, maximum_compression = 0;
};

/** Holmquist-Johnson-Cook concrete, with a corotational incremental update.
 * Use with solid_dynamics::HJCIntegration1stHalf and Integration2ndHalf.
 * State-free elastic stress interfaces deliberately reject accidental use.
 */
class HJCSolid : public ElasticSolid
{
  public:
    HJCSolid(Real rho0, const HJCParameters &parameters);
    void initializeLocalParameters(BaseParticles *particles) override;
    const HJCParameters &Parameters() const { return parameters_; }
    Real LockCompression() const { return lock_compression_; }
    Real FractureStrain(Real pressure) const;
    Real YieldStrength(Real pressure, Real damage, Real strain_rate) const;
    /** Uncapped EOS pressure; updates irreversible crushing history. */
    Real Pressure(Real compression, HJCState &state) const;
    /** Symmetric logarithmic increment in an objective frame; compression = 1/J-1. */
    void Integrate(const Mat3d &strain_increment, Real compression, Real dt, HJCState &state) const;
    /** Update particle history from successive deformation gradients. */
    Matd UpdateStress(const Matd &deformation, size_t index_i, Real dt);
    Real AcousticModulus(Real compression) const;

    Matd StressPK1(Matd &deformation, size_t index_i) override;
    Matd StressPK2(Matd &deformation, size_t index_i) override;
    Matd StressCauchy(Matd &strain, size_t index_i) override;
    Real VolumetricKirchhoff(Real J) override;
    std::string getRelevantStressMeasureName() override { return "Cauchy"; }

  private:
    HJCParameters parameters_;
    Real lock_compression_, transition_slope_;
    Mat3d *stress_;
    Matd *previous_deformation_;
    Real *damage_, *plastic_strain_, *plastic_volume_, *maximum_compression_;
    Real *pressure_, *equivalent_stress_, *acoustic_modulus_;
    Real DensePressure(Real compression) const;
};
} // namespace SPH
#endif
