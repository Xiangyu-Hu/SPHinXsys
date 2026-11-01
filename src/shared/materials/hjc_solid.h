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

/** Device updates report errors for the host to check before advancing the step. */
enum class HJCIntegrationStatus : int
{
    success = 0,
    invalid_deformation = 1,
    invalid_increment = 2,
    singular_increment = 3,
    nonfinite_state = 4
};

/** Value-only numerical model shared by the CPU and computing-kernel interfaces. */
class HJCConstitutiveModel
{
  public:
    HJCConstitutiveModel(const HJCParameters &parameters, Real bulk_modulus,
                         Real shear_modulus, Real lock_compression, Real transition_slope)
        : parameters_(parameters), K0_(bulk_modulus), G0_(shear_modulus),
          lock_compression_(lock_compression), transition_slope_(transition_slope) {}

    Real DensePressure(Real compression) const;
    Real Pressure(Real compression, HJCState &state) const;
    Real FractureStrain(Real pressure) const;
    Real YieldStrength(Real pressure, Real damage, Real strain_rate) const;
    HJCIntegrationStatus Integrate(const Mat3d &increment, Real compression, Real dt, HJCState &state) const;
    HJCIntegrationStatus UpdateStress(const Matd &deformation, const Matd &previous_deformation,
                                      Real dt, HJCState &state) const;
    Real AcousticModulus(Real compression) const;

  private:
    HJCParameters parameters_;
    Real K0_, G0_, lock_compression_, transition_slope_;
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

    class ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, HJCSolid &encloser);
        Matd UpdateStress(const Matd &deformation, size_t index_i, Real dt);
        int Status(size_t index_i) const { return status_[index_i]; }

      private:
        HJCConstitutiveModel model_;
        Mat3d *stress_;
        Matd *previous_deformation_;
        Real *damage_, *plastic_strain_, *plastic_volume_, *maximum_compression_;
        Real *pressure_, *equivalent_stress_, *acoustic_modulus_;
        int *status_;
    };

  private:
    HJCParameters parameters_;
    Real lock_compression_, transition_slope_;
    Mat3d *stress_;
    Matd *previous_deformation_;
    Real *damage_, *plastic_strain_, *plastic_volume_, *maximum_compression_;
    Real *pressure_, *equivalent_stress_, *acoustic_modulus_;
    DiscreteVariable<Mat3d> *dv_stress_ = nullptr;
    DiscreteVariable<Matd> *dv_previous_deformation_ = nullptr;
    DiscreteVariable<Real> *dv_damage_ = nullptr, *dv_plastic_strain_ = nullptr;
    DiscreteVariable<Real> *dv_plastic_volume_ = nullptr, *dv_maximum_compression_ = nullptr;
    DiscreteVariable<Real> *dv_pressure_ = nullptr, *dv_equivalent_stress_ = nullptr, *dv_acoustic_modulus_ = nullptr;
    DiscreteVariable<int> *dv_status_ = nullptr;
    HJCConstitutiveModel ConstitutiveModel() const
    {
        return HJCConstitutiveModel(parameters_, K0_, G0_, lock_compression_, transition_slope_);
    }
    Real DensePressure(Real compression) const;
};
} // namespace SPH
#endif
