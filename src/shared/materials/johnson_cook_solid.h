#ifndef JOHNSON_COOK_SOLID_H
#define JOHNSON_COOK_SOLID_H

#include "elastic_solid.h"
#include <memory>

namespace SPH
{
/** Johnson-Cook flow, thermal and dynamic-failure parameters in consistent units. */
struct JohnsonCookParameters
{
    Real yield_stress, hardening_coefficient, hardening_exponent;
    Real strain_rate_coefficient, reference_strain_rate;
    Real thermal_softening_exponent, reference_temperature, melting_temperature;
    Real specific_heat;
    Real failure_D1, failure_D2, failure_D3, failure_D4, failure_D5;
    bool failure_enabled = true;
};

/** Three-dimensional stress history; two-dimensional bodies use plane strain. */
struct JohnsonCookState
{
    Mat3d stress = Mat3d::Zero();
    Real plastic_strain = 0;
    Real plastic_strain_rate = 0;
    Real temperature = 0;
    Real damage = 0;
    Real equivalent_stress = 0;
    Real mean_stress = 0;
    Real plastic_work_specific = 0;
    bool active = true;
};

enum class JohnsonCookIntegrationStatus
{
    success = 0,
    invalid_increment,
    invalid_deformation,
    singular_increment,
    nonfinite_state
};

/** Value-only model shared by native CPU and device computing kernels. */
class JohnsonCookConstitutiveModel
{
  public:
    JohnsonCookConstitutiveModel(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                 const JohnsonCookParameters &parameters)
        : parameters_(parameters), rho0_(rho0), G0_(youngs_modulus / (2 * (1 + poisson_ratio))),
          K0_(youngs_modulus / (3 * (1 - 2 * poisson_ratio))),
          inverse_heat_capacity_(1 / (rho0 * parameters.specific_heat)),
          inverse_temperature_interval_(1 / (parameters.melting_temperature - parameters.reference_temperature)) {}
    Real NormalizedTemperature(Real temperature) const;
    Real FractureStrain(Real mean_stress_ratio, Real strain_rate, Real temperature) const;
    JohnsonCookIntegrationStatus Integrate(const Mat3d &increment, Real J, Real dt, JohnsonCookState &state) const;
    JohnsonCookIntegrationStatus UpdateStress(const Matd &F, const Matd &previous_F, Real dt, JohnsonCookState &state) const;

  protected:
    JohnsonCookParameters parameters_;
    Real rho0_, G0_, K0_, inverse_heat_capacity_, inverse_temperature_interval_;
};

/** Incremental Johnson-Cook solid with adiabatic heating and complete failure.
 * Flow stress uses the preceding plastic rate. The failure law uses the current
 * deviatoric total rate and preceding temperature. Mean stress is tension-positive.
 * Use with solid_dynamics::JohnsonCookIntegration1stHalf and its second half.
 */
class JohnsonCookSolid : public ElasticSolid
{
  public:
    JohnsonCookSolid(Real rho0, Real youngs_modulus, Real poisson_ratio,
                     const JohnsonCookParameters &parameters);
    void initializeLocalParameters(BaseParticles *particles) override;
    const JohnsonCookParameters &Parameters() const { return parameters_; }
    JohnsonCookState InitialState() const;
    Real NormalizedTemperature(Real temperature) const;
    Real FractureStrain(Real mean_stress_ratio, Real strain_rate, Real temperature) const;
    /** Explicit radial return; J is the current relative volume. */
    void Integrate(const Mat3d &strain_increment, Real J, Real dt, JohnsonCookState &state) const;
    Matd UpdateStress(const Matd &deformation, size_t index_i, Real dt);

    const JohnsonCookConstitutiveModel &ConstitutiveModel() const { return model_; }

    class ConstituteKernel
    {
      public:
        template <class ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, JohnsonCookSolid &encloser);
        Matd UpdateStress(const Matd &F, size_t i, Real dt);
        int Status(size_t i) const { return status_[i]; }

      protected:
        JohnsonCookConstitutiveModel model_;
        Mat3d *stress_;
        Matd *previous_deformation_;
        Real *previous_time_step_, *plastic_strain_, *plastic_strain_rate_;
        Real *temperature_, *damage_, *equivalent_stress_, *mean_stress_, *plastic_work_specific_;
        int *material_active_, *status_;
    };

    Matd StressPK1(Matd &deformation, size_t index_i) override;
    Matd StressPK2(Matd &deformation, size_t index_i) override;
    Matd StressCauchy(Matd &strain, size_t index_i) override;
    Real VolumetricKirchhoff(Real J) override;
    std::string getRelevantStressMeasureName() override { return "Cauchy"; }

  protected:
    JohnsonCookParameters parameters_;
    JohnsonCookConstitutiveModel model_;
    BaseParticles *particles_ = nullptr;
    std::unique_ptr<ConstituteKernel> host_kernel_;
    Mat3d *stress_;
    Matd *previous_deformation_;
    Real *previous_time_step_, *plastic_strain_, *plastic_strain_rate_;
    Real *temperature_, *damage_, *equivalent_stress_, *mean_stress_, *plastic_work_specific_;
    int *material_active_;
};
} // namespace SPH
#endif // JOHNSON_COOK_SOLID_H
