#ifndef ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_H
#define ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace electromagnetics
{
/**
 * @brief Register and initialize A-phi electromagnetic variables for TEAM7 style cases.
 * This class is a simple state updater and can be called once before time stepping.
 */
class InitializeAphiElectromagneticVariables : public LocalDynamics
{
  public:
    explicit InitializeAphiElectromagneticVariables(SPHBody &sph_body,
                                                    Real default_conductivity,
                                                    Real default_rho_cp,
                                                    Real default_magnetic_reluctivity = 1.0);
    virtual ~InitializeAphiElectromagneticVariables() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real default_conductivity_;
    Real default_rho_cp_;
    Real default_magnetic_reluctivity_;
    Vecd *vector_potential_, *vector_potential_prev_, *vector_potential_dt_;
    Vecd *electric_field_, *current_density_, *electric_potential_gradient_;
    Vecd *source_current_density_, *vector_potential_change_rate_, *curl_nu_b_;
    AngularVecd *vector_potential_curl_;
    Real *electric_potential_, *electric_potential_source_, *electric_potential_change_rate_;
    Real *joule_heat_source_, *temperature_change_rate_by_joule_;
    Real *electrical_conductivity_, *rho_cp_, *magnetic_reluctivity_;
};

/**
 * @brief Compute A_dot = dA/dt using a first-order backward difference.
 */
class UpdateVectorPotentialTimeDerivative : public LocalDynamics
{
  public:
    explicit UpdateVectorPotentialTimeDerivative(SPHBody &sph_body);
    virtual ~UpdateVectorPotentialTimeDerivative() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vector_potential_, *vector_potential_prev_, *vector_potential_dt_;
};

/**
 * @brief Compute source term for electric potential equation:
 *        S_phi = -div(sigma * A_dot)
 */
class ElectricPotentialSourceTermInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ElectricPotentialSourceTermInner(BaseInnerRelation &inner_relation);
    virtual ~ElectricPotentialSourceTermInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_, *electrical_conductivity_, *electric_potential_source_;
    Vecd *vector_potential_dt_;
};

/**
 * @brief Contact contribution of source term for electric potential equation:
 *        S_phi += -div(sigma * A_dot)|contact
 * This class should be executed after ElectricPotentialSourceTermInner.
 */
class ElectricPotentialSourceTermContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ElectricPotentialSourceTermContact(BaseContactRelation &contact_relation);
    virtual ~ElectricPotentialSourceTermContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *electrical_conductivity_, *electric_potential_source_;
    Vecd *vector_potential_dt_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electrical_conductivity_;
    StdVec<Vecd *> contact_vector_potential_dt_;
};

/**
 * @brief Enforce reference potential (Dirichlet) on a body part to avoid singularity of scalar potential equation.
 */
class ConstrainElectricPotential : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainElectricPotential(BodyPartByParticle &body_part, Real reference_value = 0.0);
    virtual ~ConstrainElectricPotential() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real reference_value_;
    Real *electric_potential_;
};

/**
 * @brief Enforce vector potential on a body part (e.g. far-field truncation A=0).
 */
class ConstrainVectorPotential : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainVectorPotential(BodyPartByParticle &body_part,
                                      const Vecd &reference_value = ZeroData<Vecd>::value);
    virtual ~ConstrainVectorPotential() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd reference_value_;
    Vecd *vector_potential_;
};

/**
 * @brief Set constant electromagnetic material parameters on one body.
 */
class SetConstantElectromagneticMaterialProperties : public LocalDynamics
{
  public:
    explicit SetConstantElectromagneticMaterialProperties(SPHBody &sph_body,
                                                          Real electrical_conductivity,
                                                          Real rho_cp,
                                                          Real magnetic_reluctivity);
    virtual ~SetConstantElectromagneticMaterialProperties() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real electrical_conductivity_value_;
    Real rho_cp_value_;
    Real magnetic_reluctivity_value_;
    Real *electrical_conductivity_, *rho_cp_, *magnetic_reluctivity_;
};

/**
 * @brief Optional temperature-dependent conductivity model:
 *        sigma(T) = sigma_ref * (1 - alpha * (T - T_ref)).
 * This class is optional for one-way coupling and can be skipped.
 */
class UpdateElectricalConductivityByLinearTemperature : public LocalDynamics
{
  public:
    explicit UpdateElectricalConductivityByLinearTemperature(SPHBody &sph_body,
                                                             Real sigma_ref,
                                                             Real alpha,
                                                             Real T_ref);
    virtual ~UpdateElectricalConductivityByLinearTemperature() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real sigma_ref_, alpha_, T_ref_;
    Real *temperature_, *electrical_conductivity_;
};

/**
 * @brief Uniform impressed source current density J_s on one body.
 */
class PrescribedSourceCurrentDensity : public LocalDynamics
{
  public:
    explicit PrescribedSourceCurrentDensity(SPHBody &sph_body, const Vecd &source_current_density);
    virtual ~PrescribedSourceCurrentDensity() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd source_current_density_value_;
    Vecd *source_current_density_;
};

/**
 * @brief Uniform impressed source current density J_s on selected body part.
 */
class PrescribedSourceCurrentDensityByBodyPart : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit PrescribedSourceCurrentDensityByBodyPart(BodyPartByParticle &body_part,
                                                      const Vecd &source_current_density);
    virtual ~PrescribedSourceCurrentDensityByBodyPart() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd source_current_density_value_;
    Vecd *source_current_density_;
};

/**
 * @brief Pseudo-time relaxation for scalar potential equation:
 *        dphi/dtau = div(sigma grad(phi)) + S_phi
 */
class ElectricPotentialRelaxationInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ElectricPotentialRelaxationInner(BaseInnerRelation &inner_relation,
                                              Real relaxation_scaling = 1.0);
    virtual ~ElectricPotentialRelaxationInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    Real *Vol_, *electrical_conductivity_;
    Real *electric_potential_, *electric_potential_source_, *electric_potential_change_rate_;
};

/**
 * @brief Contact contribution for pseudo-time relaxation of scalar potential:
 *        dphi/dtau += div(sigma grad(phi))|contact
 * This class should be executed after ElectricPotentialRelaxationInner.
 */
class ElectricPotentialRelaxationContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ElectricPotentialRelaxationContact(BaseContactRelation &contact_relation);
    virtual ~ElectricPotentialRelaxationContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    Real *electrical_conductivity_;
    Real *electric_potential_change_rate_, *electric_potential_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electrical_conductivity_;
    StdVec<Real *> contact_electric_potential_;
};

/**
 * @brief Update scalar potential from the computed pseudo-time change rate.
 */
class UpdateElectricPotentialByRelaxationRate : public LocalDynamics
{
  public:
    explicit UpdateElectricPotentialByRelaxationRate(SPHBody &sph_body,
                                                     Real relaxation_scaling = 1.0,
                                                     Real max_abs_potential = 1.0e6);
    virtual ~UpdateElectricPotentialByRelaxationRate() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real relaxation_scaling_;
    Real max_abs_potential_;
    Real *electric_potential_, *electric_potential_change_rate_;
};

/**
 * @brief Compute B = curl(A) with standard SPH gradient form.
 * Output variable is "VectorPotentialCurl" by default.
 */
class VectorPotentialCurlInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialCurlInner(BaseInnerRelation &inner_relation,
                                      const std::string &vector_potential_name = "VectorPotential",
                                      const std::string &output_name = "VectorPotentialCurl",
                                      Real curl_scaling = 1.0);
    virtual ~VectorPotentialCurlInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real curl_scaling_;
    Real *Vol_;
    Vecd *vector_potential_;
    AngularVecd *vector_potential_curl_;
};

/**
 * @brief Contact contribution to B = curl(A), accumulated on receiver body.
 * This class should be executed after VectorPotentialCurlInner.
 */
class VectorPotentialCurlContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit VectorPotentialCurlContact(BaseContactRelation &contact_relation,
                                        const std::string &vector_potential_name = "VectorPotential",
                                        const std::string &output_name = "VectorPotentialCurl",
                                        Real curl_scaling = 1.0);
    virtual ~VectorPotentialCurlContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real curl_scaling_;
    Vecd *vector_potential_;
    AngularVecd *vector_potential_curl_;
    StdVec<Real *> contact_vol_;
    StdVec<Vecd *> contact_vector_potential_;
};

/**
 * @brief Compute C = curl(nu * B), where B = curl(A), nu = magnetic reluctivity.
 * Requires B variable "VectorPotentialCurl" produced by VectorPotentialCurlInner
 * (or by AphiCurlCK if CK execution is used).
 *
 * If symmetric_volume_weight is true, each pair uses 0.5*(V_i+V_j) in place of V_j in the
 * kernel-gradient weight. This only changes results when V_i and V_j differ (e.g. adaptive
 * resolution); for uniform volumetric measures it matches the legacy weight.
 */
class CurlNuBInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit CurlNuBInner(BaseInnerRelation &inner_relation,
                          const std::string &magnetic_flux_density_name = "VectorPotentialCurl",
                          const std::string &output_name = "CurlNuB",
                          Real curl_scaling = 1.0,
                          bool symmetric_volume_weight = false);
    virtual ~CurlNuBInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real curl_scaling_;
    bool symmetric_volume_weight_;
    Real *Vol_, *magnetic_reluctivity_;
    AngularVecd *magnetic_flux_density_;
    Vecd *curl_nu_b_;
};

/**
 * @brief Contact contribution to C = curl(nu * B), accumulated on receiver body.
 * This class should be executed after CurlNuBInner.
 *
 * Split-vs-merged consistency also depends on how ContactRelation builds neighborhoods
 * (target-body mesh and search depth) versus InnerRelation on a single merged body; the
 * discrete pairing here follows the same grad-W form as CurlNuBInner when symmetric_volume_weight
 * is enabled (0.5*(V_i+V_j^contact) for cross-body pairs).
 */
class CurlNuBContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit CurlNuBContact(BaseContactRelation &contact_relation,
                            const std::string &magnetic_flux_density_name = "VectorPotentialCurl",
                            const std::string &output_name = "CurlNuB",
                            Real curl_scaling = 1.0,
                            bool symmetric_volume_weight = false);
    virtual ~CurlNuBContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real curl_scaling_;
    bool symmetric_volume_weight_;
    Real *Vol_;
    Real *magnetic_reluctivity_;
    AngularVecd *magnetic_flux_density_;
    Vecd *curl_nu_b_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
    StdVec<AngularVecd *> contact_magnetic_flux_density_;
};

/**
 * @brief Compute grad(phi) with inner-neighbor contribution.
 */
class ElectricPotentialGradientInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ElectricPotentialGradientInner(BaseInnerRelation &inner_relation);
    virtual ~ElectricPotentialGradientInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_, *electric_potential_;
    Vecd *electric_potential_gradient_;
};

/**
 * @brief Add contact-neighbor contribution to grad(phi).
 * This class should be executed after ElectricPotentialGradientInner.
 */
class ElectricPotentialGradientContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ElectricPotentialGradientContact(BaseContactRelation &contact_relation);
    virtual ~ElectricPotentialGradientContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *electric_potential_;
    Vecd *electric_potential_gradient_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electric_potential_;
};

/**
 * @brief A-equation update (time-domain explicit form):
 *        sigma * dA/dt = J_s - curl(nu * B) - sigma * grad(phi)
 */
class VectorPotentialEquationInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialEquationInner(BaseInnerRelation &inner_relation,
                                          Real relaxation_scaling = 1.0,
                                          Real max_change_rate = 1.0e6);
    virtual ~VectorPotentialEquationInner() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real *electrical_conductivity_;
    Vecd *vector_potential_, *source_current_density_, *electric_potential_gradient_;
    Vecd *curl_nu_b_, *vector_potential_change_rate_;
};

/**
 * @brief Magnetic-only A update for very low-conductivity regions (e.g. air):
 *        dA/dtau = -curl(nu * B)
 * This avoids dividing by near-zero sigma in pseudo-time iterations.
 */
class VectorPotentialEquationMagneticOnlyInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialEquationMagneticOnlyInner(BaseInnerRelation &inner_relation,
                                                      Real relaxation_scaling = 1.0,
                                                      Real max_change_rate = 1.0e6);
    virtual ~VectorPotentialEquationMagneticOnlyInner() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real relaxation_scaling_;
    Real max_change_rate_;
    Vecd *vector_potential_, *curl_nu_b_, *vector_potential_change_rate_;
};

/**
 * @brief Compute electromagnetic fields and Joule source:
 *        E = -A_dot - grad(phi), J = sigma E, Q = J dot E.
 * Also provides temperature source rate: dT/dt|Joule = Q / (rho*cp).
 */
class ElectricFieldCurrentAndJouleHeatInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ElectricFieldCurrentAndJouleHeatInner(BaseInnerRelation &inner_relation,
                                                   bool recompute_potential_gradient = true);
    virtual ~ElectricFieldCurrentAndJouleHeatInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    bool recompute_potential_gradient_;
    Real *Vol_, *electrical_conductivity_, *rho_cp_;
    Real *electric_potential_, *joule_heat_source_, *temperature_change_rate_by_joule_;
    Vecd *vector_potential_dt_, *electric_field_, *current_density_, *electric_potential_gradient_;
};

/**
 * @brief Add Joule-induced temperature change rate into an existing diffusion species change rate.
 * For temperature diffusion species "Temperature", target rate variable is typically "TemperatureChangeRate".
 */
class AddJouleHeatToDiffusionSpeciesRate : public LocalDynamics
{
  public:
    explicit AddJouleHeatToDiffusionSpeciesRate(SPHBody &sph_body,
                                                const std::string &target_species_change_rate_name = "TemperatureChangeRate");
    virtual ~AddJouleHeatToDiffusionSpeciesRate() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *temperature_change_rate_by_joule_;
    Real *target_species_change_rate_;
};
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_H
