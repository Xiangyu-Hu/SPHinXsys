#ifndef ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_H
#define ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace electromagnetics
{
/**
 * @brief Register and initialize variables for harmonic (frequency-domain) A-phi solve.
 * The complex fields are represented by real/imaginary parts.
 */
class InitializeAphiFrequencyElectromagneticVariables : public LocalDynamics
{
  public:
    explicit InitializeAphiFrequencyElectromagneticVariables(SPHBody &sph_body,
                                                             Real default_conductivity,
                                                             Real default_rho_cp,
                                                             Real default_magnetic_reluctivity = 1.0);
    virtual ~InitializeAphiFrequencyElectromagneticVariables() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real default_conductivity_;
    Real default_rho_cp_;
    Real default_magnetic_reluctivity_;

    Vecd *vector_potential_real_, *vector_potential_imag_;
    Vecd *vector_potential_change_rate_real_, *vector_potential_change_rate_imag_;
    Vecd *electric_potential_gradient_real_, *electric_potential_gradient_imag_;
    Vecd *source_current_density_real_, *source_current_density_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    AngularVecd *vector_potential_curl_real_, *vector_potential_curl_imag_;

    Real *electric_potential_real_, *electric_potential_imag_;
    Real *electric_potential_source_real_, *electric_potential_source_imag_;
    Real *electric_potential_change_rate_real_, *electric_potential_change_rate_imag_;

    Vecd *electric_field_real_, *electric_field_imag_;
    Vecd *current_density_real_, *current_density_imag_;
    Real *joule_heat_source_, *temperature_change_rate_by_joule_;

    Real *electrical_conductivity_, *rho_cp_, *magnetic_reluctivity_;
};

/**
 * @brief Prescribe constant complex source current density.
 * J_s = J_s_real + j * J_s_imag.
 */
class PrescribedComplexSourceCurrentDensity : public LocalDynamics
{
  public:
    explicit PrescribedComplexSourceCurrentDensity(SPHBody &sph_body,
                                                   const Vecd &source_current_density_real,
                                                   const Vecd &source_current_density_imag);
    virtual ~PrescribedComplexSourceCurrentDensity() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd source_current_density_real_value_;
    Vecd source_current_density_imag_value_;
    Vecd *source_current_density_real_;
    Vecd *source_current_density_imag_;
};

/**
 * @brief Generic scalar Dirichlet constraint on a named field.
 */
class ConstrainScalarFieldByName : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainScalarFieldByName(BodyPartByParticle &body_part,
                                        const std::string &field_name,
                                        Real reference_value = 0.0);
    virtual ~ConstrainScalarFieldByName() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real reference_value_;
    Real *scalar_field_;
};

/**
 * @brief Generic vector Dirichlet constraint on a named field.
 */
class ConstrainVectorFieldByName : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainVectorFieldByName(BodyPartByParticle &body_part,
                                        const std::string &field_name,
                                        const Vecd &reference_value = ZeroData<Vecd>::value);
    virtual ~ConstrainVectorFieldByName() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd reference_value_;
    Vecd *vector_field_;
};

/**
 * @brief Constrain a named vector field to the local outward normal direction
 * on an axis-aligned box shell. This gives a simple n x A = 0 style boundary
 * approximation for outer air truncation without forcing the full vector to zero.
 */
class ConstrainVectorFieldToAxisAlignedBoxNormalByName : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit ConstrainVectorFieldToAxisAlignedBoxNormalByName(BodyPartByParticle &body_part,
                                                              const std::string &field_name,
                                                              const Vecd &box_center,
                                                              const Vecd &box_halfsize);
    virtual ~ConstrainVectorFieldToAxisAlignedBoxNormalByName() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd box_center_;
    Vecd box_halfsize_;
    Vecd *positions_;
    Vecd *vector_field_;
};

/**
 * @brief Source term builder for scalar potential equation from a named vector field:
 * S_phi = -div(sigma * scale * V)
 */
class ElectricPotentialSourceFromVectorFieldInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ElectricPotentialSourceFromVectorFieldInner(BaseInnerRelation &inner_relation,
                                                         const std::string &vector_field_name,
                                                         const std::string &source_name,
                                                         Real vector_field_scaling = 1.0,
                                                         Real divergence_scaling = 1.0);
    virtual ~ElectricPotentialSourceFromVectorFieldInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real vector_field_scaling_;
    Real divergence_scaling_;
    Real *Vol_, *electrical_conductivity_, *electric_potential_source_;
    Vecd *vector_field_;
};

/**
 * @brief Contact contribution for ElectricPotentialSourceFromVectorFieldInner.
 */
class ElectricPotentialSourceFromVectorFieldContact : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ElectricPotentialSourceFromVectorFieldContact(BaseContactRelation &contact_relation,
                                                           const std::string &vector_field_name,
                                                           const std::string &source_name,
                                                           Real vector_field_scaling = 1.0,
                                                           Real divergence_scaling = 1.0);
    virtual ~ElectricPotentialSourceFromVectorFieldContact() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real vector_field_scaling_;
    Real divergence_scaling_;
    Real *electrical_conductivity_, *electric_potential_source_;
    Vecd *vector_field_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electrical_conductivity_;
    StdVec<Vecd *> contact_vector_field_;
};

/**
 * @brief Generic scalar relaxation inner operator:
 * dphi/dtau = div(sigma grad(phi)) + source.
 */
class ScalarRelaxationInnerByName : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ScalarRelaxationInnerByName(BaseInnerRelation &inner_relation,
                                         const std::string &potential_name,
                                         const std::string &source_name,
                                         const std::string &change_rate_name,
                                         Real laplacian_scaling = 1.0);
    virtual ~ScalarRelaxationInnerByName() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    Real laplacian_scaling_;
    Real *Vol_, *electrical_conductivity_;
    Real *scalar_potential_, *scalar_source_, *scalar_change_rate_;
};

/**
 * @brief Generic scalar relaxation contact operator.
 */
class ScalarRelaxationContactByName : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ScalarRelaxationContactByName(BaseContactRelation &contact_relation,
                                           const std::string &potential_name,
                                           const std::string &change_rate_name,
                                           Real laplacian_scaling = 1.0);
    virtual ~ScalarRelaxationContactByName() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    Real laplacian_scaling_;
    Real *electrical_conductivity_;
    Real *scalar_change_rate_, *scalar_potential_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electrical_conductivity_;
    StdVec<Real *> contact_scalar_potential_;
};

/**
 * @brief Jacobi-type scalar relaxation using the full inner/contact diagonal:
 * dphi/dtau = [div(sigma grad(phi)) + source] / diag(L_phi)
 * The change rate is normalized by the current pseudo-time step so that
 * UpdateScalarByRelaxationRateByName still applies phi += dt * change_rate.
 */
class ScalarRelaxationComplexByName
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit ScalarRelaxationComplexByName(BaseInnerRelation &inner_relation,
                                           BaseContactRelation &contact_relation,
                                           const std::string &potential_name,
                                           const std::string &source_name,
                                           const std::string &change_rate_name,
                                           Real laplacian_scaling = 1.0,
                                           bool use_contact = true);
    virtual ~ScalarRelaxationComplexByName() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    Real laplacian_scaling_;
    bool use_contact_;
    Real *Vol_, *electrical_conductivity_;
    Real *scalar_potential_, *scalar_source_, *scalar_change_rate_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_electrical_conductivity_;
    StdVec<Real *> contact_scalar_potential_;
};

/**
 * @brief Generic scalar update from pseudo-time change rate.
 */
class UpdateScalarByRelaxationRateByName : public LocalDynamics
{
  public:
    explicit UpdateScalarByRelaxationRateByName(SPHBody &sph_body,
                                                const std::string &scalar_name,
                                                const std::string &change_rate_name,
                                                Real relaxation_scaling = 1.0,
                                                Real max_abs_value = 1.0e6);
    virtual ~UpdateScalarByRelaxationRateByName() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real relaxation_scaling_;
    Real max_abs_value_;
    Real *scalar_field_, *scalar_change_rate_;
};

/**
 * @brief Generic scalar gradient inner operator.
 */
class ScalarGradientInnerByName : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ScalarGradientInnerByName(BaseInnerRelation &inner_relation,
                                       const std::string &scalar_name,
                                       const std::string &gradient_name,
                                       Real gradient_scaling = 1.0);
    virtual ~ScalarGradientInnerByName() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real gradient_scaling_;
    Real *Vol_, *scalar_field_;
    Vecd *scalar_gradient_;
};

/**
 * @brief Generic scalar gradient contact operator.
 */
class ScalarGradientContactByName : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ScalarGradientContactByName(BaseContactRelation &contact_relation,
                                         const std::string &scalar_name,
                                         const std::string &gradient_name,
                                         Real gradient_scaling = 1.0);
    virtual ~ScalarGradientContactByName() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real gradient_scaling_;
    Real *scalar_field_;
    Vecd *scalar_gradient_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_scalar_field_;
};

/**
 * @brief Frequency-domain A-equation relaxation for one complex component:
 * dA/dtau = [ J_s - curl(nu curl A) - sigma*grad(phi) + s_omega * sigma*omega*A_coupled ] / sigma
 * where s_omega = +1 for real equation and -1 for imaginary equation.
 */
class VectorPotentialFrequencyEquationInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialFrequencyEquationInner(BaseInnerRelation &inner_relation,
                                                   Real angular_frequency,
                                                   Real omega_coupling_sign,
                                                   const std::string &vector_potential_name,
                                                   const std::string &coupled_vector_potential_name,
                                                   const std::string &source_current_density_name,
                                                   const std::string &electric_potential_gradient_name,
                                                   const std::string &curl_nu_b_name,
                                                   const std::string &change_rate_name,
                                                   Real sigma_relaxation_scaling = 1.0,
                                                   Real sigma_relaxation_floor = TinyReal,
                                                   Real magnetic_diagonal_scaling = 1.0,
                                                   Real reference_pseudo_time_step = 1.0,
                                                   Real relaxation_scaling = 1.0,
                                                   Real max_change_rate = 1.0e6);
    virtual ~VectorPotentialFrequencyEquationInner() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real omega_coupling_sign_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;

    Real *Vol_, *electrical_conductivity_, *magnetic_reluctivity_;
    Vecd *vector_potential_, *coupled_vector_potential_;
    Vecd *source_current_density_, *electric_potential_gradient_;
    Vecd *curl_nu_b_, *vector_potential_change_rate_;
};

/**
 * @brief Frequency-domain coupled A-equation relaxation for conductor regions:
 * Solve (A_re, A_im) together with local implicit omega-coupling:
 * dA_re/dtau = [J_re - curl_re - sigma*grad(phi_re) + sigma*omega*A_im] / sigma_eff
 * dA_im/dtau = [J_im - curl_im - sigma*grad(phi_im) - sigma*omega*A_re] / sigma_eff
 * The explicit curl/grad/source terms are kept as in SPH discretization, while
 * the +/- sigma*omega cross-coupling is handled implicitly per particle.
 */
class VectorPotentialFrequencyCoupledEquationInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialFrequencyCoupledEquationInner(BaseInnerRelation &inner_relation,
                                                          Real angular_frequency,
                                                          const std::string &vector_potential_real_name,
                                                          const std::string &vector_potential_imag_name,
                                                          const std::string &source_current_density_real_name,
                                                          const std::string &source_current_density_imag_name,
                                                          const std::string &electric_potential_gradient_real_name,
                                                          const std::string &electric_potential_gradient_imag_name,
                                                          const std::string &curl_nu_b_real_name,
                                                          const std::string &curl_nu_b_imag_name,
                                                          const std::string &change_rate_real_name,
                                                          const std::string &change_rate_imag_name,
                                                          Real sigma_relaxation_scaling = 1.0,
                                                          Real sigma_relaxation_floor = TinyReal,
                                                          Real magnetic_diagonal_scaling = 1.0,
                                                          Real reference_pseudo_time_step = 1.0,
                                                          Real relaxation_scaling = 1.0,
                                                          Real max_change_rate = 1.0e6);
    virtual ~VectorPotentialFrequencyCoupledEquationInner() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;

    Real *Vol_, *electrical_conductivity_, *magnetic_reluctivity_;
    Vecd *vector_potential_real_, *vector_potential_imag_;
    Vecd *source_current_density_real_, *source_current_density_imag_;
    Vecd *electric_potential_gradient_real_, *electric_potential_gradient_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    Vecd *vector_potential_change_rate_real_, *vector_potential_change_rate_imag_;
};

/**
 * @brief Frequency-domain conductor A-equation with inner + contact diagonal estimate.
 * The residual still uses the already-assembled CurlNuB / grad(phi) fields, while
 * the local magnetic diagonal includes both inner and contact neighbors.
 */
class VectorPotentialFrequencyEquationComplex
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit VectorPotentialFrequencyEquationComplex(BaseInnerRelation &inner_relation,
                                                     BaseContactRelation &contact_relation,
                                                     Real angular_frequency,
                                                     Real omega_coupling_sign,
                                                     const std::string &vector_potential_name,
                                                     const std::string &coupled_vector_potential_name,
                                                     const std::string &source_current_density_name,
                                                     const std::string &electric_potential_gradient_name,
                                                     const std::string &curl_nu_b_name,
                                                     const std::string &change_rate_name,
                                                     Real sigma_relaxation_scaling = 1.0,
                                                     Real sigma_relaxation_floor = TinyReal,
                                                     Real magnetic_diagonal_scaling = 1.0,
                                                     Real reference_pseudo_time_step = 1.0,
                                                     Real relaxation_scaling = 1.0,
                                                     Real max_change_rate = 1.0e6,
                                                     Real balanced_magnetic_diagonal_weight = 0.0,
                                                     Real contact_diagonal_ratio_cap = 0.0,
                                                     Real adaptive_contact_diagonal_cap_strength = 0.0);
    virtual ~VectorPotentialFrequencyEquationComplex() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real omega_coupling_sign_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real balanced_magnetic_diagonal_weight_;
    Real contact_diagonal_ratio_cap_;
    Real adaptive_contact_diagonal_cap_strength_;

    Real *Vol_, *electrical_conductivity_, *magnetic_reluctivity_;
    Vecd *vector_potential_, *coupled_vector_potential_;
    Vecd *source_current_density_, *electric_potential_gradient_;
    Vecd *curl_nu_b_, *vector_potential_change_rate_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
};

/**
 * @brief Frequency-domain coupled conductor A-equation with inner + contact diagonal estimate.
 */
class VectorPotentialFrequencyCoupledEquationComplex
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit VectorPotentialFrequencyCoupledEquationComplex(BaseInnerRelation &inner_relation,
                                                            BaseContactRelation &contact_relation,
                                                            Real angular_frequency,
                                                            const std::string &vector_potential_real_name,
                                                            const std::string &vector_potential_imag_name,
                                                            const std::string &source_current_density_real_name,
                                                            const std::string &source_current_density_imag_name,
                                                            const std::string &electric_potential_gradient_real_name,
                                                            const std::string &electric_potential_gradient_imag_name,
                                                            const std::string &curl_nu_b_real_name,
                                                            const std::string &curl_nu_b_imag_name,
                                                            const std::string &change_rate_real_name,
                                                            const std::string &change_rate_imag_name,
                                                            Real sigma_relaxation_scaling = 1.0,
                                                            Real sigma_relaxation_floor = TinyReal,
                                                            Real magnetic_diagonal_scaling = 1.0,
                                                            Real reference_pseudo_time_step = 1.0,
                                                            Real relaxation_scaling = 1.0,
                                                            Real max_change_rate = 1.0e6,
                                                            Real balanced_magnetic_diagonal_weight = 0.0,
                                                            Real contact_diagonal_ratio_cap = 0.0,
                                                            Real adaptive_contact_diagonal_cap_strength = 0.0);
    virtual ~VectorPotentialFrequencyCoupledEquationComplex() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real balanced_magnetic_diagonal_weight_;
    Real contact_diagonal_ratio_cap_;
    Real adaptive_contact_diagonal_cap_strength_;

    Real *Vol_, *electrical_conductivity_, *magnetic_reluctivity_;
    Vecd *vector_potential_real_, *vector_potential_imag_;
    Vecd *source_current_density_real_, *source_current_density_imag_;
    Vecd *electric_potential_gradient_real_, *electric_potential_gradient_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    Vecd *vector_potential_change_rate_real_, *vector_potential_change_rate_imag_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
};

/**
 * @brief Frequency-domain coupled conductor A-equation with a local vector magnetic self-block.
 * The curl(nu curl A) contribution is treated with a local Matd block, while the
 * +/- omega*sigma coupling between real/imag parts is solved in one small 2*Dim system.
 * A conservative scalar coupled solve is kept as fallback and as a norm limiter.
 */
class VectorPotentialFrequencyCoupledBlockEquationComplex
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit VectorPotentialFrequencyCoupledBlockEquationComplex(BaseInnerRelation &inner_relation,
                                                                 BaseContactRelation &contact_relation,
                                                                 Real angular_frequency,
                                                                 const std::string &vector_potential_real_name,
                                                                 const std::string &vector_potential_imag_name,
                                                                 const std::string &source_current_density_real_name,
                                                                 const std::string &source_current_density_imag_name,
                                                                 const std::string &electric_potential_gradient_real_name,
                                                                 const std::string &electric_potential_gradient_imag_name,
                                                                 const std::string &curl_nu_b_real_name,
                                                                 const std::string &curl_nu_b_imag_name,
                                                                 const std::string &change_rate_real_name,
                                                                 const std::string &change_rate_imag_name,
                                                                 Real sigma_relaxation_scaling = 1.0,
                                                                 Real sigma_relaxation_floor = TinyReal,
                                                                 Real magnetic_diagonal_scaling = 1.0,
                                                                 Real reference_pseudo_time_step = 1.0,
                                                                 Real relaxation_scaling = 1.0,
                                                                 Real max_change_rate = 1.0e6,
                                                                 Real balanced_magnetic_diagonal_weight = 0.0,
                                                                 Real contact_diagonal_ratio_cap = 0.0,
                                                                 Real adaptive_contact_diagonal_cap_strength = 0.0);
    virtual ~VectorPotentialFrequencyCoupledBlockEquationComplex() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real balanced_magnetic_diagonal_weight_;
    Real contact_diagonal_ratio_cap_;
    Real adaptive_contact_diagonal_cap_strength_;

    Real *Vol_, *electrical_conductivity_, *magnetic_reluctivity_;
    Vecd *vector_potential_real_, *vector_potential_imag_;
    Vecd *source_current_density_real_, *source_current_density_imag_;
    Vecd *electric_potential_gradient_real_, *electric_potential_gradient_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    Vecd *vector_potential_change_rate_real_, *vector_potential_change_rate_imag_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
};

/**
 * @brief Frequency-domain magnetic-only A-equation with inner + contact diagonal estimate.
 * Intended for impressed-current multiturn coils where the coil current is prescribed
 * and no isotropic conductor reaction term is solved inside the coil domain.
 */
class VectorPotentialFrequencyMagneticOnlyEquationComplex
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit VectorPotentialFrequencyMagneticOnlyEquationComplex(BaseInnerRelation &inner_relation,
                                                                 BaseContactRelation &contact_relation,
                                                                 const std::string &vector_potential_name,
                                                                 const std::string &source_current_density_name,
                                                                 const std::string &curl_nu_b_name,
                                                                 const std::string &change_rate_name,
                                                                 Real magnetic_diagonal_scaling = 1.0,
                                                                 Real reference_pseudo_time_step = 1.0,
                                                                 Real relaxation_scaling = 1.0,
                                                                 Real max_change_rate = 1.0e6,
                                                                 Real balanced_magnetic_diagonal_weight = 0.0,
                                                                 Real contact_diagonal_ratio_cap = 0.0,
                                                                 Real adaptive_contact_diagonal_cap_strength = 0.0);
    virtual ~VectorPotentialFrequencyMagneticOnlyEquationComplex() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real balanced_magnetic_diagonal_weight_;
    Real contact_diagonal_ratio_cap_;
    Real adaptive_contact_diagonal_cap_strength_;

    Real *Vol_, *magnetic_reluctivity_;
    Vecd *vector_potential_;
    Vecd *source_current_density_;
    Vecd *curl_nu_b_, *vector_potential_change_rate_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
};

/**
 * @brief Frequency-domain magnetic-only A-equation with a local vector self-block solve.
 * The update uses the self-block induced by the current curl(nu curl A) stencil and
 * falls back to the conservative scalar diagonal when the local block solve is not safer.
 * This is intended for insulating regions such as air, where point-wise scalar Jacobi
 * can be too weak for the contact-coupled magnetic propagation.
 */
class VectorPotentialFrequencyMagneticOnlyBlockEquationComplex
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    explicit VectorPotentialFrequencyMagneticOnlyBlockEquationComplex(BaseInnerRelation &inner_relation,
                                                                      BaseContactRelation &contact_relation,
                                                                      const std::string &vector_potential_name,
                                                                      const std::string &source_current_density_name,
                                                                      const std::string &curl_nu_b_name,
                                                                      const std::string &change_rate_name,
                                                                      Real magnetic_diagonal_scaling = 1.0,
                                                                      Real reference_pseudo_time_step = 1.0,
                                                                      Real relaxation_scaling = 1.0,
                                                                      Real max_change_rate = 1.0e6,
                                                                      Real balanced_magnetic_diagonal_weight = 0.0,
                                                                      Real contact_diagonal_ratio_cap = 0.0,
                                                                      Real adaptive_contact_diagonal_cap_strength = 0.0);
    virtual ~VectorPotentialFrequencyMagneticOnlyBlockEquationComplex() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;
    Real balanced_magnetic_diagonal_weight_;
    Real contact_diagonal_ratio_cap_;
    Real adaptive_contact_diagonal_cap_strength_;

    Real *Vol_, *magnetic_reluctivity_;
    Vecd *vector_potential_;
    Vecd *source_current_density_;
    Vecd *curl_nu_b_, *vector_potential_change_rate_;
    StdVec<Real *> contact_vol_;
    StdVec<Real *> contact_magnetic_reluctivity_;
};

/**
 * @brief Frequency-domain magnetic-only A-equation relaxation for insulating regions:
 * dA/dtau = J_s - curl(nu curl A)
 * This avoids dividing by sigma and avoids sigma*omega / sigma*grad(phi) terms in air.
 */
class VectorPotentialFrequencyMagneticOnlyEquationInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialFrequencyMagneticOnlyEquationInner(BaseInnerRelation &inner_relation,
                                                               const std::string &vector_potential_name,
                                                               const std::string &source_current_density_name,
                                                               const std::string &curl_nu_b_name,
                                                               const std::string &change_rate_name,
                                                               Real magnetic_diagonal_scaling = 1.0,
                                                               Real reference_pseudo_time_step = 1.0,
                                                               Real relaxation_scaling = 1.0,
                                                               Real max_change_rate = 1.0e6);
    virtual ~VectorPotentialFrequencyMagneticOnlyEquationInner() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real magnetic_diagonal_scaling_;
    Real reference_pseudo_time_step_;
    Real relaxation_scaling_;
    Real max_change_rate_;

    Real *Vol_, *magnetic_reluctivity_;
    Vecd *vector_potential_;
    Vecd *source_current_density_;
    Vecd *curl_nu_b_;
    Vecd *vector_potential_change_rate_;
};

/**
 * @brief Diagnostic residual evaluation for frequency-domain A equations.
 * Computes per-particle residual vector and a normalized residual magnitude.
 */
class FrequencyAEquationResidualDiagnostic : public LocalDynamics
{
  public:
    explicit FrequencyAEquationResidualDiagnostic(
        SPHBody &sph_body,
        Real angular_frequency,
        Real omega_coupling_sign,
        const std::string &a_name,
        const std::string &coupled_a_name,
        const std::string &source_name,
        const std::string &grad_phi_name,
        const std::string &curl_nu_b_name,
        const std::string &residual_name,
        const std::string &relative_residual_name);
    virtual ~FrequencyAEquationResidualDiagnostic() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real omega_coupling_sign_;
    Real *electrical_conductivity_;
    Vecd *a_, *coupled_a_;
    Vecd *source_, *grad_phi_, *curl_nu_b_;
    Vecd *residual_;
    Real *relative_residual_;
};

/**
 * @brief Compute RMS Joule source from complex fields:
 * E_re =  omega * A_im - grad(phi_re)
 * E_im = -omega * A_re - grad(phi_im)
 * J_re = sigma * E_re
 * J_im = sigma * E_im
 * Q = 0.5 * (J_re dot E_re + J_im dot E_im)
 */
class FrequencyElectricFieldCurrentAndJouleHeatInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit FrequencyElectricFieldCurrentAndJouleHeatInner(BaseInnerRelation &inner_relation,
                                                            Real angular_frequency,
                                                            const std::string &a_real_name = "VectorPotentialReal",
                                                            const std::string &a_imag_name = "VectorPotentialImag",
                                                            const std::string &grad_phi_real_name = "ElectricPotentialGradientReal",
                                                            const std::string &grad_phi_imag_name = "ElectricPotentialGradientImag",
                                                            const std::string &e_real_name = "ElectricFieldReal",
                                                            const std::string &e_imag_name = "ElectricFieldImag",
                                                            const std::string &j_real_name = "CurrentDensityReal",
                                                            const std::string &j_imag_name = "CurrentDensityImag",
                                                            const std::string &joule_name = "JouleHeatSource",
                                                            const std::string &temperature_change_rate_name = "TemperatureChangeRateByJoule");
    virtual ~FrequencyElectricFieldCurrentAndJouleHeatInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real *electrical_conductivity_, *rho_cp_;
    Vecd *a_real_, *a_imag_;
    Vecd *grad_phi_real_, *grad_phi_imag_;
    Vecd *electric_field_real_, *electric_field_imag_;
    Vecd *current_density_real_, *current_density_imag_;
    Real *joule_heat_source_, *temperature_change_rate_by_joule_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_H
