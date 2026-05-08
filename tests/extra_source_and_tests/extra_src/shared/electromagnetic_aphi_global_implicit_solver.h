#ifndef ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_H
#define ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace electromagnetics
{
/**
 * @brief Local block preconditioner for the coupled frequency-domain A equation.
 * It consumes externally assembled residual fields and produces a search direction.
 */
class VectorPotentialFrequencyCoupledPreconditionerInner
    : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VectorPotentialFrequencyCoupledPreconditionerInner(
        BaseInnerRelation &inner_relation,
        Real angular_frequency,
        const std::string &residual_real_name,
        const std::string &residual_imag_name,
        const std::string &search_real_name,
        const std::string &search_imag_name,
        Real sigma_relaxation_scaling,
        Real sigma_relaxation_floor,
        Real magnetic_diagonal_scaling,
        Real max_search_norm);
    virtual ~VectorPotentialFrequencyCoupledPreconditionerInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real sigma_relaxation_scaling_;
    Real sigma_relaxation_floor_;
    Real magnetic_diagonal_scaling_;
    Real max_search_norm_;

    Real *Vol_;
    Real *electrical_conductivity_;
    Real *magnetic_reluctivity_;
    Vecd *residual_real_;
    Vecd *residual_imag_;
    Vecd *search_real_;
    Vecd *search_imag_;
};

/**
 * @brief Assemble the linear operator action L(z) for the coupled frequency-domain A system.
 * Residual is defined as r = b - L(A), so the operator action is used in
 * r_new = r - alpha * L(z).
 */
class FrequencyVectorPotentialLinearOperatorComplex : public LocalDynamics
{
  public:
    explicit FrequencyVectorPotentialLinearOperatorComplex(
        SPHBody &sph_body,
        Real angular_frequency,
        const std::string &search_real_name,
        const std::string &search_imag_name,
        const std::string &curl_nu_b_real_name,
        const std::string &curl_nu_b_imag_name,
        const std::string &operator_action_real_name,
        const std::string &operator_action_imag_name);
    virtual ~FrequencyVectorPotentialLinearOperatorComplex() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real *electrical_conductivity_;
    Vecd *search_real_;
    Vecd *search_imag_;
    Vecd *curl_nu_b_real_;
    Vecd *curl_nu_b_imag_;
    Vecd *operator_action_real_;
    Vecd *operator_action_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_H
