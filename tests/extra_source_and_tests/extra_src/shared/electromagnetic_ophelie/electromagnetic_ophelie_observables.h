#ifndef ELECTROMAGNETIC_OPHELIE_OBSERVABLES_H
#define ELECTROMAGNETIC_OPHELIE_OBSERVABLES_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real hostVolWeightedSum(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real weighted_sum = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        weighted_sum += vol[i] * values[i];
    }
    return weighted_sum;
}

inline Real hostVecdFieldMax(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_value = std::max(max_value, values[i].norm());
    }
    return max_value;
}

inline Real hostScalarFieldMax(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_value = std::max(max_value, std::abs(values[i]));
    }
    return max_value;
}

inline Real hostVolWeightedDot(BaseParticles &particles, const Real *lhs_values, const Real *rhs_values,
                               size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real dot_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        dot_value += vol[i] * lhs_values[i] * rhs_values[i];
    }
    return dot_value;
}

inline Real hostVolWeightedNorm(BaseParticles &particles, const Real *values, size_t total_real_particles)
{
    return std::sqrt(std::max(hostVolWeightedDot(particles, values, values, total_real_particles), Real(0)));
}

inline void hostAssignScalarField(BaseParticles &particles, const std::string &variable_name, const Real *values,
                                 size_t total_real_particles)
{
    Real *field = particles.getVariableDataByName<Real>(variable_name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        field[i] = values[i];
    }
    syncVariableToDevice<Real>(particles, variable_name);
}

inline void hostReadScalarField(BaseParticles &particles, const std::string &variable_name, Real *values,
                              size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *field = particles.getVariableDataByName<Real>(variable_name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        values[i] = field[i];
    }
}

inline void hostSubtractScaledVector(Real *values, const Real *direction, Real scale, size_t total_real_particles)
{
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        values[i] -= scale * direction[i];
    }
}

inline void hostScaledAdd(Real *values, Real scale, const Real *direction, size_t total_real_particles)
{
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        values[i] += scale * direction[i];
    }
}

struct OphelieRunMetrics
{
    size_t n_glass = 0;
    size_t n_coil = 0;
    Real joule_power_raw = 0.0;
    Real joule_power_scaled = 0.0;
    Real power_scale = 1.0;
    Real field_scale = 1.0;
    Real effective_current_amplitude = 0.0;
    Real div_j_rel_level0 = 0.0;
    Real div_j_rel_phi = 0.0;
    Real div_j_reduction = 0.0;
    Real max_a_coil = 0.0;
    Real max_a_ind = 0.0;
    Real max_a_src = 0.0;
    Real max_b_src = 0.0;
    Real self_induction_j_rel_change = 0.0;
    size_t self_induction_iterations_used = 0;
    Real max_e_imag = 0.0;
    Real max_j_imag = 0.0;
    Real max_joule_heat = 0.0;
    Real min_joule_heat = 0.0;
    Real phi_solver_rel_residual = 0.0;
    OpheliePhiSolverKind phi_solver_kind = OpheliePhiSolverKind::GMRES;
    Real max_phi_imag = 0.0;
    Real max_phi_rhs_imag = 0.0;
    bool with_phi_correction = false;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_OBSERVABLES_H
