#ifndef APHI_GMRES_RESIDUAL_DECOMPOSITION_HELPERS_H
#define APHI_GMRES_RESIDUAL_DECOMPOSITION_HELPERS_H

#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{

struct AphiTrueResidualBlockDecomposition
{
    Real res_a_real_norm = 0.0;
    Real res_a_imag_norm = 0.0;
    Real res_phi_real_norm = 0.0;
    Real res_phi_imag_norm = 0.0;
    Real res_a_total_norm = 0.0;
    Real res_phi_total_norm = 0.0;
    Real res_a_fraction = 0.0;
    Real res_phi_fraction = 0.0;
};

namespace test
{

inline Real hostComponentVolWeightedL2(BaseParticles &particles, const std::string &field_name,
                                       size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, field_name);
    const Real *values = particles.getVariableDataByName<Real>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        squared += vol[i] * values[i] * values[i];
    }
    return std::sqrt(squared);
}

inline Real hostVectorComponentVolWeightedL2(BaseParticles &particles, const std::string &field_name,
                                             size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, field_name);
    const Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        squared += vol[i] * values[i].squaredNorm();
    }
    return std::sqrt(squared);
}

inline AphiTrueResidualBlockDecomposition hostTrueResidualBlockDecomposition(
    BaseParticles &particles, const AphiBlockNames &true_residual_block, size_t total_real_particles)
{
    AphiTrueResidualBlockDecomposition decomposition;
    decomposition.res_a_real_norm =
        hostVectorComponentVolWeightedL2(particles, true_residual_block.a_real, total_real_particles);
    decomposition.res_a_imag_norm =
        hostVectorComponentVolWeightedL2(particles, true_residual_block.a_imag, total_real_particles);
    decomposition.res_phi_real_norm =
        hostComponentVolWeightedL2(particles, true_residual_block.phi_real, total_real_particles);
    decomposition.res_phi_imag_norm =
        hostComponentVolWeightedL2(particles, true_residual_block.phi_imag, total_real_particles);
    decomposition.res_a_total_norm = std::sqrt(decomposition.res_a_real_norm * decomposition.res_a_real_norm +
                                               decomposition.res_a_imag_norm * decomposition.res_a_imag_norm);
    decomposition.res_phi_total_norm = std::sqrt(decomposition.res_phi_real_norm * decomposition.res_phi_real_norm +
                                                 decomposition.res_phi_imag_norm * decomposition.res_phi_imag_norm);
    const Real total_norm = std::sqrt(decomposition.res_a_total_norm * decomposition.res_a_total_norm +
                                      decomposition.res_phi_total_norm * decomposition.res_phi_total_norm);
    if (total_norm > TinyReal)
    {
        decomposition.res_a_fraction = decomposition.res_a_total_norm / total_norm;
        decomposition.res_phi_fraction = decomposition.res_phi_total_norm / total_norm;
    }
    return decomposition;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_RESIDUAL_DECOMPOSITION_HELPERS_H
