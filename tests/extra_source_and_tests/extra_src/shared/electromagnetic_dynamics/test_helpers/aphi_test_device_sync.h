#ifndef APHI_TEST_DEVICE_SYNC_H
#define APHI_TEST_DEVICE_SYNC_H

#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "execution_policy.h"
#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Sync one particle variable from device to host before getVariableDataByName (SYCL tests only). */
template <typename DataType>
inline void syncVariableToHost(BaseParticles &particles, const std::string &name)
{
#if SPHINXSYS_USE_SYCL
    particles.template getVariableByName<DataType>(name)->prepareForOutput(execution::par_device);
#else
    (void)particles;
    (void)name;
#endif
}

/** Push host particle data to device before dynamics reads it (SYCL tests only). */
template <typename DataType>
inline void syncVariableToDevice(BaseParticles &particles, const std::string &name)
{
#if SPHINXSYS_USE_SYCL
    particles.template getVariableByName<DataType>(name)->finalizeLoadIn(execution::par_device);
#else
    (void)particles;
    (void)name;
#endif
}

inline void syncAphiBlockToHost(BaseParticles &particles, const AphiBlockNames &block_names)
{
    syncVariableToHost<Vecd>(particles, block_names.a_real);
    syncVariableToHost<Vecd>(particles, block_names.a_imag);
    syncVariableToHost<Real>(particles, block_names.phi_real);
    syncVariableToHost<Real>(particles, block_names.phi_imag);
}

inline void syncAphiBlockToDevice(BaseParticles &particles, const AphiBlockNames &block_names)
{
    syncVariableToDevice<Vecd>(particles, block_names.a_real);
    syncVariableToDevice<Vecd>(particles, block_names.a_imag);
    syncVariableToDevice<Real>(particles, block_names.phi_real);
    syncVariableToDevice<Real>(particles, block_names.phi_imag);
}

inline void syncAphiMaterialToHost(BaseParticles &particles, const AphiMaterialNames &material_names)
{
    syncVariableToHost<Real>(particles, material_names.sigma);
    syncVariableToHost<Real>(particles, material_names.nu);
}

inline void syncAphiLhsRegressionToHost(BaseParticles &particles, const AphiVariableNames &names)
{
    syncAphiBlockToHost(particles, names.lhs);
    syncAphiBlockToHost(particles, names.residual);
    syncAphiBlockToHost(particles, names.rhs);
}

/** Pull Joule-heating / thermal state from device after CK dynamics (SYCL). */
inline void syncThermalCouplingFieldsToHost(BaseParticles &particles, const std::string &temperature_name,
                                            const std::string &rho_cp_name, const std::string &joule_source_name)
{
    syncVariableToHost<Real>(particles, temperature_name);
    syncVariableToHost<Real>(particles, rho_cp_name);
    syncVariableToHost<Real>(particles, joule_source_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
}

/** Push host-initialized thermal inputs to device before CK thermal step (SYCL). Never push Joule without device→host first. */
inline void syncThermalCouplingInitialFieldsToDevice(BaseParticles &particles, const std::string &temperature_name,
                                                     const std::string &rho_cp_name)
{
    syncVariableToDevice<Real>(particles, temperature_name);
    syncVariableToDevice<Real>(particles, rho_cp_name);
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEST_DEVICE_SYNC_H
