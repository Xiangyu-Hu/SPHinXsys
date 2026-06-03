#ifndef ELECTROMAGNETIC_OPHELIE_DEVICE_SYNC_H
#define ELECTROMAGNETIC_OPHELIE_DEVICE_SYNC_H

#include "electromagnetic_ophelie_field_names.h"
#include "execution_policy.h"
#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

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

inline void syncCoilSourceFieldsToDevice(BaseParticles &particles, const OphelieCoilFieldNames &names)
{
    syncVariableToDevice<Vecd>(particles, names.j_src_real);
    syncVariableToDevice<Vecd>(particles, names.j_src_imag);
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Real>(particles, "VolumetricMeasure");
}

inline void syncGlassElectromagneticFieldsToDevice(BaseParticles &particles, const OphelieGlassFieldNames &names)
{
    syncVariableToDevice<Real>(particles, names.sigma);
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_ind_real);
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
    syncVariableToDevice<Vecd>(particles, names.a_src_imag);
    syncVariableToDevice<Vecd>(particles, names.b_src_real);
    syncVariableToDevice<Vecd>(particles, names.b_src_imag);
    syncVariableToDevice<Vecd>(particles, names.j_imag);
    syncVariableToDevice<Vecd>(particles, "Position");
    syncVariableToDevice<Real>(particles, "VolumetricMeasure");
}

inline void syncGlassElectromagneticFieldsToHost(BaseParticles &particles, const OphelieGlassFieldNames &names)
{
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.a_ind_real);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Vecd>(particles, names.a_src_imag);
    syncVariableToHost<Vecd>(particles, names.b_src_real);
    syncVariableToHost<Vecd>(particles, names.b_src_imag);
    syncVariableToHost<Vecd>(particles, names.e_real);
    syncVariableToHost<Vecd>(particles, names.e_imag);
    syncVariableToHost<Vecd>(particles, names.j_real);
    syncVariableToHost<Vecd>(particles, names.j_imag);
    syncVariableToHost<Real>(particles, names.joule_heat);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DEVICE_SYNC_H
