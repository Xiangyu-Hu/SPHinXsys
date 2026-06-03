#ifndef ELECTROMAGNETIC_OPHELIE_CURRENT_MOMENT_SOURCE_H
#define ELECTROMAGNETIC_OPHELIE_CURRENT_MOMENT_SOURCE_H

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

/** Unified Biot–Savart source sample: current moment in A·m (real/imag phasor components). */
struct OphelieCurrentMomentSample
{
    Vecd position = Vecd::Zero();
    Vecd moment_real = Vecd::Zero();
    Vecd moment_imag = Vecd::Zero();
};

inline OphelieCurrentMomentSample makeVolumeCurrentMoment(const Vecd &position, const Vecd &j_src_real,
                                                          const Vecd &j_src_imag, Real volume)
{
    OphelieCurrentMomentSample sample;
    sample.position = position;
    sample.moment_real = j_src_real * volume;
    sample.moment_imag = j_src_imag * volume;
    return sample;
}

inline OphelieCurrentMomentSample makeFilamentCurrentMoment(const Vecd &position, const Vecd &dl, Real current_real,
                                                            Real current_imag = 0.0)
{
    OphelieCurrentMomentSample sample;
    sample.position = position;
    sample.moment_real = current_real * dl;
    sample.moment_imag = current_imag * dl;
    return sample;
}

inline void accumulateBiotSavartFromMoments(const StdVec<OphelieCurrentMomentSample> &moments, const Vecd &observer,
                                            Real mu0, Real softening_length, Vecd &a_real, Vecd &a_imag, Vecd &b_real,
                                            Vecd &b_imag)
{
    a_real = Vecd::Zero();
    a_imag = Vecd::Zero();
    b_real = Vecd::Zero();
    b_imag = Vecd::Zero();
    const Real coeff = mu0 / (4.0 * Pi);
    const Real eps2 = softening_length * softening_length;
    for (const OphelieCurrentMomentSample &moment : moments)
    {
        const Vecd r = observer - moment.position;
        const Real r2 = r.squaredNorm() + eps2;
        const Real inv_r = 1.0 / std::sqrt(r2);
        const Real inv_r3 = inv_r / r2;
        a_real += coeff * moment.moment_real * inv_r;
        a_imag += coeff * moment.moment_imag * inv_r;
        b_real += coeff * moment.moment_real.cross(r) * inv_r3;
        b_imag += coeff * moment.moment_imag.cross(r) * inv_r3;
    }
}

inline void applyFilamentBiotSavartToGlassHost(BaseParticles &glass_particles, const OphelieGlassFieldNames &glass_names,
                                               const StdVec<OphelieCurrentMomentSample> &moments, Real mu0,
                                               Real softening_length)
{
    const size_t n = glass_particles.TotalRealParticles();
    syncVariableToHost<Vecd>(glass_particles, "Position");
    syncVariableToHost<Vecd>(glass_particles, glass_names.a_coil_real);
    syncVariableToHost<Vecd>(glass_particles, glass_names.a_coil_imag);
    syncVariableToHost<Vecd>(glass_particles, glass_names.b_coil_real);
    syncVariableToHost<Vecd>(glass_particles, glass_names.b_coil_imag);

    Vecd *a_coil_real = glass_particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = glass_particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Vecd *b_coil_real = glass_particles.getVariableDataByName<Vecd>(glass_names.b_coil_real);
    Vecd *b_coil_imag = glass_particles.getVariableDataByName<Vecd>(glass_names.b_coil_imag);
    const Vecd *pos = glass_particles.getVariableDataByName<Vecd>("Position");

    for (size_t i = 0; i < n; ++i)
    {
        Vecd a_r = Vecd::Zero();
        Vecd a_i = Vecd::Zero();
        Vecd b_r = Vecd::Zero();
        Vecd b_i = Vecd::Zero();
        accumulateBiotSavartFromMoments(moments, pos[i], mu0, softening_length, a_r, a_i, b_r, b_i);
        a_coil_real[i] = a_r;
        a_coil_imag[i] = a_i;
        b_coil_real[i] = b_r;
        b_coil_imag[i] = b_i;
    }
    syncVariableToDevice<Vecd>(glass_particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(glass_particles, glass_names.a_coil_imag);
    syncVariableToDevice<Vecd>(glass_particles, glass_names.b_coil_real);
    syncVariableToDevice<Vecd>(glass_particles, glass_names.b_coil_imag);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CURRENT_MOMENT_SOURCE_H
