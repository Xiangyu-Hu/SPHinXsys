#ifndef APHI_CONTACT_LEFT_FIELD_ERROR_HELPERS_H
#define APHI_CONTACT_LEFT_FIELD_ERROR_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline bool isLeftHalfCoreParticle(const Vecd &position, Real body_length, Real body_height, Real body_width,
                                   Real core_shell, Real x_interface)
{
    return isCoreParticle(position, body_length, body_height, body_width, core_shell) &&
           position[0] < x_interface - TinyReal;
}

inline bool isLeftBodyCoreInteriorParticle(const Vecd &position, Real body_length, Real body_height, Real body_width,
                                           Real core_shell, Real x_interface, Real interface_band_half_width)
{
    return isLeftHalfCoreParticle(position, body_length, body_height, body_width, core_shell, x_interface) &&
           isAwayFromInterface(position, x_interface, interface_band_half_width);
}

inline bool isLeftBodyInterfaceBandParticle(const Vecd &position, Real body_length, Real body_height, Real body_width,
                                            Real core_shell, Real x_interface, Real interface_band_half_width)
{
    return isCoreParticle(position, body_length, body_height, body_width, core_shell) &&
           position[0] < x_interface - TinyReal &&
           std::abs(position[0] - x_interface) <= interface_band_half_width + TinyReal;
}

enum class AphiLeftBodyErrorZone
{
    FullLeftHalfCore,
    CoreInterior,
    InterfaceBand
};

inline bool isLeftBodyErrorZoneParticle(AphiLeftBodyErrorZone zone, const Vecd &position, Real body_length,
                                        Real body_height, Real body_width, Real core_shell, Real x_interface,
                                        Real interface_band_half_width)
{
    switch (zone)
    {
    case AphiLeftBodyErrorZone::FullLeftHalfCore:
        return isLeftHalfCoreParticle(position, body_length, body_height, body_width, core_shell, x_interface);
    case AphiLeftBodyErrorZone::CoreInterior:
        return isLeftBodyCoreInteriorParticle(position, body_length, body_height, body_width, core_shell, x_interface,
                                              interface_band_half_width);
    case AphiLeftBodyErrorZone::InterfaceBand:
        return isLeftBodyInterfaceBandParticle(position, body_length, body_height, body_width, core_shell, x_interface,
                                               interface_band_half_width);
    }
    return false;
}

struct AphiLeftBodyGaugeErrorMetrics
{
    Real block_L2_rel = 0.0;
    Real block_Linf_rel = 0.0;
    Real A_real_mean_error = 0.0;
    Real A_imag_mean_error = 0.0;
    Real phi_real_mean_error = 0.0;
    Real phi_imag_mean_error = 0.0;
    Real block_L2_after_mean_subtraction_rel = 0.0;
    Real block_Linf_after_mean_subtraction_rel = 0.0;
    size_t particle_count = 0;
};

struct AphiLeftBodyPhiComponentErrorMetrics
{
    Real phi_L2_rel = 0.0;
    Real phi_Linf_rel = 0.0;
    Real phi_real_mean_error = 0.0;
    Real phi_imag_mean_error = 0.0;
    Real phi_L2_after_mean_subtraction_rel = 0.0;
    Real phi_Linf_after_mean_subtraction_rel = 0.0;
    size_t particle_count = 0;
};

struct AphiLeftBodySplitComparisonMetrics
{
    Real contact_vs_exact_L2_rel = 0.0;
    Real contact_vs_exact_Linf_rel = 0.0;
    Real contact_vs_mono_L2_rel = 0.0;
    Real contact_vs_mono_Linf_rel = 0.0;
    Real mono_vs_exact_L2_rel = 0.0;
    Real mono_vs_exact_Linf_rel = 0.0;
    size_t matched_mono_particles = 0;
    size_t missing_mono_particles = 0;
};

struct AphiLeftBodyZonePartitionMetrics
{
    AphiLeftBodyGaugeErrorMetrics full_left_half{};
    AphiLeftBodyGaugeErrorMetrics core_interior{};
    AphiLeftBodyGaugeErrorMetrics interface_band{};
    AphiLeftBodySplitComparisonMetrics core_interior_split{};
    AphiLeftBodySplitComparisonMetrics interface_band_split{};
};

inline void collectLeftBodyZoneBlockByPosition(AphiBlockMapByPosition &block_map, BaseParticles &particles,
                                               const AphiBlockNames &block, AphiLeftBodyErrorZone zone,
                                               Real body_length, Real body_height, Real body_width, Real core_shell,
                                               Real x_interface, Real interface_band_half_width)
{
    syncAphiBlockToHost(particles, block);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block.phi_imag);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isLeftBodyErrorZoneParticle(zone, positions[i], body_length, body_height, body_width, core_shell,
                                         x_interface, interface_band_half_width))
        {
            continue;
        }

        AphiContactBlockByPosition value;
        value.a_real = a_real[i];
        value.a_imag = a_imag[i];
        value.phi_real = phi_real[i];
        value.phi_imag = phi_imag[i];
        block_map[makeContactPositionKey(positions[i])] = value;
    }
}

inline void collectLeftHalfBlockByPosition(AphiBlockMapByPosition &block_map, BaseParticles &particles,
                                           const AphiBlockNames &block, Real body_length, Real body_height,
                                           Real body_width, Real core_shell, Real x_interface)
{
    collectLeftBodyZoneBlockByPosition(block_map, particles, block, AphiLeftBodyErrorZone::FullLeftHalfCore,
                                       body_length, body_height, body_width, core_shell, x_interface, 0.0);
}

inline AphiLeftBodyGaugeErrorMetrics hostLeftBodyGaugeErrorMetrics(
    BaseParticles &particles, const AphiBlockNames &approx_block, const AphiBlockNames &exact_block,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real x_interface, AphiLeftBodyErrorZone zone = AphiLeftBodyErrorZone::FullLeftHalfCore,
    Real interface_band_half_width = 0.0)
{
    syncAphiBlockToHost(particles, approx_block);
    syncAphiBlockToHost(particles, exact_block);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Vecd *approx_a_real = particles.getVariableDataByName<Vecd>(approx_block.a_real);
    const Vecd *approx_a_imag = particles.getVariableDataByName<Vecd>(approx_block.a_imag);
    const Real *approx_phi_real = particles.getVariableDataByName<Real>(approx_block.phi_real);
    const Real *approx_phi_imag = particles.getVariableDataByName<Real>(approx_block.phi_imag);
    const Vecd *exact_a_real = particles.getVariableDataByName<Vecd>(exact_block.a_real);
    const Vecd *exact_a_imag = particles.getVariableDataByName<Vecd>(exact_block.a_imag);
    const Real *exact_phi_real = particles.getVariableDataByName<Real>(exact_block.phi_real);
    const Real *exact_phi_imag = particles.getVariableDataByName<Real>(exact_block.phi_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    AphiLeftBodyGaugeErrorMetrics metrics;
    Real exact_norm_squared = 0.0;
    Real error_norm_squared = 0.0;
    Real mean_subtracted_error_squared = 0.0;
    Real total_volume = 0.0;

    Real sum_vol_a_real_x = 0.0;
    Real sum_vol_a_real_y = 0.0;
    Real sum_vol_a_real_z = 0.0;
    Real sum_vol_a_imag_x = 0.0;
    Real sum_vol_a_imag_y = 0.0;
    Real sum_vol_a_imag_z = 0.0;
    Real sum_vol_phi_real = 0.0;
    Real sum_vol_phi_imag = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isLeftBodyErrorZoneParticle(zone, positions[i], body_length, body_height, body_width, core_shell,
                                         x_interface, interface_band_half_width))
        {
            continue;
        }

        const Vecd diff_a_real = approx_a_real[i] - exact_a_real[i];
        const Vecd diff_a_imag = approx_a_imag[i] - exact_a_imag[i];
        const Real diff_phi_real = approx_phi_real[i] - exact_phi_real[i];
        const Real diff_phi_imag = approx_phi_imag[i] - exact_phi_imag[i];
        const Real vol_i = vol[i];

        metrics.particle_count += 1;
        total_volume += vol_i;
        exact_norm_squared += vol_i * (exact_a_real[i].squaredNorm() + exact_a_imag[i].squaredNorm() +
                                       exact_phi_real[i] * exact_phi_real[i] + exact_phi_imag[i] * exact_phi_imag[i]);
        error_norm_squared += vol_i * (diff_a_real.squaredNorm() + diff_a_imag.squaredNorm() + diff_phi_real * diff_phi_real +
                                       diff_phi_imag * diff_phi_imag);

        sum_vol_a_real_x += vol_i * diff_a_real[0];
        sum_vol_a_real_y += vol_i * diff_a_real[1];
        sum_vol_a_real_z += vol_i * diff_a_real[2];
        sum_vol_a_imag_x += vol_i * diff_a_imag[0];
        sum_vol_a_imag_y += vol_i * diff_a_imag[1];
        sum_vol_a_imag_z += vol_i * diff_a_imag[2];
        sum_vol_phi_real += vol_i * diff_phi_real;
        sum_vol_phi_imag += vol_i * diff_phi_imag;

        metrics.block_Linf_rel = std::max(metrics.block_Linf_rel, diff_a_real.norm());
        metrics.block_Linf_rel = std::max(metrics.block_Linf_rel, diff_a_imag.norm());
        metrics.block_Linf_rel = std::max(metrics.block_Linf_rel, std::abs(diff_phi_real));
        metrics.block_Linf_rel = std::max(metrics.block_Linf_rel, std::abs(diff_phi_imag));
    }

    if (metrics.particle_count == 0 || total_volume <= TinyReal)
    {
        return metrics;
    }

    const Real inv_vol = Real(1) / total_volume;
    const Vecd mean_diff_a_real(sum_vol_a_real_x * inv_vol, sum_vol_a_real_y * inv_vol, sum_vol_a_real_z * inv_vol);
    const Vecd mean_diff_a_imag(sum_vol_a_imag_x * inv_vol, sum_vol_a_imag_y * inv_vol, sum_vol_a_imag_z * inv_vol);
    const Real mean_diff_phi_real = sum_vol_phi_real * inv_vol;
    const Real mean_diff_phi_imag = sum_vol_phi_imag * inv_vol;

    metrics.A_real_mean_error = mean_diff_a_real.norm();
    metrics.A_imag_mean_error = mean_diff_a_imag.norm();
    metrics.phi_real_mean_error = std::abs(mean_diff_phi_real);
    metrics.phi_imag_mean_error = std::abs(mean_diff_phi_imag);

    const Real exact_norm = std::sqrt(exact_norm_squared);
    metrics.block_L2_rel = std::sqrt(error_norm_squared) / (exact_norm + TinyReal);
    metrics.block_Linf_rel = metrics.block_Linf_rel / (exact_norm + TinyReal);

    Real mean_subtracted_linf = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isLeftBodyErrorZoneParticle(zone, positions[i], body_length, body_height, body_width, core_shell,
                                         x_interface, interface_band_half_width))
        {
            continue;
        }

        const Vecd centered_a_real = approx_a_real[i] - exact_a_real[i] - mean_diff_a_real;
        const Vecd centered_a_imag = approx_a_imag[i] - exact_a_imag[i] - mean_diff_a_imag;
        const Real centered_phi_real = approx_phi_real[i] - exact_phi_real[i] - mean_diff_phi_real;
        const Real centered_phi_imag = approx_phi_imag[i] - exact_phi_imag[i] - mean_diff_phi_imag;
        const Real vol_i = vol[i];

        mean_subtracted_error_squared += vol_i * (centered_a_real.squaredNorm() + centered_a_imag.squaredNorm() +
                                                centered_phi_real * centered_phi_real +
                                                centered_phi_imag * centered_phi_imag);
        mean_subtracted_linf = std::max(mean_subtracted_linf, centered_a_real.norm());
        mean_subtracted_linf = std::max(mean_subtracted_linf, centered_a_imag.norm());
        mean_subtracted_linf = std::max(mean_subtracted_linf, std::abs(centered_phi_real));
        mean_subtracted_linf = std::max(mean_subtracted_linf, std::abs(centered_phi_imag));
    }

    metrics.block_L2_after_mean_subtraction_rel = std::sqrt(mean_subtracted_error_squared) / (exact_norm + TinyReal);
    metrics.block_Linf_after_mean_subtraction_rel = mean_subtracted_linf / (exact_norm + TinyReal);
    return metrics;
}

inline AphiLeftBodyPhiComponentErrorMetrics hostLeftBodyPhiComponentErrorMetrics(
    BaseParticles &particles, const AphiBlockNames &approx_block, const AphiBlockNames &exact_block,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real x_interface, AphiLeftBodyErrorZone zone = AphiLeftBodyErrorZone::FullLeftHalfCore,
    Real interface_band_half_width = 0.0)
{
    syncAphiBlockToHost(particles, approx_block);
    syncAphiBlockToHost(particles, exact_block);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Real *approx_phi_real = particles.getVariableDataByName<Real>(approx_block.phi_real);
    const Real *approx_phi_imag = particles.getVariableDataByName<Real>(approx_block.phi_imag);
    const Real *exact_phi_real = particles.getVariableDataByName<Real>(exact_block.phi_real);
    const Real *exact_phi_imag = particles.getVariableDataByName<Real>(exact_block.phi_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    AphiLeftBodyPhiComponentErrorMetrics metrics;
    Real exact_norm_squared = 0.0;
    Real error_norm_squared = 0.0;
    Real mean_subtracted_error_squared = 0.0;
    Real total_volume = 0.0;
    Real sum_vol_phi_real = 0.0;
    Real sum_vol_phi_imag = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isLeftBodyErrorZoneParticle(zone, positions[i], body_length, body_height, body_width, core_shell,
                                         x_interface, interface_band_half_width))
        {
            continue;
        }

        const Real diff_phi_real = approx_phi_real[i] - exact_phi_real[i];
        const Real diff_phi_imag = approx_phi_imag[i] - exact_phi_imag[i];
        const Real vol_i = vol[i];

        metrics.particle_count += 1;
        total_volume += vol_i;
        exact_norm_squared +=
            vol_i * (exact_phi_real[i] * exact_phi_real[i] + exact_phi_imag[i] * exact_phi_imag[i]);
        error_norm_squared += vol_i * (diff_phi_real * diff_phi_real + diff_phi_imag * diff_phi_imag);
        sum_vol_phi_real += vol_i * diff_phi_real;
        sum_vol_phi_imag += vol_i * diff_phi_imag;

        metrics.phi_Linf_rel = std::max(metrics.phi_Linf_rel, std::abs(diff_phi_real));
        metrics.phi_Linf_rel = std::max(metrics.phi_Linf_rel, std::abs(diff_phi_imag));
    }

    if (metrics.particle_count == 0 || total_volume <= TinyReal)
    {
        return metrics;
    }

    const Real inv_vol = Real(1) / total_volume;
    const Real mean_diff_phi_real = sum_vol_phi_real * inv_vol;
    const Real mean_diff_phi_imag = sum_vol_phi_imag * inv_vol;
    metrics.phi_real_mean_error = std::abs(mean_diff_phi_real);
    metrics.phi_imag_mean_error = std::abs(mean_diff_phi_imag);

    const Real exact_norm = std::sqrt(exact_norm_squared);
    metrics.phi_L2_rel = std::sqrt(error_norm_squared) / (exact_norm + TinyReal);
    metrics.phi_Linf_rel = metrics.phi_Linf_rel / (exact_norm + TinyReal);

    Real mean_subtracted_linf = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isLeftBodyErrorZoneParticle(zone, positions[i], body_length, body_height, body_width, core_shell,
                                         x_interface, interface_band_half_width))
        {
            continue;
        }

        const Real centered_phi_real = approx_phi_real[i] - exact_phi_real[i] - mean_diff_phi_real;
        const Real centered_phi_imag = approx_phi_imag[i] - exact_phi_imag[i] - mean_diff_phi_imag;
        const Real vol_i = vol[i];

        mean_subtracted_error_squared +=
            vol_i * (centered_phi_real * centered_phi_real + centered_phi_imag * centered_phi_imag);
        mean_subtracted_linf = std::max(mean_subtracted_linf, std::abs(centered_phi_real));
        mean_subtracted_linf = std::max(mean_subtracted_linf, std::abs(centered_phi_imag));
    }

    metrics.phi_L2_after_mean_subtraction_rel = std::sqrt(mean_subtracted_error_squared) / (exact_norm + TinyReal);
    metrics.phi_Linf_after_mean_subtraction_rel = mean_subtracted_linf / (exact_norm + TinyReal);
    return metrics;
}

inline Real relativeBlockMapL2Error(const AphiBlockMapByPosition &approx, const AphiBlockMapByPosition &reference,
                                    size_t &matched_particles, size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real error_squared = 0.0;
    Real reference_squared = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = approx.find(entry.first);
        if (it == approx.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        const AphiContactBlockByPosition &ref = entry.second;
        const AphiContactBlockByPosition &app = it->second;
        const Vecd diff_a_real = app.a_real - ref.a_real;
        const Vecd diff_a_imag = app.a_imag - ref.a_imag;
        const Real diff_phi_real = app.phi_real - ref.phi_real;
        const Real diff_phi_imag = app.phi_imag - ref.phi_imag;

        error_squared += diff_a_real.squaredNorm() + diff_a_imag.squaredNorm() + diff_phi_real * diff_phi_real +
                         diff_phi_imag * diff_phi_imag;
        reference_squared += ref.a_real.squaredNorm() + ref.a_imag.squaredNorm() + ref.phi_real * ref.phi_real +
                             ref.phi_imag * ref.phi_imag;
    }

    return std::sqrt(error_squared) / (std::sqrt(reference_squared) + TinyReal);
}

inline Real relativeBlockMapLinfError(const AphiBlockMapByPosition &approx, const AphiBlockMapByPosition &reference,
                                      size_t &matched_particles, size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real max_diff = 0.0;
    Real reference_scale = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = approx.find(entry.first);
        if (it == approx.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        const AphiContactBlockByPosition &ref = entry.second;
        const AphiContactBlockByPosition &app = it->second;
        const Vecd diff_a_real = app.a_real - ref.a_real;
        const Vecd diff_a_imag = app.a_imag - ref.a_imag;
        const Real diff_phi_real = app.phi_real - ref.phi_real;
        const Real diff_phi_imag = app.phi_imag - ref.phi_imag;

        max_diff = std::max(max_diff, diff_a_real.norm());
        max_diff = std::max(max_diff, diff_a_imag.norm());
        max_diff = std::max(max_diff, std::abs(diff_phi_real));
        max_diff = std::max(max_diff, std::abs(diff_phi_imag));
        reference_scale = std::max(reference_scale, ref.a_real.norm());
        reference_scale = std::max(reference_scale, ref.a_imag.norm());
        reference_scale = std::max(reference_scale, std::abs(ref.phi_real));
        reference_scale = std::max(reference_scale, std::abs(ref.phi_imag));
    }

    return max_diff / (reference_scale + TinyReal);
}

inline AphiLeftBodySplitComparisonMetrics hostLeftBodySplitComparisonMetrics(
    BaseParticles &contact_left_particles, BaseParticles &mono_particles, const AphiBlockNames &contact_solution,
    const AphiBlockNames &mono_solution, const AphiBlockNames &exact_block, Real body_length, Real body_height,
    Real body_width, Real core_shell, Real x_interface, AphiLeftBodyErrorZone zone = AphiLeftBodyErrorZone::FullLeftHalfCore,
    Real interface_band_half_width = 0.0)
{
    AphiBlockMapByPosition contact_map;
    AphiBlockMapByPosition mono_map;
    AphiBlockMapByPosition exact_map;

    collectLeftBodyZoneBlockByPosition(contact_map, contact_left_particles, contact_solution, zone, body_length,
                                       body_height, body_width, core_shell, x_interface, interface_band_half_width);
    collectLeftBodyZoneBlockByPosition(mono_map, mono_particles, mono_solution, zone, body_length, body_height,
                                       body_width, core_shell, x_interface, interface_band_half_width);
    collectLeftBodyZoneBlockByPosition(exact_map, contact_left_particles, exact_block, zone, body_length, body_height,
                                       body_width, core_shell, x_interface, interface_band_half_width);

    AphiLeftBodySplitComparisonMetrics metrics;
    size_t matched = 0;
    size_t missing = 0;
    metrics.contact_vs_exact_Linf_rel = relativeBlockMapLinfError(contact_map, exact_map, matched, missing);
    metrics.contact_vs_exact_L2_rel = relativeBlockMapL2Error(contact_map, exact_map, matched, missing);
    metrics.contact_vs_mono_Linf_rel = relativeBlockMapLinfError(contact_map, mono_map, matched, missing);
    metrics.contact_vs_mono_L2_rel = relativeBlockMapL2Error(contact_map, mono_map, matched, missing);
    metrics.mono_vs_exact_Linf_rel = relativeBlockMapLinfError(mono_map, exact_map, matched, missing);
    metrics.mono_vs_exact_L2_rel = relativeBlockMapL2Error(mono_map, exact_map, matched, missing);
    metrics.matched_mono_particles = matched;
    metrics.missing_mono_particles = missing;
    return metrics;
}

inline AphiLeftBodyZonePartitionMetrics hostLeftBodyZonePartitionMetrics(
    BaseParticles &contact_left_particles, BaseParticles &mono_particles, const AphiBlockNames &contact_solution,
    const AphiBlockNames &mono_solution, const AphiBlockNames &exact_block, const Vecd *left_positions,
    size_t left_total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell,
    Real x_interface, Real interface_band_half_width)
{
    AphiLeftBodyZonePartitionMetrics metrics;
    metrics.full_left_half = hostLeftBodyGaugeErrorMetrics(
        contact_left_particles, contact_solution, exact_block, left_positions, left_total_real_particles, body_length,
        body_height, body_width, core_shell, x_interface, AphiLeftBodyErrorZone::FullLeftHalfCore,
        interface_band_half_width);
    metrics.core_interior = hostLeftBodyGaugeErrorMetrics(
        contact_left_particles, contact_solution, exact_block, left_positions, left_total_real_particles, body_length,
        body_height, body_width, core_shell, x_interface, AphiLeftBodyErrorZone::CoreInterior,
        interface_band_half_width);
    metrics.interface_band = hostLeftBodyGaugeErrorMetrics(
        contact_left_particles, contact_solution, exact_block, left_positions, left_total_real_particles, body_length,
        body_height, body_width, core_shell, x_interface, AphiLeftBodyErrorZone::InterfaceBand,
        interface_band_half_width);
    metrics.core_interior_split = hostLeftBodySplitComparisonMetrics(
        contact_left_particles, mono_particles, contact_solution, mono_solution, exact_block, body_length, body_height,
        body_width, core_shell, x_interface, AphiLeftBodyErrorZone::CoreInterior, interface_band_half_width);
    metrics.interface_band_split = hostLeftBodySplitComparisonMetrics(
        contact_left_particles, mono_particles, contact_solution, mono_solution, exact_block, body_length, body_height,
        body_width, core_shell, x_interface, AphiLeftBodyErrorZone::InterfaceBand, interface_band_half_width);
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_LEFT_FIELD_ERROR_HELPERS_H
