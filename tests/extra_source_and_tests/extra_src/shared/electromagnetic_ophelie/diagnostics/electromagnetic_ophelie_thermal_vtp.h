#ifndef ELECTROMAGNETIC_OPHELIE_THERMAL_VTP_H
#define ELECTROMAGNETIC_OPHELIE_THERMAL_VTP_H

#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_progress.h"
#include "io_environment.h"
#include "io_vtk.h"
#include "predefined_bodies.h"
#include "sph_system.h"

#include <iomanip>
#include <sstream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Optional per-step VTP dumps during explicit thermal integration. */
struct OphelieThermalVtpRecordingOptions
{
    bool enabled = false;
    SPHSystem *sph_system = nullptr;
    SolidBody *glass_body = nullptr;
    const OphelieGlassFieldNames *names = nullptr;
    const OphelieParameters *params = nullptr;
    /** Write every N steps when > 0; when 0, write only at end of integration. */
    size_t record_interval = 0;
    size_t iteration_base = 0;
    bool include_diffusion_fields = false;
};

/** Pull device thermal / Q fields to host before BodyStatesRecordingToVtp. */
inline void syncOphelieThermalFieldsToHostForVtp(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                const OphelieParameters &params, bool include_diffusion_fields)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, names.joule_heat);
    if (params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
    }
    if (include_diffusion_fields)
    {
        syncVariableToHost<Real>(particles, kOphelieThermalLaplaceTField);
        syncVariableToHost<Real>(particles, kOphelieThermalConductivityField);
        syncVariableToHost<Real>(particles, kOphelieThermalBoundaryMaskField);
    }
}

inline void registerOphelieThermalBodyStatesRecording(BodyStatesRecordingToVtp &recorder, SolidBody &glass_body,
                                                    const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                                    bool include_diffusion_fields)
{
    recorder.addToWrite<Real>(glass_body, kOphelieTemperatureField);
    recorder.addToWrite<Real>(glass_body, kOphelieThermalDeltaTField);
    recorder.addToWrite<Real>(glass_body, names.joule_heat);
    if (params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        recorder.addToWrite<Real>(glass_body, names.joule_heat_edge_recon_complex);
    }
    if (include_diffusion_fields)
    {
        recorder.addToWrite<Real>(glass_body, kOphelieThermalLaplaceTField);
        recorder.addToWrite<Real>(glass_body, kOphelieThermalConductivityField);
        recorder.addToWrite<Real>(glass_body, kOphelieThermalBoundaryMaskField);
    }
}

inline void writeOphelieThermalBodyStatesVtp(SPHSystem &sph_system, SolidBody &glass_body,
                                           const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                           size_t iteration, bool include_diffusion_fields)
{
    syncOphelieThermalFieldsToHostForVtp(glass_body.getBaseParticles(), names, params, include_diffusion_fields);
    BodyStatesRecordingToVtp recorder(sph_system);
    registerOphelieThermalBodyStatesRecording(recorder, glass_body, names, params, include_diffusion_fields);
    glass_body.setNewlyUpdated();
    recorder.writeToFile(iteration);
    std::ostringstream filename;
    filename << "GlassBody_ite_" << std::setw(10) << std::setfill('0') << iteration << ".vtp";
    logOphelieOutputArtifact(IO::getEnvironment().OutputFolder() + "/" + filename.str());
}

inline void writeOphelieThermalVtpIfDue(const OphelieThermalVtpRecordingOptions *recording, size_t step_index,
                                        size_t n_steps)
{
    if (recording == nullptr || !recording->enabled || recording->sph_system == nullptr ||
        recording->glass_body == nullptr || recording->names == nullptr || recording->params == nullptr)
    {
        return;
    }
    const size_t one_based = step_index + 1;
    const bool at_end = one_based == n_steps;
    const bool on_interval =
        recording->record_interval > 0 && (one_based % recording->record_interval) == 0;
    const bool write_now = (recording->record_interval == 0 && at_end) || on_interval || at_end;
    if (!write_now)
    {
        return;
    }
    const size_t iteration = recording->iteration_base + one_based;
    writeOphelieThermalBodyStatesVtp(*recording->sph_system, *recording->glass_body, *recording->names,
                                     *recording->params, iteration, recording->include_diffusion_fields);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_THERMAL_VTP_H
