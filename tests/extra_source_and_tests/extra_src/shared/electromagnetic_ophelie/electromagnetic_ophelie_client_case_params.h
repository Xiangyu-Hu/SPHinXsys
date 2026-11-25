/**
 * @file electromagnetic_ophelie_client_case_params.h
 * @brief Client NC / stirring case presets (Doc229).
 *
 * Power ledger: generator_input_power_W=30000 is NOT absorbed_power_W.
 * Current A / peak-vs-RMS still UNKNOWN — do not invent 1000 A as confirmed.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_CLIENT_CASE_PARAMS_H
#define ELECTROMAGNETIC_OPHELIE_CLIENT_CASE_PARAMS_H

#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_material_laws.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"

#include <fstream>
#include <iomanip>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OphelieClientConfirmation
{
    Confirmed = 0,
    Provisional = 1,
    Borrowed = 2,
    Unknown = 3
};

inline const char *clientConfirmationName(OphelieClientConfirmation c)
{
    switch (c)
    {
    case OphelieClientConfirmation::Confirmed:
        return "confirmed";
    case OphelieClientConfirmation::Provisional:
        return "provisional";
    case OphelieClientConfirmation::Borrowed:
        return "borrowed";
    default:
        return "unknown";
    }
}

enum class OphelieClientSigmaPreset
{
    ClientLnExpProvisional = 0,
    FrenchTable1Sigma16Borrowed = 1
};

enum class OphelieClientMuPreset
{
    /** mu(T) formula table installed; hydro may still freeze later. */
    ClientFormula = 0,
    /** Constant mu = formula(T_ref); mark as frozen development path. */
    FrozenAtReferenceT = 1
};

struct OphelieClientPowerLedger
{
    Real generator_input_power_w = Real(30000);
    OphelieClientConfirmation generator_power_status = OphelieClientConfirmation::Confirmed;
    Real coil_current_a = Real(0); // unset
    const char *current_convention = "unknown_peak_or_rms";
    OphelieClientConfirmation current_status = OphelieClientConfirmation::Unknown;
    Real p_em_absorbed_raw_w = Real(0); // filled after EM solve
    Real absorbed_power_target_w = Real(0);
    bool absorbed_power_target_active = false;
    Real coupling_efficiency = Real(0); // estimate only if explicitly set
    bool coupling_efficiency_active = false;
};

struct OphelieClientCaseParams
{
    OphelieClientAssemblyKind assembly = OphelieClientAssemblyKind::NaturalConvection;
    const char *case_id = "";
    const char *schema_version = "client-case-2026-11-20";

    Real frequency_hz = Real(100000);
    OphelieClientConfirmation frequency_status = OphelieClientConfirmation::Provisional;

    Real stirrer_rpm = Real(0);
    OphelieClientConfirmation rpm_status = OphelieClientConfirmation::Confirmed;

    Real t_initial_k = Real(1273.15);
    OphelieClientConfirmation t_initial_status = OphelieClientConfirmation::Provisional;
    const char *t_initial_note = "development IC; must NOT equal probe target as forced BC";

    Real probe_target_k = Real(1273.15);
    OphelieClientConfirmation probe_target_status = OphelieClientConfirmation::Provisional;
    const char *probe_target_note = "thermocouple magnitude target; missing time window / error bar";

    OphelieClientPowerLedger power{};

    OphelieClientSigmaPreset sigma_preset = OphelieClientSigmaPreset::FrenchTable1Sigma16Borrowed;
    OphelieClientMuPreset mu_preset = OphelieClientMuPreset::FrozenAtReferenceT;
    Real mu_freeze_t_k = kClientNcProbeTargetK;

    Real rho_kg_m3 = Real(2700);
    OphelieClientConfirmation rho_status = OphelieClientConfirmation::Provisional;
    const char *rho_note = "interpreted 2.7 g/cm3; unit pending client confirm";

    Real cp = Real(1150);
    Real k_th = Real(4);
    Real beta = Real(1.0e-5);
    OphelieClientConfirmation cp_k_beta_status = OphelieClientConfirmation::Borrowed;

    bool stirring_enabled = false;
    bool elastic_paddle_required = false;
    const char *mechanics_mode = "none"; // none | rigid_prescribed | elastic_driven

    const char *material_id = kClientMaterialIdFrenchSigma16Borrowed;

    /** Parametric melt fill height [m]; provisional until client confirms liquid depth. */
    Real melt_height_m = Real(0.180);
    OphelieClientConfirmation melt_height_status = OphelieClientConfirmation::Provisional;
    Real particle_spacing_m = Real(0.005); // match client case.json default; override via --dp=
};

inline OphelieClientCaseParams makeClientNaturalConvectionDevelopmentParams()
{
    OphelieClientCaseParams p;
    p.assembly = OphelieClientAssemblyKind::NaturalConvection;
    p.case_id = "client_nc_dev_frozen_mu_sigma16";
    p.frequency_hz = Real(100000);
    p.stirrer_rpm = Real(0);
    p.stirring_enabled = false;
    p.elastic_paddle_required = false;
    p.mechanics_mode = "none";
    p.t_initial_k = kClientNcProbeTargetK;
    p.probe_target_k = kClientNcProbeTargetK;
    p.sigma_preset = OphelieClientSigmaPreset::FrenchTable1Sigma16Borrowed;
    p.mu_preset = OphelieClientMuPreset::FrozenAtReferenceT;
    p.mu_freeze_t_k = kClientNcProbeTargetK;
    p.material_id = kClientMaterialIdFrenchSigma16Borrowed;
    p.power.generator_input_power_w = Real(30000);
    p.power.absorbed_power_target_active = false;
    return p;
}

/** Same geometry / power / freq; client ln-exp sigma for side-by-side EM Q compare. */
inline OphelieClientCaseParams makeClientNaturalConvectionClientSigmaContrastParams()
{
    OphelieClientCaseParams p = makeClientNaturalConvectionDevelopmentParams();
    p.case_id = "client_nc_dev_frozen_mu_client_ln_sigma";
    p.sigma_preset = OphelieClientSigmaPreset::ClientLnExpProvisional;
    p.material_id = kClientMaterialIdMuSigmaProvisional;
    return p;
}

inline OphelieClientCaseParams makeClientStirringDevelopmentParams()
{
    OphelieClientCaseParams p;
    p.assembly = OphelieClientAssemblyKind::MechanicallyStirred;
    p.case_id = "client_stir_dev_10rpm";
    p.frequency_hz = Real(100000);
    p.stirrer_rpm = Real(10);
    p.rpm_status = OphelieClientConfirmation::Confirmed;
    p.stirring_enabled = true;
    p.elastic_paddle_required = true;
    p.mechanics_mode = "elastic_driven"; // Doc228; not yet implemented as PASS
    p.t_initial_k = kClientNcProbeTargetK;
    p.probe_target_k = kClientStirPeakTrendK;
    p.probe_target_status = OphelieClientConfirmation::Provisional;
    p.probe_target_note = "qualitative 1400C trend only; not a hard acceptance gate";
    p.sigma_preset = OphelieClientSigmaPreset::FrenchTable1Sigma16Borrowed;
    p.mu_preset = OphelieClientMuPreset::FrozenAtReferenceT;
    p.mu_freeze_t_k = kClientNcProbeTargetK;
    p.material_id = kClientMaterialIdFrenchSigma16Borrowed;
    p.power.generator_input_power_w = Real(30000);
    return p;
}

/**
 * Map client NC parametric melt into French-reduced geometry knobs so existing
 * frozen_q / natural hydro paths can run on the client Reload without a fork.
 */
inline void applyClientNcMeltToFrenchReduced(OphelieFrenchReducedCaseParams &french,
                                            const OphelieClientMeltGeometryM &melt, Real dp = Real(0.010))
{
    french.glass_radius = melt.melt_radius_m;
    french.glass_half_height = melt.glass_half_height;
    french.glass_center = melt.glass_center;
    french.dp = dp;
    french.require_canonical_bottom_at_zero = false; // bottom at inner_bottom_z, not z=0
    french.coil.loop_radius = melt.coil_loop_radius_m;
    french.coil.num_loops = static_cast<size_t>(melt.coil_num_loops);
    french.coil.z_min = melt.coil_z_min_m;
    french.coil.z_max = melt.coil_z_max_m;
    french.coil_z_user_set = true;
    french.auto_coil_z = false;
    french.frequency_hz = Real(100000); // Doc229 provisional
    french.coil.use_cell_centered_loops = false;
    french.ampere_turns = french.coil.current_per_loop * static_cast<Real>(french.coil.num_loops);
}

inline OphelieFrenchGlassMaterialLaws buildClientCaseMaterialLaws(const OphelieClientCaseParams &p)
{
    OphelieFrenchGlassMaterialLaws laws =
        (p.sigma_preset == OphelieClientSigmaPreset::ClientLnExpProvisional)
            ? makeClientProvisionalMuSigmaMaterialLaws()
            : makeClientGeometryFrenchSigma16BorrowedMaterialLaws();

    // Keep density / borrowed thermal constants aligned with preset numbers.
    laws.rho.clear();
    laws.rho.setPoint(0, p.t_initial_k, p.rho_kg_m3);
    laws.cp.clear();
    laws.cp.setPoint(0, p.t_initial_k, p.cp);
    laws.k.clear();
    laws.k.setPoint(0, p.t_initial_k, p.k_th);
    laws.beta.clear();
    laws.beta.setPoint(0, p.t_initial_k, p.beta);

    if (p.mu_preset == OphelieClientMuPreset::FrozenAtReferenceT)
        freezeClientDynamicViscosityAtTemperature(laws, p.mu_freeze_t_k);
    return laws;
}

/** Install client development materials into literature knobs (does not touch French makers). */
inline void applyClientNcDevelopmentToLiterature(FrenchLiteratureParameters &lit, const OphelieClientCaseParams &p)
{
    lit.regime = FrenchLiteratureRegime::NaturalConvection;
    lit.frequency_hz = p.frequency_hz;
    lit.t_initial_k = p.t_initial_k;
    lit.generator_power_w = p.power.generator_input_power_w; // supply-side; not absorbed target
    lit.target_glass_absorbed_power_w = Real(0);             // unknown until current arrives
    lit.stirrer_rpm = Real(0);
    // Leave enable_boussinesq to literature default / CLI (--zero-buoyancy for isothermal hydro).
    lit.constant_thermal_rheology = true;
    lit.material = buildClientCaseMaterialLaws(p);
    lit.constant_sigma = (p.sigma_preset == OphelieClientSigmaPreset::FrenchTable1Sigma16Borrowed);
}

inline bool writeClientEffectiveCaseJson(const OphelieClientCaseParams &p, const std::string &path)
{
    const auto laws = buildClientCaseMaterialLaws(p);
    const Real t_ref = p.mu_freeze_t_k;
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "{\n";
    out << "  \"schema_version\": \"" << p.schema_version << "\",\n";
    out << "  \"case_id\": \"" << p.case_id << "\",\n";
    out << "  \"material_id\": \"" << p.material_id << "\",\n";
    out << "  \"stirring_enabled\": " << (p.stirring_enabled ? "true" : "false") << ",\n";
    out << "  \"mechanics_mode\": \"" << p.mechanics_mode << "\",\n";
    out << "  \"elastic_paddle_required\": " << (p.elastic_paddle_required ? "true" : "false") << ",\n";
    out << "  \"melt_height_m\": " << p.melt_height_m << ",\n";
    out << "  \"melt_height_status\": \"" << clientConfirmationName(p.melt_height_status) << "\",\n";
    out << "  \"particle_spacing_m\": " << p.particle_spacing_m << ",\n";
    out << "  \"frequency_hz\": " << p.frequency_hz << ",\n";
    out << "  \"frequency_status\": \"" << clientConfirmationName(p.frequency_status) << "\",\n";
    out << "  \"stirrer_rpm\": " << p.stirrer_rpm << ",\n";
    out << "  \"rpm_status\": \"" << clientConfirmationName(p.rpm_status) << "\",\n";
    out << "  \"t_initial_k\": " << p.t_initial_k << ",\n";
    out << "  \"t_initial_status\": \"" << clientConfirmationName(p.t_initial_status) << "\",\n";
    out << "  \"t_initial_note\": \"" << p.t_initial_note << "\",\n";
    out << "  \"probe_target_k\": " << p.probe_target_k << ",\n";
    out << "  \"probe_target_status\": \"" << clientConfirmationName(p.probe_target_status) << "\",\n";
    out << "  \"probe_target_note\": \"" << p.probe_target_note << "\",\n";
    out << "  \"power\": {\n";
    out << "    \"generator_input_power_W\": " << p.power.generator_input_power_w << ",\n";
    out << "    \"generator_power_status\": \"" << clientConfirmationName(p.power.generator_power_status)
        << "\",\n";
    out << "    \"coil_current_A\": " << p.power.coil_current_a << ",\n";
    out << "    \"current_convention\": \"" << p.power.current_convention << "\",\n";
    out << "    \"current_status\": \"" << clientConfirmationName(p.power.current_status) << "\",\n";
    out << "    \"P_em_absorbed_raw_W\": " << p.power.p_em_absorbed_raw_w << ",\n";
    out << "    \"absorbed_power_target_W\": " << p.power.absorbed_power_target_w << ",\n";
    out << "    \"absorbed_power_target_active\": "
        << (p.power.absorbed_power_target_active ? "true" : "false") << ",\n";
    out << "    \"note\": \"30 kW is generator input; do not copy French 50/60 kW absorbed presets\"\n";
    out << "  },\n";
    out << "  \"material_effective\": {\n";
    out << "    \"rho_kg_m3\": " << laws.rho.evaluate(t_ref) << ",\n";
    out << "    \"rho_status\": \"" << clientConfirmationName(p.rho_status) << "\",\n";
    out << "    \"rho_note\": \"" << p.rho_note << "\",\n";
    out << "    \"mu_Pa_s_at_ref\": " << laws.mu.evaluate(t_ref) << ",\n";
    out << "    \"mu_preset\": \""
        << (p.mu_preset == OphelieClientMuPreset::FrozenAtReferenceT ? "frozen_at_reference_T" : "client_formula")
        << "\",\n";
    out << "    \"mu_freeze_T_K\": " << p.mu_freeze_t_k << ",\n";
    out << "    \"mu_source\": \"" << laws.mu_source_ << "\",\n";
    out << "    \"sigma_S_per_m_at_ref\": " << laws.sigma.evaluate(t_ref) << ",\n";
    out << "    \"sigma_source\": \"" << laws.sigma_source_ << "\",\n";
    out << "    \"cp\": " << laws.cp.evaluate(t_ref) << ",\n";
    out << "    \"k\": " << laws.k.evaluate(t_ref) << ",\n";
    out << "    \"beta\": " << laws.beta.evaluate(t_ref) << ",\n";
    out << "    \"cp_k_beta_status\": \"" << clientConfirmationName(p.cp_k_beta_status) << "\"\n";
    out << "  },\n";
    out << "  \"warnings\": ["
        << "\"old pow(10,4.5-8500/T) client sigma interpretation is retired\","
        << "\"do not nail probe_target as Dirichlet BC or domain IC then claim PASS\","
        << "\"elastic_driven stirring case requires Doc228 E0-F0 before software mechanical project\""
        << "]\n";
    out << "}\n";
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLIENT_CASE_PARAMS_H
