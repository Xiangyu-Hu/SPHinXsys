/**
 * @file test_3d_ophelie_doc231_b3_adiabatic_heat.cpp
 * @brief Doc231 B3 — adiabatic prescribed volumetric heat on client NC Reload (no flow).
 *
 * Gate: short-window mean ΔT ≈ P Δt / (ρ cp V); eps_E ≤ 1%.
 * Boundaries closed (E_side=E_bottom=E_top=0). Not EM; not 30 kW generator input.
 */
#include "electromagnetic_ophelie_client_case_params.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_material_laws.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{
inline bool reloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

inline std::string resolveDefaultReload()
{
    const StdVec<std::string> candidates = {
        "./reload",
        "../../test_3d_ophelie_client_nc_glass_relax/bin/reload",
        "../../../test_3d_ophelie_client_nc_glass_relax/bin/reload",
    };
    for (const std::string &c : candidates)
        if (reloadXmlExists(c))
            return c;
    return candidates.front();
}

struct Cli
{
    std::string reload_dir = resolveDefaultReload();
    Real dp = Real(0.010);
    Real melt_height = Real(0.180);
    Real heat_power_w = Real(1.0e4); // development absorbed W (not generator 30 kW)
    Real end_time = Real(0.2);
    Real dt = Real(1.0e-3);
    Real eps_e_max = Real(0.01);
    Real dT_rel_max = Real(0.01);
    std::string out_dir = ".";
    std::string ledger_csv = "thermal_ledger.csv";
    std::string summary_json = "summary.json";
};

void applyCli(int ac, char *av[], Cli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
            cli.reload_dir = av[i] + 13;
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        else if (std::strncmp(av[i], "--melt-height=", 14) == 0)
            cli.melt_height = static_cast<Real>(std::atof(av[i] + 14));
        else if (std::strncmp(av[i], "--heat-power=", 13) == 0)
            cli.heat_power_w = static_cast<Real>(std::atof(av[i] + 13));
        else if (std::strncmp(av[i], "--end-time=", 11) == 0)
            cli.end_time = static_cast<Real>(std::atof(av[i] + 11));
        else if (std::strncmp(av[i], "--dt=", 5) == 0)
            cli.dt = static_cast<Real>(std::atof(av[i] + 5));
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            cli.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--ledger-csv=", 13) == 0)
            cli.ledger_csv = av[i] + 13;
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            cli.summary_json = av[i] + 15;
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc231 B3 adiabatic heat: --reload-dir= --heat-power= --end-time= --out-dir=\n";
            std::exit(0);
        }
    }
}

inline Real hostMeltVolume(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real V = 0;
    for (size_t i = 0; i < n; ++i)
        V += vol[i];
    return V;
}

inline Real hostMeanDeltaT(BaseParticles &particles)
{
    // Prefer OphelieThermalDeltaT: absolute Temperature at ~1273 K loses small ΔT in float32.
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real num = 0, den = 0;
    for (size_t i = 0; i < n; ++i)
    {
        num += dT[i] * vol[i];
        den += vol[i];
    }
    return den > TinyReal ? num / den : Real(0);
}

inline void hostResetThermalState(BaseParticles &particles, Real t0)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        dT[i] = Real(0);
        T[i] = t0;
    }
}

/** Host adiabatic heat step: accumulate ΔT; T = t0 + ΔT (ΔT is ledger truth for float32). */
inline void hostApplyUniformHeat(BaseParticles &particles, Real q_w_m3, Real rho, Real cp, Real t0, Real dt)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Real *T = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real inv = Real(1) / (rho * cp + TinyReal);
    const Real heating = q_w_m3 * dt * inv;
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        dT[i] += heating;
        T[i] = t0 + dT[i];
    }
}

/** H from ΔT field (avoids float32 cancellation in T - t0 at ~1273 K). */
inline Real hostThermalEnergyFromDeltaT(BaseParticles &particles, Real rho, Real cp)
{
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *dT = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real E = 0;
    for (size_t i = 0; i < n; ++i)
        E += rho * cp * dT[i] * vol[i];
    return E;
}
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    std::cerr << "test_3d_ophelie_doc231_b3_adiabatic_heat requires SYCL build.\n";
    return 1;
#else
    Cli cli;
    applyCli(ac, av, cli);
    fs::create_directories(cli.out_dir);

    if (!reloadXmlExists(cli.reload_dir))
    {
        std::cerr << "[ophelie][doc231][B3] FATAL: missing Reload under " << cli.reload_dir << "\n";
        return 2;
    }
    if (!(cli.heat_power_w > TinyReal) || !(cli.end_time > TinyReal) || !(cli.dt > TinyReal))
    {
        std::cerr << "[ophelie][doc231][B3] FATAL: need heat-power>0, end-time>0, dt>0\n";
        return 2;
    }

    OphelieClientCaseParams preset = makeClientNaturalConvectionDevelopmentParams();
    preset.melt_height_m = cli.melt_height;
    preset.particle_spacing_m = cli.dp;
    const auto laws = buildClientCaseMaterialLaws(preset);
    const Real rho = laws.rho.evaluate(preset.t_initial_k);
    const Real cp = laws.cp.evaluate(preset.t_initial_k);
    const Real t0 = preset.t_initial_k;

    OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
    melt.melt_height_m = cli.melt_height;
    melt.glass_half_height = Real(0.5) * melt.melt_height_m;
    melt.glass_center_z_m = melt.inner_bottom_z_m + melt.glass_half_height;
    melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);
    const Real R = melt.melt_radius_m;
    const Real pad = Real(4) * cli.dp;
    BoundingBoxd bb(Vecd(-R, -R, melt.inner_bottom_z_m) - Vecd(pad, pad, pad),
                    Vecd(R, R, melt.inner_bottom_z_m + melt.melt_height_m) + Vecd(pad, pad, pad));

    SPHSystem system(bb, cli.dp);
    system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    FluidBody glass(system, makeShared<OphelieClientNcGlassCylinderShape>("GlassBody", melt, 20));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho, Real(8));
    glass.generateParticles<BaseParticles, Reload>(glass.Name());
    glass.getBaseParticles().registerStateVariable<Real>(kOphelieTemperatureField, t0);
    glass.getBaseParticles().registerStateVariable<Real>(kOphelieThermalDeltaTField, Real(0));
    glass.getBaseParticles().addEvolvingVariable<Real>(kOphelieTemperatureField);
    glass.getBaseParticles().addEvolvingVariable<Real>(kOphelieThermalDeltaTField);

    const size_t n = glass.getBaseParticles().TotalRealParticles();
    if (n == 0)
    {
        std::cerr << "[ophelie][doc231][B3] FATAL: zero particles\n";
        return 2;
    }

    BaseParticles &particles = glass.getBaseParticles();
    hostResetThermalState(particles, t0);
    const Real V = hostMeltVolume(particles);
    const Real q = cli.heat_power_w / std::max(V, TinyReal);
    const Real M_cp = rho * cp * V;

    // Establish H baseline after IC (Doc231: don't count IC replace as time-advance error).
    FrenchAppliedEnergyLedger ledger;
    const Real E0 = hostThermalEnergyFromDeltaT(particles, rho, cp); // ~0 after reset
    ledger.setH0(E0);
    Real t = 0;
    Real E_Q = 0;
    const Real E_side = 0, E_bottom = 0, E_top = 0; // adiabatic

    const fs::path ledger_path = fs::path(cli.out_dir) / cli.ledger_csv;
    std::ofstream led(ledger_path);
    led << std::setprecision(10);
    led << "t_s,dt_s,P_Q_W,E_Q_J,E_side_J,E_bottom_J,E_top_J,H_J,DeltaH_J,R_E_J,eps_E,dT_mean_K\n";

    int step = 0;
    while (t < cli.end_time - TinyReal)
    {
        Real dt_use = cli.dt;
        if (t + dt_use > cli.end_time)
            dt_use = cli.end_time - t;
        hostApplyUniformHeat(particles, q, rho, cp, t0, dt_use);
        E_Q += cli.heat_power_w * dt_use;
        ledger.accumulateQ(cli.heat_power_w, dt_use);
        t += dt_use;
        ++step;

        const Real H = hostThermalEnergyFromDeltaT(particles, rho, cp);
        const Real DeltaH = H - E0;
        // Doc231: DeltaH = E_Q - E_side - E_bottom - E_top + R_E  →  R_E = DeltaH - E_Q (+losses)
        const Real R_E = DeltaH - E_Q + E_side + E_bottom + E_top;
        const Real scale = std::max({std::abs(E_Q), std::abs(E_side) + std::abs(E_bottom) + std::abs(E_top),
                                     std::abs(DeltaH), Real(1)});
        const Real eps_E = std::abs(R_E) / scale;
        const Real dT_mean = hostMeanDeltaT(particles);
        led << t << "," << dt_use << "," << cli.heat_power_w << "," << E_Q << "," << E_side << "," << E_bottom
            << "," << E_top << "," << H << "," << DeltaH << "," << R_E << "," << eps_E << "," << dT_mean << "\n";
    }

    const Real H_end = hostThermalEnergyFromDeltaT(particles, rho, cp);
    const Real DeltaH = H_end - E0;
    const Real R_E = ledger.residualH(H_end);
    const Real eps_E = ledger.residualHRel(H_end);
    const Real dT_num = hostMeanDeltaT(particles);
    const Real dT_from_H = DeltaH / std::max(M_cp, TinyReal);
    const Real dT_analytic = cli.heat_power_w * t / std::max(M_cp, TinyReal);
    const Real dT_rel =
        std::abs(dT_analytic) > TinyReal ? std::abs(dT_num - dT_analytic) / std::abs(dT_analytic) : Real(1);

    const bool finite_ok = std::isfinite(static_cast<double>(dT_num)) && std::isfinite(static_cast<double>(eps_E));
    const bool energy_ok = (eps_E <= cli.eps_e_max) && !ledger.invalid_accumulate;
    const bool dT_ok = (dT_rel <= cli.dT_rel_max) && (dT_num > TinyReal);
    const bool pass = finite_ok && energy_ok && dT_ok && (step > 0);

    std::cout << "[ophelie][doc231][B3] adiabatic prescribed heat n=" << n << " V=" << V
              << " P=" << cli.heat_power_w << " Q=" << q << " t=" << t << " steps=" << step << "\n";
    std::cout << "[ophelie][doc231][B3] dT_num=" << dT_num << " dT_from_H=" << dT_from_H
              << " dT_analytic=" << dT_analytic << " dT_rel=" << dT_rel << " E_Q=" << E_Q
              << " DeltaH=" << DeltaH << " R_E=" << R_E << " eps_E=" << eps_E << "\n";
    std::cout << "[ophelie][doc231][B3] THERMAL_LEDGER_PASS=" << (pass ? 1 : 0)
              << " (adiabatic inlet only; Robin/EM not in this gate)\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc231_b3_adiabatic_heat\",\n";
        js << "  \"tag\": \"THERMAL_LEDGER\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"THERMAL_LEDGER_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"mode\": \"adiabatic_prescribed_uniform_volumetric\",\n";
        js << "  \"note\": \"not EM; not generator 30kW; flow off; boundaries closed; dT from OphelieThermalDeltaT\",\n";
        js << "  \"reload_dir\": \"" << cli.reload_dir << "\",\n";
        js << "  \"material_id\": \"" << preset.material_id << "\",\n";
        js << "  \"n_particles\": " << n << ",\n";
        js << "  \"melt_volume_m3\": " << V << ",\n";
        js << "  \"rho\": " << rho << ",\n";
        js << "  \"cp\": " << cp << ",\n";
        js << "  \"T0_K\": " << t0 << ",\n";
        js << "  \"heat_power_W\": " << cli.heat_power_w << ",\n";
        js << "  \"Q_uniform_Wpm3\": " << q << ",\n";
        js << "  \"physical_time_s\": " << t << ",\n";
        js << "  \"E_Q_J\": " << E_Q << ",\n";
        js << "  \"E_side_J\": 0,\n";
        js << "  \"E_bottom_J\": 0,\n";
        js << "  \"E_top_J\": 0,\n";
        js << "  \"DeltaH_J\": " << DeltaH << ",\n";
        js << "  \"R_E_J\": " << R_E << ",\n";
        js << "  \"eps_E\": " << eps_E << ",\n";
        js << "  \"eps_E_max\": " << cli.eps_e_max << ",\n";
        js << "  \"dT_mean_K\": " << dT_num << ",\n";
        js << "  \"dT_from_H_K\": " << dT_from_H << ",\n";
        js << "  \"dT_analytic_K\": " << dT_analytic << ",\n";
        js << "  \"dT_rel_err\": " << dT_rel << ",\n";
        js << "  \"ledger_csv\": \"" << ledger_path.string() << "\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][doc231][B3] ledger=" << ledger_path.string()
              << " summary=" << summary_path.string() << "\n";
    return pass ? 0 : 1;
#endif
}
