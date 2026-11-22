/**
 * @file test_3d_ophelie_doc231_b1_probe_sampling.cpp
 * @brief Doc231 B1 — fixed-space probe sampling gates (NC glass Reload).
 *
 * 1) Uniform T → kernel sample equals T.
 * 2) Linear T(z) → relative error ≤ 1% of amplitude (interior probes).
 * Probe location status remains unconfirmed_sensor_location (separate from algo gate).
 */
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_fixed_space_probe.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

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
        "../reload",
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
    Real err_rel_max = Real(0.01); // Doc231 B1 development gate
    std::string out_dir = ".";
    std::string csv_name = "probe_sampling_check.csv";
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
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            cli.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--csv=", 6) == 0)
            cli.csv_name = av[i] + 6;
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            cli.summary_json = av[i] + 15;
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc231 B1 probe sampling: --reload-dir= --dp= --melt-height= --out-dir=\n";
            std::exit(0);
        }
    }
}

using ProbeSample = OphelieFixedSpaceProbeSample;
using ProbeDef = OphelieFixedSpaceProbeDef;

inline ProbeSample kernelSampleT(BaseParticles &particles, const Vecd &probe, Real kernel_h)
{
    syncFixedSpaceProbeFieldsToHost(particles, "Temperature");
    return hostSampleFixedSpaceTemperature(particles, probe, kernel_h, "Temperature");
}
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    std::cerr << "test_3d_ophelie_doc231_b1_probe_sampling requires SYCL build.\n";
    return 1;
#else
    Cli cli;
    applyCli(ac, av, cli);
    fs::create_directories(cli.out_dir);

    if (!reloadXmlExists(cli.reload_dir))
    {
        std::cerr << "[ophelie][doc231][B1] FATAL: missing Reload.xml under " << cli.reload_dir << "\n";
        return 2;
    }

    OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
    melt.melt_height_m = cli.melt_height;
    melt.glass_half_height = Real(0.5) * melt.melt_height_m;
    melt.glass_center_z_m = melt.inner_bottom_z_m + melt.glass_half_height;
    melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);

    const OphelieClientProbeDefinition thermo =
        makeClientNominalThermocoupleProbe(makeClientNaturalConvectionAssemblyCatalog());
    const Vecd tip(thermo.solver_m.x, thermo.solver_m.y, thermo.solver_m.z);

    // Interior development probes (axis / mid-radius) — not claimed as sensor locations.
    const Real z_bot = melt.inner_bottom_z_m;
    const Real z_top = z_bot + melt.melt_height_m;
    const Real z_mid = Real(0.5) * (z_bot + z_top);
    const Real R = melt.melt_radius_m;
    StdVec<ProbeDef> probes = {
        {"thermocouple_nominal_tip", tip, thermo.status, true,
         std::string("CAD AABB tip candidate")},
        {"dev_axis_mid", Vecd(0, 0, z_mid), "development", true, "development interior axis mid-height"},
        {"dev_axis_lower_quartile", Vecd(0, 0, z_bot + Real(0.25) * melt.melt_height_m), "development", true,
         "development interior"},
        {"dev_r_half_mid", Vecd(Real(0.5) * R, 0, z_mid), "development", true, "development mid-radius"},
    };

    const Real pad = Real(4) * cli.dp;
    BoundingBoxd bb = BoundingBoxd(Vecd(-R, -R, z_bot) - Vecd(pad, pad, pad),
                                   Vecd(R, R, z_top) + Vecd(pad, pad, pad));

    SPHSystem system(bb, cli.dp);
    system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    // Body name must match Reload "GlassBody"; particles from client NC relax.
    FluidBody glass(system, makeShared<OphelieClientNcGlassCylinderShape>("GlassBody", melt, 20));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(Real(2700), Real(8));
    glass.generateParticles<BaseParticles, Reload>(glass.Name());
    glass.getBaseParticles().registerStateVariable<Real>("Temperature", Real(0));
    glass.getBaseParticles().addEvolvingVariable<Real>("Temperature");

    const size_t n = glass.getBaseParticles().TotalRealParticles();
    if (n == 0)
    {
        std::cerr << "[ophelie][doc231][B1] FATAL: zero glass particles in Reload\n";
        return 2;
    }

    const Real kernel_h = Real(1.3) * cli.dp; // match typical adaptation ratio
    Vecd *pos = glass.getBaseParticles().getVariableDataByName<Vecd>("Position");
    Real *temp = glass.getBaseParticles().getVariableDataByName<Real>("Temperature");

    const Real T_uniform = Real(1273.15);
    for (size_t i = 0; i < n; ++i)
        temp[i] = T_uniform;

    const Real T_bot = Real(1200);
    const Real T_top = Real(1300);
    const Real amp = T_top - T_bot;
    auto T_linear = [&](Real z) {
        const Real a = (z - z_bot) / std::max(melt.melt_height_m, TinyReal);
        return T_bot + amp * std::min(Real(1), std::max(Real(0), a));
    };

    const fs::path csv_path = fs::path(cli.out_dir) / cli.csv_name;
    std::ofstream csv(csv_path);
    csv << std::setprecision(10);
    csv << "probe_id,x,y,z,field,T_exact,T_sample,abs_err,rel_err_vs_amp,n_neighbors,weight_sum,"
           "nearest_dist_m,expect_interior,pass_row,note\n";

    bool uniform_ok = true;
    bool linear_ok = true;
    bool support_ok = true;

    for (const ProbeDef &p : probes)
    {
        const ProbeSample s = kernelSampleT(glass.getBaseParticles(), p.pos, kernel_h);
        const Real abs_err = std::abs(s.T - T_uniform);
        const Real rel = abs_err / std::max(amp, Real(1)); // vs Doc231 amp scale even for uniform
        const bool row_pass = std::isfinite(static_cast<double>(s.T)) && (abs_err <= Real(1.0e-3)) &&
                              (s.n_neighbors >= 4) && (s.weight_sum > TinyReal);
        if (!row_pass)
            uniform_ok = false;
        if (s.n_neighbors < 4)
            support_ok = false;
        csv << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ",uniform," << T_uniform
            << "," << s.T << "," << abs_err << "," << rel << "," << s.n_neighbors << "," << s.weight_sum
            << "," << s.nearest_dist << "," << (p.expect_interior ? 1 : 0) << "," << (row_pass ? 1 : 0)
            << "," << p.note << "\n";
    }

    for (size_t i = 0; i < n; ++i)
        temp[i] = T_linear(pos[i][2]);

    for (const ProbeDef &p : probes)
    {
        if (!p.expect_interior)
            continue;
        const ProbeSample s = kernelSampleT(glass.getBaseParticles(), p.pos, kernel_h);
        const Real T_ex = T_linear(p.pos[2]);
        const Real abs_err = std::abs(s.T - T_ex);
        const Real rel = abs_err / std::max(amp, TinyReal);
        const bool row_pass = std::isfinite(static_cast<double>(s.T)) && (rel <= cli.err_rel_max) &&
                              (s.n_neighbors >= 4) && (s.weight_sum > TinyReal);
        if (!row_pass)
            linear_ok = false;
        if (s.n_neighbors < 4)
            support_ok = false;
        csv << p.id << "," << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ",linear_z," << T_ex << ","
            << s.T << "," << abs_err << "," << rel << "," << s.n_neighbors << "," << s.weight_sum << ","
            << s.nearest_dist << ",1," << (row_pass ? 1 : 0) << "," << p.note << "\n";
    }

    // Algo gate only; sensor confirmation is separate.
    const bool probe_algo_pass = uniform_ok && linear_ok && support_ok;
    const bool sensor_confirmed = false; // status unconfirmed_sensor_location

    std::cout << "[ophelie][doc231][B1] n_glass=" << n << " reload=" << cli.reload_dir
              << " tip=(" << tip[0] << "," << tip[1] << "," << tip[2] << ") status=" << thermo.status
              << " kernel_h=" << kernel_h << "\n";
    std::cout << "[ophelie][doc231][B1] uniform_ok=" << (uniform_ok ? 1 : 0)
              << " linear_ok=" << (linear_ok ? 1 : 0) << " support_ok=" << (support_ok ? 1 : 0)
              << " PROBE_SAMPLING_PASS=" << (probe_algo_pass ? 1 : 0)
              << " sensor_confirmed=" << (sensor_confirmed ? 1 : 0) << "\n";
    std::cout << "[ophelie][doc231][B1] csv=" << csv_path.string() << "\n";

    const fs::path summary_path = fs::path(cli.summary_json).is_absolute()
                                      ? fs::path(cli.summary_json)
                                      : fs::path(cli.out_dir) / cli.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc231_b1_probe_sampling\",\n";
        js << "  \"tag\": \"PROBE_SAMPLING\",\n";
        js << "  \"passed\": " << (probe_algo_pass ? 1 : 0) << ",\n";
        js << "  \"PROBE_SAMPLING_PASS\": " << (probe_algo_pass ? 1 : 0) << ",\n";
        js << "  \"sensor_location_confirmed\": " << (sensor_confirmed ? 1 : 0) << ",\n";
        js << "  \"probe_status\": \"" << thermo.status << "\",\n";
        js << "  \"reload_dir\": \"" << cli.reload_dir << "\",\n";
        js << "  \"n_glass\": " << n << ",\n";
        js << "  \"dp_m\": " << cli.dp << ",\n";
        js << "  \"kernel_h_m\": " << kernel_h << ",\n";
        js << "  \"err_rel_max\": " << cli.err_rel_max << ",\n";
        js << "  \"uniform_ok\": " << (uniform_ok ? 1 : 0) << ",\n";
        js << "  \"linear_ok\": " << (linear_ok ? 1 : 0) << ",\n";
        js << "  \"support_ok\": " << (support_ok ? 1 : 0) << ",\n";
        js << "  \"tip_solver_m\": [" << tip[0] << ", " << tip[1] << ", " << tip[2] << "],\n";
        js << "  \"csv\": \"" << csv_path.string() << "\",\n";
        js << "  \"note\": \"Algo gate only; sensor junction still unconfirmed\"\n";
        js << "}\n";
    }
    return probe_algo_pass ? 0 : 1;
#endif
}
