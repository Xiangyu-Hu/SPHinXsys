/**
 * @file test_3d_ophelie_robin_plate_halfspace.cpp
 * @brief Doc207 A: frozen-flow plate Robin vs semi-infinite analytic (production thermal kernels).
 *
 * Domain box [0,L]^3 (default L=0.12). Face x=0: Robin via
 *   --thermal-wall-model=legacy-shell (default): ApplyOphelieThermalFrenchNaturalBcCK
 *   --thermal-wall-model=wall-temperature (Doc207 B): ApplyOphelieThermalWallTemperatureBcCK
 *     with geometric gap d=x-x_wall (series k(Ti-Tw)/d = h(Tw-Tinf)).
 * Planar shell area weights (Vol/shell + scale). Other faces adiabatic.
 * Laplace = OpheliePairwiseLaplaceCK + ApplyLaplaceDiffusionDtCK. Frozen flow; no Q/radiation.
 *
 * NOT FRENCH_VALIDATED. Compares to independent analytic (scripts/robin_halfspace_analytic.py).
 */
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
struct Cli
{
    Real dp = 0.015;
    Real L = 0.12;
    Real rho = 2750.0;
    Real cp = 1150.0;
    Real k = 4.0;
    Real Ti = 1473.0;
    Real Tinf = 303.0;
    Real h = 300.0;
    Real end_time = 50.0;
    Real dt = 0.05;
    Real pair_weight_regularization = Real(0.01);
    bool area_normalize = true;
    std::string thermal_wall_model = "legacy-shell"; // or wall-temperature
    std::string output_dir;
    std::string dump_times = "1,5,10,20,50";
    std::string summary_json = "summary.json";
};

bool parseCli(int argc, char **argv, Cli &cli)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto take = [&](const std::string &key, Real &dst) {
            const std::string p = key + "=";
            if (arg.rfind(p, 0) == 0)
            {
                dst = static_cast<Real>(std::atof(arg.c_str() + p.size()));
                return true;
            }
            return false;
        };
        auto take_s = [&](const std::string &key, std::string &dst) {
            const std::string p = key + "=";
            if (arg.rfind(p, 0) == 0)
            {
                dst = arg.substr(p.size());
                return true;
            }
            return false;
        };
        if (take("--dp", cli.dp) || take("--L", cli.L) || take("--rho", cli.rho) || take("--cp", cli.cp) ||
            take("--k", cli.k) || take("--Ti", cli.Ti) || take("--Tinf", cli.Tinf) || take("--h", cli.h) ||
            take("--end-time", cli.end_time) || take("--dt", cli.dt) ||
            take("--pair-weight-regularization", cli.pair_weight_regularization) ||
            take_s("--output-dir", cli.output_dir) || take_s("--dump-times", cli.dump_times) ||
            take_s("--summary-json", cli.summary_json) || take_s("--thermal-wall-model", cli.thermal_wall_model))
            continue;
        if (arg == "--area-normalize")
        {
            cli.area_normalize = true;
            continue;
        }
        if (arg == "--no-area-normalize")
        {
            cli.area_normalize = false;
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            std::cout
                << "Doc207 A/B robin plate halfspace\n"
                << "  --dp=0.015 --L=0.12 --rho=2750 --cp=1150 --k=4 --Ti=1473 --Tinf=303 --h=300\n"
                << "  --end-time=50 --dt=0.05 --area-normalize|--no-area-normalize\n"
                << "  --thermal-wall-model=legacy-shell|wall-temperature\n"
                << "  --dump-times=1,5,10,20,50 --output-dir= --summary-json=\n";
            return false;
        }
        std::cerr << "[doc207a] unknown arg: " << arg << std::endl;
        return false;
    }
    if (!(cli.dp > 0 && cli.L > 0 && cli.rho > 0 && cli.cp > 0 && cli.k > 0 && cli.h >= 0 && cli.dt > 0 &&
          cli.end_time >= 0))
    {
        std::cerr << "[doc207a] invalid numeric args\n";
        return false;
    }
    if (cli.thermal_wall_model != "legacy-shell" && cli.thermal_wall_model != "wall-temperature" &&
        cli.thermal_wall_model != "wall_temperature" && cli.thermal_wall_model != "series-gap")
    {
        std::cerr << "[doc207a] --thermal-wall-model must be legacy-shell|wall-temperature\n";
        return false;
    }
    return true;
}

std::vector<Real> parseTimes(const std::string &s)
{
    std::vector<Real> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (tok.empty())
            continue;
        out.push_back(static_cast<Real>(std::atof(tok.c_str())));
    }
    std::sort(out.begin(), out.end());
    return out;
}

struct PlateAreaAudit
{
    double A_target = 0;
    double A_raw = 0;
    double scale = 1;
    int n_shell = 0;
};

/**
 * Production B3 convention on planar face x=0.
 * Kernel ledger area = Vol/shell; w = scale on shell.
 * Also writes geometric wall gap d = max(x - x0, 0) for Doc207 B wall-temperature model.
 */
inline PlateAreaAudit hostBuildPlateFaceWeights(BaseParticles &particles, Real x0, Real shell, Real A_target,
                                                bool normalize)
{
    PlateAreaAudit a;
    a.A_target = static_cast<double>(A_target);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, kOphelieThermalWallGapField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mask = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryMaskField);
    Real *face = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryFaceField);
    Real *w_side = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightSideField);
    Real *w_bottom = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightBottomField);
    Real *w_free = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryWeightFreeField);
    Real *gap = particles.getVariableDataByName<Real>(kOphelieThermalWallGapField);
    const size_t n = particles.TotalRealParticles();
    std::vector<char> on(n, 0);
    const double inv_shell = 1.0 / (static_cast<double>(shell) + 1e-30);
    for (size_t i = 0; i < n; ++i)
    {
        gap[i] = Real(0);
        const Real x = pos[i][0] - x0;
        if (x < shell)
        {
            on[i] = 1;
            a.A_raw += static_cast<double>(vol[i]) * inv_shell;
            ++a.n_shell;
            gap[i] = std::max(x, Real(0)); // geometric particle→face distance
        }
    }
    a.scale = (normalize && a.A_raw > 1e-30) ? a.A_target / a.A_raw : 1.0;
    for (size_t i = 0; i < n; ++i)
    {
        w_side[i] = w_bottom[i] = w_free[i] = Real(0);
        mask[i] = on[i] ? Real(1) : Real(0);
        face[i] = on[i] ? Real(1) : Real(0); // 1 = side channel
        if (on[i])
            w_side[i] = static_cast<Real>(a.scale);
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryFaceField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightSideField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightBottomField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryWeightFreeField);
    syncVariableToDevice<Real>(particles, kOphelieThermalWallGapField);
    return a;
}

inline void hostInitUniformT(BaseParticles &particles, Real Ti, Real t_ref)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        temp[i] = Ti;
        delta[i] = Ti - t_ref;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

struct DumpStats
{
    double T_shell_mean = 0;
    double Tw_recon_mean = 0;
    double T_far_mean = 0;
    double H_J = 0; // rho cp integral (T-Ti) * V  (negative when cooling)
    int n_shell = 0;
    int n_far = 0;
};

inline DumpStats hostDumpProfile(BaseParticles &particles, const std::string &path, Real x0, Real L, Real shell,
                                 Real rho, Real cp, Real Ti, bool have_tw)
{
    DumpStats st;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    if (have_tw)
        syncVariableToHost<Real>(particles, kOphelieThermalWallTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalWallGapField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *tw = have_tw ? particles.getVariableDataByName<Real>(kOphelieThermalWallTemperatureField) : nullptr;
    const Real *gap = particles.getVariableDataByName<Real>(kOphelieThermalWallGapField);
    const size_t n = particles.TotalRealParticles();
    std::ofstream ofs(path);
    ofs << "x_m,y_m,z_m,T_K,vol_m3,Tw_recon_K,wall_gap_m\n";
    double ms = 0, mTs = 0, mTw = 0, mf = 0, mTf = 0, H = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - x0;
        const Real tw_i = (have_tw && tw) ? tw[i] : Real(0);
        ofs << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << temp[i] << "," << vol[i] << "," << tw_i
            << "," << gap[i] << "\n";
        H += static_cast<double>(rho * cp * (temp[i] - Ti) * vol[i]);
        if (x < shell)
        {
            ms += vol[i];
            mTs += temp[i] * vol[i];
            mTw += tw_i * vol[i];
            ++st.n_shell;
        }
        if (x > L - 2 * shell)
        {
            mf += vol[i];
            mTf += temp[i] * vol[i];
            ++st.n_far;
        }
    }
    st.T_shell_mean = ms > 0 ? mTs / ms : 0;
    st.Tw_recon_mean = ms > 0 ? mTw / ms : 0;
    st.T_far_mean = mf > 0 ? mTf / mf : 0;
    st.H_J = H;
    return st;
}
} // namespace

int main(int argc, char **argv)
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
        return 1;

    if (!cli.output_dir.empty())
    {
        std::filesystem::create_directories(cli.output_dir);
        if (cli.summary_json.find('/') == std::string::npos)
            cli.summary_json = cli.output_dir + "/" + cli.summary_json;
    }

    const Real alpha = cli.k / (cli.rho * cli.cp + TinyReal);
    const Real x0 = 0.0;
    const Real L = cli.L;
    const Real dp = cli.dp;
    const Real shell = dp; // production French B3 shell = dp
    const Real A_face = L * L;
    const Real dt_limit = Real(0.1) * dp * dp / (alpha + TinyReal); // conservative 3D-ish
    if (cli.dt > dt_limit)
    {
        std::cerr << "[doc207a] FATAL: --dt=" << cli.dt << " exceeds ~0.1 dp^2/alpha=" << dt_limit << std::endl;
        return 1;
    }

    const Vecd center(Real(0.5) * L, Real(0.5) * L, Real(0.5) * L);
    const Vecd halfsize(Real(0.5) * L, Real(0.5) * L, Real(0.5) * L);
    const BoundingBoxd bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));
    SPHSystem sph_system(bounds, dp);
    SolidBody body(sph_system, makeShared<OphelieTestGlassBoxShape>("RobinPlate", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMatterMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const Real t_ref = cli.Ti; // DeltaT=0 initially
    registerOphelieJouleHeatTemperatureField(particles, t_ref);
    registerOphelieThermalDiffusionAuxFields(particles, cli.k);
    UniquePtr<Inner<>> inner = makeUnique<Inner<>>(body);

    hostInitUniformT(particles, cli.Ti, t_ref);
    const PlateAreaAudit area = hostBuildPlateFaceWeights(particles, x0, shell, A_face, cli.area_normalize);

    OphelieThermalDiffusionOneWayOptions bc;
    bc.enable_french_natural_bc = true;
    bc.h_side = cli.h;
    bc.h_bottom = Real(0);
    bc.h_free = Real(0);
    bc.emissivity = Real(0);
    bc.t_cool = cli.Tinf;
    bc.t_ambient = cli.Tinf;
    bc.t_rad_ambient = cli.Tinf;
    bc.boundary_shell_thickness = shell;
    bc.pair_weight_regularization = cli.pair_weight_regularization;
    bc.thermal_wall_model = cli.thermal_wall_model;
    const bool use_wall_T = ophelieThermalWallModelIsWallTemperature(bc);

    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace(
        *inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
        cli.pair_weight_regularization);
    StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK> apply_diff(body, cli.rho, cli.cp, t_ref);
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK>> apply_legacy;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcCK>> apply_wall;
    UniquePtr<ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK>> loss_legacy;
    UniquePtr<ReduceDynamicsCK<MainExecutionPolicy, OphelieWallTemperatureBoundaryLossPowerReduceCK>> loss_wall;
    if (use_wall_T)
    {
        apply_wall = makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalWallTemperatureBcCK>>(
            body, kOphelieThermalDeltaTField, kOphelieTemperatureField, t_ref, cli.rho, cli.cp, shell, bc);
        loss_wall =
            makeUnique<ReduceDynamicsCK<MainExecutionPolicy, OphelieWallTemperatureBoundaryLossPowerReduceCK>>(
                body, t_ref, shell, bc);
    }
    else
    {
        apply_legacy = makeUnique<StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK>>(
            body, kOphelieThermalDeltaTField, kOphelieTemperatureField, t_ref, cli.rho, cli.cp, shell, bc);
        loss_legacy = makeUnique<ReduceDynamicsCK<MainExecutionPolicy, OphelieFrenchBoundaryLossPowerReduceCK>>(
            body, t_ref, shell, bc);
    }
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cells(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner(*inner);

    std::vector<Real> dumps = parseTimes(cli.dump_times);
    std::vector<char> dump_done(dumps.size(), 0);

    const std::string mon_path =
        cli.output_dir.empty() ? std::string("doc207a_monitor.csv") : (cli.output_dir + "/monitor.csv");
    std::ofstream mon(mon_path);
    mon << "t_s,P_loss_W,H_J,T_shell_mean,Tw_recon_mean,T_far_mean,E_loss_J\n";

    const char *bc_path =
        use_wall_T ? "ApplyOphelieThermalWallTemperatureBcCK(series-gap)" : "ApplyOphelieThermalFrenchNaturalBcCK";
    std::cout << "[doc207a] start n=" << particles.TotalRealParticles() << " L=" << L << " dp=" << dp
              << " shell=" << shell << " A_target=" << area.A_target << " A_raw=" << area.A_raw
              << " scale=" << area.scale << " n_shell=" << area.n_shell << " alpha=" << alpha << " dt=" << cli.dt
              << " dt_limit~" << dt_limit << " end=" << cli.end_time
              << " area_normalize=" << (cli.area_normalize ? 1 : 0)
              << " thermal_wall_model=" << cli.thermal_wall_model
              << " path=OpheliePairwiseLaplaceCK+ApplyLaplaceDiffusionDtCK+" << bc_path << std::endl
              << std::flush;

    Real t = 0.0;
    double E_loss = 0.0; // integral P dt
    Real P_last = 0.0;
    DumpStats st0 =
        hostDumpProfile(particles,
                        (cli.output_dir.empty() ? std::string("profile_t0.csv") : cli.output_dir + "/profile_t0.csv"),
                        x0, L, shell, cli.rho, cli.cp, cli.Ti, use_wall_T);
    mon << "0," << 0 << "," << st0.H_J << "," << st0.T_shell_mean << "," << st0.Tw_recon_mean << "," << st0.T_far_mean
        << ",0\n";

    auto maybe_dump = [&](Real tnow) {
        for (size_t k = 0; k < dumps.size(); ++k)
        {
            if (dump_done[k])
                continue;
            if (tnow + Real(1e-9) >= dumps[k])
            {
                std::ostringstream fn;
                if (!cli.output_dir.empty())
                    fn << cli.output_dir << "/";
                fn << "profile_t" << dumps[k] << ".csv";
                DumpStats st =
                    hostDumpProfile(particles, fn.str(), x0, L, shell, cli.rho, cli.cp, cli.Ti, use_wall_T);
                std::cout << "[doc207a] dump t=" << tnow << " -> " << fn.str() << " T_shell=" << st.T_shell_mean
                          << " Tw_recon=" << st.Tw_recon_mean << " T_far=" << st.T_far_mean << " H=" << st.H_J
                          << std::endl
                          << std::flush;
                dump_done[k] = 1;
            }
        }
    };
    maybe_dump(0.0);

    while (t + TinyReal < cli.end_time)
    {
        const Real dt = std::min(cli.dt, cli.end_time - t);
        update_cells.exec();
        update_inner.exec();
        laplace.exec();
        apply_diff.exec(dt);
        // F02: integrate loss on the same T the BC operator applies (before update).
        if (use_wall_T)
        {
            P_last = loss_wall->exec();
            apply_wall->exec(dt);
        }
        else
        {
            P_last = loss_legacy->exec();
            apply_legacy->exec(dt);
        }
        E_loss += static_cast<double>(P_last) * static_cast<double>(dt);
        t += dt;

        syncVariableToHost<Real>(particles, kOphelieTemperatureField);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<Vecd>(particles, "Position");
        if (use_wall_T)
            syncVariableToHost<Real>(particles, kOphelieThermalWallTemperatureField);
        const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Real *tw =
            use_wall_T ? particles.getVariableDataByName<Real>(kOphelieThermalWallTemperatureField) : nullptr;
        const size_t n = particles.TotalRealParticles();
        double H = 0, ms = 0, mTs = 0, mTw = 0, mf = 0, mTf = 0;
        for (size_t i = 0; i < n; ++i)
        {
            H += static_cast<double>(cli.rho * cli.cp * (temp[i] - cli.Ti) * vol[i]);
            const Real x = pos[i][0] - x0;
            if (x < shell)
            {
                ms += vol[i];
                mTs += temp[i] * vol[i];
                if (tw)
                    mTw += tw[i] * vol[i];
            }
            if (x > L - 2 * shell)
            {
                mf += vol[i];
                mTf += temp[i] * vol[i];
            }
        }
        const double Ts = ms > 0 ? mTs / ms : 0;
        const double Twm = ms > 0 ? mTw / ms : 0;
        const double Tf = mf > 0 ? mTf / mf : 0;
        mon << t << "," << P_last << "," << H << "," << Ts << "," << Twm << "," << Tf << "," << E_loss << "\n";
        maybe_dump(t);
        if (static_cast<int>(t / cli.dt + Real(0.5)) % 50 == 0)
        {
            std::cout << "[doc207a] t=" << t << " P=" << P_last << " T_shell=" << Ts << " Tw_recon=" << Twm
                      << " T_far=" << Tf << " H=" << H << std::endl
                      << std::flush;
        }
    }

    // final dump if end_time in list
    maybe_dump(t);

    const double H_end = [&]() {
        syncVariableToHost<Real>(particles, kOphelieTemperatureField);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        double H = 0;
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            H += static_cast<double>(cli.rho * cli.cp * (temp[i] - cli.Ti) * vol[i]);
        return H;
    }();

    // Energy: cooling → H_end < 0; loss ledger E_loss > 0; residual |H_end + E_loss|
    const double resid = std::abs(H_end + E_loss);
    const double scale = std::max(std::abs(E_loss), 1.0);
    const bool energy_ok = resid / scale <= 0.02;

    std::ofstream js(cli.summary_json);
    js << "{\n"
       << "  \"case\": \"doc207a_robin_plate_halfspace\",\n"
       << "  \"NOT_FRENCH_VALIDATED\": 1,\n"
       << "  \"dp\": " << dp << ",\n"
       << "  \"L\": " << L << ",\n"
       << "  \"shell_m\": " << shell << ",\n"
       << "  \"A_target_m2\": " << area.A_target << ",\n"
       << "  \"A_raw_m2\": " << area.A_raw << ",\n"
       << "  \"area_scale\": " << area.scale << ",\n"
       << "  \"n_shell\": " << area.n_shell << ",\n"
       << "  \"area_normalize\": " << (cli.area_normalize ? 1 : 0) << ",\n"
       << "  \"thermal_wall_model\": \"" << cli.thermal_wall_model << "\",\n"
       << "  \"rho\": " << cli.rho << ",\n"
       << "  \"cp\": " << cli.cp << ",\n"
       << "  \"k\": " << cli.k << ",\n"
       << "  \"alpha\": " << alpha << ",\n"
       << "  \"Ti\": " << cli.Ti << ",\n"
       << "  \"Tinf\": " << cli.Tinf << ",\n"
       << "  \"h\": " << cli.h << ",\n"
       << "  \"end_time\": " << t << ",\n"
       << "  \"dt\": " << cli.dt << ",\n"
       << "  \"P_loss_W_end\": " << P_last << ",\n"
       << "  \"E_loss_J_ledger\": " << E_loss << ",\n"
       << "  \"H_end_J\": " << H_end << ",\n"
       << "  \"energy_residual_rel\": " << (resid / scale) << ",\n"
       << "  \"energy_ok_2pct\": " << (energy_ok ? 1 : 0) << ",\n"
       << "  \"NO_INDEPENDENT_WALL_TEMPERATURE\": " << (use_wall_T ? 0 : 1) << ",\n"
       << "  \"RECONSTRUCTED_WALL_TEMPERATURE\": " << (use_wall_T ? 1 : 0) << ",\n"
       << "  \"wall_gap_definition\": \"d=max(x-x_wall,0) geometric\",\n"
       << "  \"conduction_path\": \"OpheliePairwiseLaplaceCK+ApplyLaplaceDiffusionDtCK+"
       << (use_wall_T ? "ApplyOphelieThermalWallTemperatureBcCK" : "ApplyOphelieThermalFrenchNaturalBcCK")
       << "\",\n"
       << "  \"note\": \""
       << (use_wall_T ? "Tw from series-gap reconstruction; compare Tw_recon and q/E to analytic."
                      : "Shell-particle T is NOT analytic Tw; compare via Doc207 scripts.")
       << "\"\n"
       << "}\n";

    std::cout << "[doc207a] end t=" << t << " P_end=" << P_last << " E_loss=" << E_loss << " H_end=" << H_end
              << " resid_rel=" << (resid / scale) << " energy_ok=" << (energy_ok ? 1 : 0)
              << " model=" << cli.thermal_wall_model << " summary=" << cli.summary_json << std::endl;
    return energy_ok ? 0 : 2;
}
