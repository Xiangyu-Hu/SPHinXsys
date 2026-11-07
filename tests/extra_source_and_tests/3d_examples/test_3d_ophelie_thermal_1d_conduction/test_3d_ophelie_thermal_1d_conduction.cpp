/**
 * @file test_3d_ophelie_thermal_1d_conduction.cpp
 * @brief Doc77-B: 1D conduction with production ΔT-authoritative Laplace diffusion path.
 *
 * Domain 0≤x≤L (thin in y,z). Ends Dirichlet Th/Tc; lateral faces insulated (no clamp).
 * Analytic: Tss = Th+(Tc−Th)x/L,
 *           T(x,t)=Tss(x)+A sin(πx/L) exp(−α π² t / L²), α=k/(ρ cp).
 *
 * Reuses: OpheliePairwiseLaplaceCK + ApplyLaplaceDiffusionDtCK (same as production natural/stir).
 * Absolute T0≈1473 K background checks small-increment accumulation on ΔT.
 */
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "sphinxsys.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
struct Cli
{
    Real dp = 0.02;
    Real L = 0.2;
    Real half_yz = 0.04;
    Real rho = 2500.0;
    Real cp = 1200.0;
    Real k = 200.0; /**< elevated for short smoke; still ΔT-authoritative path */
    Real t_h = 1483.0; /**< left */
    Real t_c = 1463.0; /**< right */
    Real t0 = 1473.0;  /**< ΔT reference (production glass) */
    Real amp = 8.0;
    Real end_time = 12.0;
    Real dt = 0.1;
    Real boundary_width_factor = 1.2;
    /** Match production / cavity default (Doc92). Was hardcoded 1.0 → under-diffusion. */
    Real pair_weight_regularization = Real(0.01);
    std::string out_csv = "output/doc77b_1d_conduction.csv";
    std::string summary_json;
};

bool parseCli(int argc, char *argv[], Cli &cli)
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
        if (take("--dp", cli.dp) || take("--L", cli.L) || take("--half-yz", cli.half_yz) || take("--rho", cli.rho) ||
            take("--cp", cli.cp) || take("--k", cli.k) || take("--Th", cli.t_h) || take("--Tc", cli.t_c) ||
            take("--T0", cli.t0) || take("--A", cli.amp) || take("--end-time", cli.end_time) || take("--dt", cli.dt) ||
            take("--pair-weight-regularization", cli.pair_weight_regularization) ||
            take_s("--out-csv", cli.out_csv) || take_s("--summary-json", cli.summary_json))
        {
            continue;
        }
        if (arg == "--help" || arg == "-h")
        {
            std::cout << "Doc77-B 1D conduction\n"
                      << "  --dp= --L= --k= --Th= --Tc= --T0= --A= --end-time= --dt=\n"
                      << "  --pair-weight-regularization= (default 0.01, production)\n"
                      << "  --out-csv= --summary-json=\n";
            return false;
        }
        std::cerr << "[doc77b] unknown arg: " << arg << std::endl;
        return false;
    }
    return true;
}

inline Real analyticT(Real x, Real t, Real L, Real th, Real tc, Real amp, Real alpha)
{
    const Real tss = th + (tc - th) * (x / (L + TinyReal));
    return tss + amp * std::sin(Pi * x / (L + TinyReal)) * std::exp(-alpha * Pi * Pi * t / (L * L + TinyReal));
}

inline void hostApplyEndDirichlet(BaseParticles &particles, Real x0, Real L, Real shell, Real th, Real tc, Real t0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - x0;
        if (x <= shell)
        {
            temp[i] = th;
            delta[i] = th - t0;
        }
        else if (x >= L - shell)
        {
            temp[i] = tc;
            delta[i] = tc - t0;
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

inline void hostInitAnalytic(BaseParticles &particles, Real x0, Real L, Real th, Real tc, Real amp, Real alpha,
                             Real t, Real t0)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = std::min(L, std::max(Real(0), pos[i][0] - x0));
        const Real Ti = analyticT(x, t, L, th, tc, amp, alpha);
        temp[i] = Ti;
        delta[i] = Ti - t0;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
}

struct ErrStats
{
    Real l2 = 0.0;
    Real linf = 0.0;
    Real amp_num = 0.0;
    Real amp_an = 0.0;
    Real energy = 0.0;
};

inline ErrStats hostErrors(BaseParticles &particles, Real x0, Real L, Real th, Real tc, Real amp, Real alpha, Real t,
                           Real t0, Real rho, Real cp, Real shell)
{
    ErrStats e;
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    double err2 = 0.0, ref2 = 0.0, amp_n = 0.0, amp_d = 0.0, E = 0.0;
    Real linf = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = std::min(L, std::max(Real(0), pos[i][0] - x0));
        // Exclude Dirichlet shells from L2 / amplitude (boundary enforced).
        if (x <= shell || x >= L - shell)
        {
            continue;
        }
        const Real Tan = analyticT(x, t, L, th, tc, amp, alpha);
        const Real Terr = temp[i] - Tan;
        err2 += static_cast<double>(Terr) * Terr * vol[i];
        ref2 += static_cast<double>(Tan - t0) * (Tan - t0) * vol[i];
        linf = std::max(linf, std::abs(Terr));
        const Real s = std::sin(Pi * x / (L + TinyReal));
        const Real tss = th + (tc - th) * (x / (L + TinyReal));
        amp_n += (temp[i] - tss) * s * vol[i];
        amp_d += s * s * vol[i];
        E += rho * cp * (temp[i] - t0) * vol[i];
    }
    e.l2 = static_cast<Real>(std::sqrt(err2) / (std::sqrt(ref2) + TinyReal));
    e.linf = linf;
    e.amp_num = static_cast<Real>(amp_n / (amp_d + TinyReal));
    e.amp_an = amp * std::exp(-alpha * Pi * Pi * t / (L * L + TinyReal));
    e.energy = static_cast<Real>(E);
    return e;
}
} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
    {
        return 1;
    }
    const Real alpha = cli.k / (cli.rho * cli.cp + TinyReal);
    const Real x0 = 0.0;
    const Vecd center(Real(0.5) * cli.L, 0.0, 0.0);
    const Vecd halfsize(Real(0.5) * cli.L, cli.half_yz, cli.half_yz);
    const Real dp = cli.dp;
    const Real shell = cli.boundary_width_factor * dp;

    const BoundingBoxd bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));
    SPHSystem sph_system(bounds, dp);
    SolidBody body(sph_system, makeShared<OphelieTestGlassBoxShape>("Slab", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMatterMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(particles, cli.t0);
    registerOphelieThermalDiffusionAuxFields(particles, cli.k);
    UniquePtr<Inner<>> inner = makeUnique<Inner<>>(body);

    hostInitAnalytic(particles, x0, cli.L, cli.t_h, cli.t_c, cli.amp, alpha, Real(0), cli.t0);
    hostApplyEndDirichlet(particles, x0, cli.L, shell, cli.t_h, cli.t_c, cli.t0);

    InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace(
        *inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
        cli.pair_weight_regularization);
    StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK> apply_diff(body, cli.rho, cli.cp, cli.t0);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cells(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner(*inner);

    std::ofstream ofs(cli.out_csv);
    ofs << "t,l2_rel,linf,amp_num,amp_an,amp_rel,energy\n";

    const Real tau = cli.L * cli.L / (alpha * Pi * Pi + TinyReal);
    std::cout << "[doc77b] start n=" << particles.TotalRealParticles() << " L=" << cli.L << " dp=" << dp
              << " alpha=" << alpha << " tau_mode=" << tau << " T0=" << cli.t0 << " Th=" << cli.t_h
              << " Tc=" << cli.t_c << " A=" << cli.amp << " dt=" << cli.dt << " end=" << cli.end_time
              << " pair_weight_reg=" << cli.pair_weight_regularization << std::endl
              << std::flush;

    ErrStats e0 = hostErrors(particles, x0, cli.L, cli.t_h, cli.t_c, cli.amp, alpha, Real(0), cli.t0, cli.rho, cli.cp,
                             shell);
    ofs << "0," << e0.l2 << "," << e0.linf << "," << e0.amp_num << "," << e0.amp_an << ",0," << e0.energy << "\n";

    Real t = 0.0;
    ErrStats e_end = e0;
    Real amp_rel_max = 0.0;
    while (t + TinyReal < cli.end_time)
    {
        update_cells.exec();
        update_inner.exec();
        laplace.exec();
        apply_diff.exec(cli.dt);
        hostApplyEndDirichlet(particles, x0, cli.L, shell, cli.t_h, cli.t_c, cli.t0);
        t += cli.dt;
        e_end = hostErrors(particles, x0, cli.L, cli.t_h, cli.t_c, cli.amp, alpha, t, cli.t0, cli.rho, cli.cp, shell);
        const Real amp_rel =
            std::abs(e_end.amp_num - e_end.amp_an) / (std::abs(e_end.amp_an) + TinyReal);
        amp_rel_max = std::max(amp_rel_max, amp_rel);
        ofs << t << "," << e_end.l2 << "," << e_end.linf << "," << e_end.amp_num << "," << e_end.amp_an << ","
            << amp_rel << "," << e_end.energy << "\n";
        if (static_cast<int>(t / cli.dt + Real(0.5)) % 20 == 0)
        {
            std::cout << "[doc77b] t=" << t << " l2=" << e_end.l2 << " linf=" << e_end.linf
                      << " amp_num=" << e_end.amp_num << " amp_an=" << e_end.amp_an << " amp_rel=" << amp_rel
                      << std::endl
                      << std::flush;
        }
    }

    // Decay-rate check: log(amp) slope vs −α π²/L² over mid window.
    // Pre-declared: L2_rel < 0.08, |amp−amp_an|/|amp_an| < 0.15 at end, amp decays.
    const bool decay_ok = e_end.amp_num < e0.amp_num && e_end.amp_an > TinyReal;
    const Real amp_rel_end = std::abs(e_end.amp_num - e_end.amp_an) / (std::abs(e_end.amp_an) + TinyReal);
    const bool amp_ok = amp_rel_end < Real(0.15);
    const bool l2_ok = e_end.l2 < Real(0.08);
    const bool linf_ok = e_end.linf < Real(1.5); // K, on ~8 K mode + 1473 background
    const bool passed = decay_ok && amp_ok && l2_ok && linf_ok;

    std::cout << "[doc77b] end t=" << t << " l2=" << e_end.l2 << " linf=" << e_end.linf
              << " amp_num=" << e_end.amp_num << " amp_an=" << e_end.amp_an << " amp_rel_end=" << amp_rel_end
              << " amp_rel_max=" << amp_rel_max << " decay_ok=" << (decay_ok ? 1 : 0) << " amp_ok=" << (amp_ok ? 1 : 0)
              << " l2_ok=" << (l2_ok ? 1 : 0) << " linf_ok=" << (linf_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    std::cout << "[doc77b] reuse_path=OpheliePairwiseLaplaceCK+ApplyLaplaceDiffusionDtCK(ΔT-authoritative)"
              << std::endl;
    if (!cli.summary_json.empty())
    {
        std::ofstream js(cli.summary_json);
        js << "{\n"
           << "  \"case\": \"doc77b\",\n"
           << "  \"dp\": " << dp << ",\n"
           << "  \"pair_weight_regularization\": " << cli.pair_weight_regularization << ",\n"
           << "  \"end_time\": " << t << ",\n"
           << "  \"l2_rel\": " << e_end.l2 << ",\n"
           << "  \"linf\": " << e_end.linf << ",\n"
           << "  \"amp_num\": " << e_end.amp_num << ",\n"
           << "  \"amp_an\": " << e_end.amp_an << ",\n"
           << "  \"amp_rel_end\": " << amp_rel_end << ",\n"
           << "  \"decay_ok\": " << (decay_ok ? 1 : 0) << ",\n"
           << "  \"amp_ok\": " << (amp_ok ? 1 : 0) << ",\n"
           << "  \"l2_ok\": " << (l2_ok ? 1 : 0) << ",\n"
           << "  \"linf_ok\": " << (linf_ok ? 1 : 0) << ",\n"
           << "  \"passed\": " << (passed ? 1 : 0) << ",\n"
           << "  \"compare\": {\"l2_rel_target\": 0.08, \"amp_rel_target\": 0.15}\n"
           << "}\n";
        std::cout << "[doc77b] summary_json=" << cli.summary_json << std::endl;
    }
    return passed ? 0 : 2;
}
