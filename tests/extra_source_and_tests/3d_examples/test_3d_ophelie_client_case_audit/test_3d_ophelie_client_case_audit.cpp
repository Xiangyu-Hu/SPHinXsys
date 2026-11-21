/**
 * @file test_3d_ophelie_client_case_audit.cpp
 * @brief Host audit for Doc229 client materials / geometry / effective presets.
 *
 * Does NOT claim NC temperature PASS or EM absorbed-power prediction.
 * Writes formula CSV + assembly JSON + effective-case JSON under --output-dir.
 */
#include "electromagnetic_ophelie_client_case_params.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{

struct Cli
{
    std::string output_dir = "client_case_audit_out";
    std::string repo_root;
    bool fail_on_missing_stl = true;
};

bool parseCli(int argc, char **argv, Cli &cli)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--output-dir") == 0 && i + 1 < argc)
            cli.output_dir = argv[++i];
        else if (std::strcmp(argv[i], "--repo-root") == 0 && i + 1 < argc)
            cli.repo_root = argv[++i];
        else if (std::strcmp(argv[i], "--allow-missing-stl") == 0)
            cli.fail_on_missing_stl = false;
        else if (std::strcmp(argv[i], "--help") == 0)
        {
            std::cout << "Usage: test_3d_ophelie_client_case_audit [--output-dir DIR] [--repo-root DIR]\n";
            return false;
        }
    }
    if (cli.repo_root.empty())
    {
        // Prefer walking up from CWD for client NC relax data/
        const fs::path marker =
            "tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/data/glass_NC.STL";
        fs::path cur = fs::current_path();
        for (int up = 0; up < 8; ++up)
        {
            if (fs::exists(cur / marker))
            {
                cli.repo_root = cur.string();
                break;
            }
            if (!cur.has_parent_path() || cur == cur.root_path())
                break;
            cur = cur.parent_path();
        }
    }
    return true;
}

std::string sha256FileHex(const fs::path &path)
{
    // Portable enough: shell out to sha256sum when available; else size+mtime token.
    const std::string cmd = "sha256sum \"" + path.string() + "\" 2>/dev/null | awk '{print $1}'";
    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        return "unavailable";
    char buf[128]{};
    std::string out;
    while (fgets(buf, sizeof(buf), pipe))
        out += buf;
    pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r' || out.back() == ' '))
        out.pop_back();
    if (out.size() == 64)
        return out;
    std::ostringstream oss;
    oss << "nostdhash_size=" << fs::file_size(path);
    return oss.str();
}

bool auditAssemblyStls(const OphelieClientAssemblyCatalog &cat, const fs::path &repo, const fs::path &out_dir,
                       bool fail_on_missing, std::ostream &log)
{
    const fs::path sub = repo / cat.stl_subdir;
    std::ofstream man(out_dir / (std::string(cat.assembly_id) + "_stl_manifest.csv"));
    man << "role,file,exists,bytes,sha256_or_token\n";
    bool ok = true;
    for (const auto &part : cat.parts)
    {
        const fs::path f = sub / part.relative_stl;
        const bool exists = fs::exists(f);
        std::uint64_t bytes = 0;
        std::string hash = "missing";
        if (exists)
        {
            bytes = static_cast<std::uint64_t>(fs::file_size(f));
            hash = sha256FileHex(f);
        }
        else
        {
            ok = false;
            log << "[FAIL] missing STL " << f << "\n";
        }
        man << part.role << ',' << part.relative_stl << ',' << (exists ? 1 : 0) << ',' << bytes << ',' << hash
            << '\n';
    }
    if (!writeClientAssemblyMetadataJson(cat, (out_dir / (std::string(cat.assembly_id) + "_geometry.json")).string()))
    {
        log << "[FAIL] cannot write geometry json for " << cat.assembly_id << "\n";
        ok = false;
    }
    if (!ok && !fail_on_missing)
    {
        log << "[WARN] STL gaps allowed by --allow-missing-stl\n";
        return true;
    }
    return ok;
}

} // namespace

int main(int argc, char *argv[])
{
    Cli cli;
    if (!parseCli(argc, argv, cli))
        return 2;

    fs::create_directories(cli.output_dir);
    std::ofstream summary(fs::path(cli.output_dir) / "AUDIT_SUMMARY.md");
    summary << "# Client case audit (Doc229)\n\n";
    summary << "- repo_root: `" << (cli.repo_root.empty() ? "(not found)" : cli.repo_root) << "`\n";
    summary << "- output: `" << cli.output_dir << "`\n\n";

    bool ok = true;

    // 1) Formula recompute
    {
        const fs::path csv = fs::path(cli.output_dir) / "client_formula_check.csv";
        std::string err;
        const bool formula_ok = writeClientFormulaCheckCsv(csv.string(), &err);
        if (!formula_ok)
        {
            std::cerr << "[FAIL] formula check " << err << "\n";
            summary << "- formula_check: **FAIL**\n";
            ok = false;
        }
        else
        {
            std::cout << "[PASS] formula check -> " << csv << "\n";
            summary << "- formula_check: PASS (`client_formula_check.csv`)\n";
            // Explicit 1473 K gate from Doc229
            const Real mu = clientDynamicViscosityPaS(Real(1473));
            const Real sg = clientConductivityLnExpSm(Real(1473));
            summary << "  - 1473 K: mu=" << mu << " (ref 5.87400286), sigma=" << sg << " (ref 0.28068105)\n";
            summary << "  - **retired**: pow(10, 4.5-8500/T) client sigma reading\n";
        }
    }

    if (cli.repo_root.empty())
    {
        std::cerr << "[FAIL] cannot locate repo root with client NC relax data/ STLs; pass --repo-root\n";
        summary << "- geometry: **FAIL** (repo root not found)\n";
        ok = false;
    }
    else
    {
        const fs::path repo = cli.repo_root;
        const auto nc = makeClientNaturalConvectionAssemblyCatalog();
        const auto st = makeClientStirringAssemblyCatalog();
        ok = auditAssemblyStls(nc, repo, cli.output_dir, cli.fail_on_missing_stl, std::cerr) && ok;
        ok = auditAssemblyStls(st, repo, cli.output_dir, cli.fail_on_missing_stl, std::cerr) && ok;
        summary << "- geometry NC/stir manifests: " << (ok ? "PASS/WARN see CSV" : "**FAIL**") << "\n";

        const auto probe_nc = makeClientNominalThermocoupleProbe(nc);
        summary << "- NC nominal probe solver_m = [" << probe_nc.solver_m.x << ", " << probe_nc.solver_m.y << ", "
                << probe_nc.solver_m.z << "] status=" << probe_nc.status << "\n";
    }

    // 2) Effective presets
    {
        const auto nc16 = makeClientNaturalConvectionDevelopmentParams();
        const auto nc_client = makeClientNaturalConvectionClientSigmaContrastParams();
        const auto stir = makeClientStirringDevelopmentParams();
        ok = writeClientEffectiveCaseJson(nc16, (fs::path(cli.output_dir) / "effective_client_nc_sigma16.json").string()) &&
             ok;
        ok = writeClientEffectiveCaseJson(
                 nc_client, (fs::path(cli.output_dir) / "effective_client_nc_client_ln_sigma.json").string()) &&
             ok;
        ok = writeClientEffectiveCaseJson(stir, (fs::path(cli.output_dir) / "effective_client_stir.json").string()) &&
             ok;

        const auto laws16 = buildClientCaseMaterialLaws(nc16);
        const auto laws_c = buildClientCaseMaterialLaws(nc_client);
        logClientMaterialIdentity(nc16.material_id, laws16, nc16.mu_freeze_t_k);
        logClientMaterialIdentity(nc_client.material_id, laws_c, nc_client.mu_freeze_t_k);

        summary << "- effective presets written: NC sigma16, NC client_ln_sigma, stirring\n";
        summary << "  - frozen mu @ " << nc16.mu_freeze_t_k << " K = " << laws16.mu.evaluate(nc16.mu_freeze_t_k)
                << " Pa·s\n";
        summary << "  - sigma16 @ ref = " << laws16.sigma.evaluate(nc16.mu_freeze_t_k)
                << " ; client_ln @ ref = " << laws_c.sigma.evaluate(nc_client.mu_freeze_t_k) << "\n";
        summary << "  - generator_input_power_W=30000 (not absorbed)\n";
        summary << "  - coil_current: **unknown** (no invented kA)\n";
    }

    summary << "\n## Explicitly NOT claimed\n\n";
    summary << "- NC thermocouple 1000°C quantitative PASS\n";
    summary << "- EM absorbed power equal to 30 kW input\n";
    summary << "- Doc228 elastic paddle E0–F1\n";
    summary << "- same-current sigma16 vs client_ln EM Q compare (needs confirmed current + EM run)\n";

    std::cout << (ok ? "[PASS] client case audit\n" : "[FAIL] client case audit\n");
    return ok ? 0 : 1;
}
