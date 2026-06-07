/** Stage 10.17 P2: vacuum plate (sigma=0) + TEAM7 coil path source; Bz must enter mT range. */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_bz_reference_probe_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_vacuum_source_b_sanity";
    Team7CoilPathSourceSpec coil_spec;
    const AphiTeam7NativeSourceRhsAuditSummary source_audit = runTeam7NativeSourceRhsAudit(ac, av, coil_spec);
    const bool source_audit_ok = team7NativeSourceRhsAuditPassed(source_audit);

    AphiTeam7NativeBzReferenceProbeSpec spec;
    spec.report_dir = "team7_native_vacuum_source_b_sanity_report";
    spec.solve_spec.plate_sigma = 0.0;
    spec.solve_spec.vacuum_plate_nonconductive = true;
    spec.solve_spec.omega = 2.0 * Pi * 50.0;
    spec.min_max_abs_bz_a1b1_mT = 0.01;
    spec.max_max_abs_bz_a1b1_mT = 100.0;
    spec.require_a1b1_profile_gate = false;
    spec.require_bz_magnitude_sanity = true;

    const AphiTeam7NativeBzReferenceProbeSummary summary = runTeam7NativeBzReferenceProbe(ac, av, spec);
    const bool bz_ok = team7NativeBzReferenceProbePassed(summary, spec);
    const bool passed = source_audit_ok && bz_ok;

    printTeam7NativeSourceRhsAudit(k_test_name, source_audit, source_audit_ok);
    printTeam7NativeBzReferenceProbeSummary(k_test_name, summary, spec, bz_ok);
    std::cout << k_test_name << " final_passed=" << (passed ? 1 : 0)
              << " vacuum_plate_sigma=0 profile_gate=off (magnitude_only)" << std::endl;
    return passed ? 0 : 1;
}
