/** Stage 10.17: native TEAM7 reload Bz probe vs muFEM reference (50 Hz, SYCL CK). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_bz_reference_probe_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_geometry_bz_reference_probe";
    AphiTeam7NativeBzReferenceProbeSpec spec;
    spec.solve_spec.omega = 2.0 * Pi * 50.0;
    spec.min_max_abs_bz_a1b1_mT = 0.05;
    spec.require_bz_magnitude_sanity = true;
    /** P4: freeze conductor TEAM7 sign; must match prior best scan on A1-B1. */
    spec.use_frozen_team7_bz_sign = true;
    spec.frozen_sign_bz_real = 1;
    spec.frozen_sign_bz_imag = -1;
    spec.run_sign_scan_diagnostic = true;
    spec.require_frozen_sign_matches_scan = true;
    spec.max_a1b1_profile_rel_err = 0.25;
    spec.write_vtp = team7NativeCliRequestsWriteVtp(ac, av);
    const AphiTeam7NativeBzReferenceProbeSummary summary = runTeam7NativeBzReferenceProbe(ac, av, spec);
    const bool passed = team7NativeBzReferenceProbePassed(summary, spec);
    printTeam7NativeBzReferenceProbeSummary(k_test_name, summary, spec, passed);
    std::cout << k_test_name << " native_reload_geometry=1 stage_10_17_p4_frozen_sign=1" << std::endl;
    return passed ? 0 : 1;
}
