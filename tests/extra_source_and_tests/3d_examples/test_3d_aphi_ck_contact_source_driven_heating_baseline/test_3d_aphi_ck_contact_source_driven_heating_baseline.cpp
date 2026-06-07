/** Stage 10.14 P4: two-body Contact source-driven baseline (air + TEAM7 right half, lambda_A off). */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h"

using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiContactSourceDrivenHeatingBaselineMetrics metrics = runTwoBodyContactSourceDrivenHeatingBaseline(ac, av);
    const bool passed = contactSourceDrivenHeatingBaselinePassed(metrics);
    printContactSourceDrivenHeatingBaselineMetrics("test_3d_aphi_ck_contact_source_driven_heating_baseline", metrics,
                                                   passed);
    return passed ? 0 : 1;
}
