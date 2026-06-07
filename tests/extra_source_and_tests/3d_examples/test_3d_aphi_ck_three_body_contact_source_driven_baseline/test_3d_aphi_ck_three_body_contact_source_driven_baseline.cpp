/** Stage 10.14 P4b: three-body Contact source-driven baseline. */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h"

using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiContactSourceDrivenHeatingBaselineMetrics metrics =
        runThreeBodyContactSourceDrivenHeatingBaseline(ac, av);
    const bool passed = contactSourceDrivenHeatingBaselinePassed(metrics);
    printContactSourceDrivenHeatingBaselineMetrics("test_3d_aphi_ck_three_body_contact_source_driven_baseline", metrics,
                                                   passed);
    return passed ? 0 : 1;
}
