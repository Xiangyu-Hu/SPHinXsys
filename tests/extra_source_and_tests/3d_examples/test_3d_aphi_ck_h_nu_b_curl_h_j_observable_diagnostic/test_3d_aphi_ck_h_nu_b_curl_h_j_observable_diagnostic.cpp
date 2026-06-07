/**
 * Stage 10.12-P5: H=nu B and curlH vs J observable diagnostic (Linear2D + Az2D, B-corrected grad).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_h_curl_h_j_observable_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real core_shell = 2.5 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real max_linear2d_b_error = 0.02;
    const Real max_linear2d_h_error = 0.02;
    const Real max_az2d_b_error = 0.08;
    const Real max_az2d_h_error = 0.08;

    const AphiHCurlHJObservableRow linear2d = runHCurlHJObservableRow(
        ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega,
        AphiDivFreeValidationFieldKind::Linear2D);
    const AphiHCurlHJObservableRow az2d = runHCurlHJObservableRow(
        ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega,
        AphiDivFreeValidationFieldKind::Az2D);

    const char *test_name = "test_3d_aphi_ck_h_nu_b_curl_h_j_observable_diagnostic";
    printHCurlHJObservableRow(test_name, linear2d);
    printHCurlHJObservableRow(test_name, az2d);

    const bool linear2d_ok = hCurlHJObservableLinear2DGateOk(linear2d, max_linear2d_b_error, max_linear2d_h_error);
    const bool az2d_ok = hCurlHJObservableOscillatoryGateOk(az2d, max_az2d_b_error, max_az2d_h_error);
    const bool curl_h_j_reported =
        std::isfinite(linear2d.curl_h_error_vs_exact) && std::isfinite(linear2d.j_error_vs_exact) &&
        std::isfinite(az2d.curl_h_error_vs_exact) && std::isfinite(az2d.j_error_vs_exact);
    const bool passed = linear2d_ok && az2d_ok && curl_h_j_reported;

    std::cout << test_name << " linear2d_ok=" << (linear2d_ok ? 1 : 0) << " az2d_ok=" << (az2d_ok ? 1 : 0)
              << " curl_h_j_reported=" << (curl_h_j_reported ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
