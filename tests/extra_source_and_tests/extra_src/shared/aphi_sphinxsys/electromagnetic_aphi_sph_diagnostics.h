#ifndef ELECTROMAGNETIC_APHI_SPH_DIAGNOSTICS_H
#define ELECTROMAGNETIC_APHI_SPH_DIAGNOSTICS_H

#include "execution_policy.h"

#include <type_traits>

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct MatrixFreeAPhiSphNativeDiagnostics
{
    bool uses_device_policy_ = std::is_same<execution::MainExecutionPolicy, execution::ParallelDevicePolicy>::value;

    size_t standard_gradient_calls_ = 0;
    size_t laplace_helmholtz_calls_ = 0;
    size_t divergence_of_gradient_calls_ = 0;
    size_t sigma_grad_phi_calls_ = 0;
    size_t harmonic_weighted_gradient_calls_ = 0;
    size_t standard_divergence_calls_ = 0;
    size_t harmonic_divergence_calls_ = 0;

    double standard_gradient_seconds_ = 0.0;
    double laplace_helmholtz_seconds_ = 0.0;
    double divergence_of_gradient_seconds_ = 0.0;
    double sigma_grad_phi_seconds_ = 0.0;
    double harmonic_weighted_gradient_seconds_ = 0.0;
    double standard_divergence_seconds_ = 0.0;
    double harmonic_divergence_seconds_ = 0.0;
};

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_DIAGNOSTICS_H
