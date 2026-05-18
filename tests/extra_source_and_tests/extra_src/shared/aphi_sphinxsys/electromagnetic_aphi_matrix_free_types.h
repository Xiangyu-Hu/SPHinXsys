#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_TYPES_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_TYPES_H

#include "sphinxsys.h"

#include <complex>
#include <Eigen/Core>

namespace SPH
{
namespace electromagnetics
{
using Complex = std::complex<Real>;
using Vec3c = Eigen::Matrix<Complex, Dimensions, 1>;

/** @brief POD complex scalar for device / SYCL kernels (avoid std::complex in kernel payloads). */
struct ComplexValue
{
    Real real_;
    Real imag_;
};

/** @brief POD complex 3-vector for device / SYCL kernels. */
struct ComplexVec3
{
    ComplexValue x_;
    ComplexValue y_;
    ComplexValue z_;
};

inline ComplexValue complexValueFromStd(const Complex &z)
{
    return ComplexValue{z.real(), z.imag()};
}

inline Complex complexToStd(const ComplexValue &v)
{
    return Complex(v.real_, v.imag_);
}

/** @brief One pseudo-iteration diagnostic row (host-side history; safe to log to CSV). */
struct MatrixFreeAPhiIterationRecord
{
    size_t outer_iteration_ = 0;
    size_t load_step_ = 0;
    Real effective_source_scale_ = 1.0;
    Real effective_load_scale_ = 1.0;
    Real field_update_l2_ = 0.0;
    Real relative_field_update_l2_ = 0.0;
    Real residual_ax_l2_ = 0.0;
    Real residual_ay_l2_ = 0.0;
    Real residual_az_l2_ = 0.0;
    Real residual_phi_l2_ = 0.0;
    Real divergence_a_l2_ = 0.0;
    Real divergence_j_l2_ = 0.0;
};

using MatrixFreeAPhiIterationHistory = StdVec<MatrixFreeAPhiIterationRecord>;

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_TYPES_H
