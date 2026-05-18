#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct LineGaugeProjectionDiagnostics
{
    Real chi_l2_error_;
    Real chi_max_error_;
    Real projection_before_l2_;
    Real projection_after_l2_;
    Real gradient_div_a_before_l2_;
    Real gradient_div_a_after_l2_;
    Real electric_field_before_l2_;
    Real electric_field_after_l2_;
    Real electric_field_change_l2_;

    LineGaugeProjectionDiagnostics();
};

StdVec<Complex> computeLineGradient(const StdVec<Complex> &field, Real dx,
                                    Complex left_boundary = Complex(0.0, 0.0),
                                    Complex right_boundary = Complex(0.0, 0.0));

StdVec<Complex> computeLineNegativeLaplace(const StdVec<Complex> &field, Real dx,
                                           Complex left_boundary = Complex(0.0, 0.0),
                                           Complex right_boundary = Complex(0.0, 0.0));

StdVec<Complex> computeLineElectricField(const StdVec<Complex> &ax, const StdVec<Complex> &phi, Real dx,
                                         Real angular_frequency,
                                         Complex left_phi_boundary = Complex(0.0, 0.0),
                                         Complex right_phi_boundary = Complex(0.0, 0.0));

void applyLineGaugeTransform(StdVec<Complex> &ax, StdVec<Complex> &phi, const StdVec<Complex> &chi, Real dx,
                             Real angular_frequency,
                             Complex left_chi_boundary = Complex(0.0, 0.0),
                             Complex right_chi_boundary = Complex(0.0, 0.0));

Real computeComplexFieldL2Norm(const StdVec<Complex> &field);

Real computeComplexFieldL2Difference(const StdVec<Complex> &lhs, const StdVec<Complex> &rhs);

LineGaugeProjectionDiagnostics evaluateLineGaugeProjectionDiagnostics(
    const StdVec<Complex> &exact_chi, const StdVec<Complex> &solved_chi, const StdVec<Complex> &projection_before,
    const StdVec<Complex> &projection_after, const StdVec<Complex> &gradient_div_a_before,
    const StdVec<Complex> &gradient_div_a_after, const StdVec<Complex> &electric_field_before,
    const StdVec<Complex> &electric_field_after);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_gauge_projection.hpp"

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_H
