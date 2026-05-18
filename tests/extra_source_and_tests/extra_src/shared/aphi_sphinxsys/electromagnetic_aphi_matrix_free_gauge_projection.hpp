#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_gauge_projection.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline LineGaugeProjectionDiagnostics::LineGaugeProjectionDiagnostics()
    : chi_l2_error_(0.0), chi_max_error_(0.0), projection_before_l2_(0.0), projection_after_l2_(0.0),
      gradient_div_a_before_l2_(0.0), gradient_div_a_after_l2_(0.0), electric_field_before_l2_(0.0),
      electric_field_after_l2_(0.0), electric_field_change_l2_(0.0)
{
}

inline StdVec<Complex> computeLineGradient(const StdVec<Complex> &field, Real dx, Complex left_boundary,
                                           Complex right_boundary)
{
    const size_t number_of_points = field.size();
    StdVec<Complex> gradient(number_of_points, Complex(0.0, 0.0));

    if (number_of_points == 0)
    {
        return gradient;
    }

    const Real inverse_two_dx = 0.5 / dx;
    for (size_t i = 0; i != number_of_points; ++i)
    {
        const Complex left_value = (i == 0) ? left_boundary : field[i - 1];
        const Complex right_value = (i + 1 == number_of_points) ? right_boundary : field[i + 1];
        gradient[i] = (right_value - left_value) * inverse_two_dx;
    }

    return gradient;
}

inline StdVec<Complex> computeLineNegativeLaplace(const StdVec<Complex> &field, Real dx, Complex left_boundary,
                                                  Complex right_boundary)
{
    const size_t number_of_points = field.size();
    StdVec<Complex> laplace(number_of_points, Complex(0.0, 0.0));

    if (number_of_points == 0)
    {
        return laplace;
    }

    const Real inverse_dx_sq = 1.0 / (dx * dx);
    for (size_t i = 0; i != number_of_points; ++i)
    {
        const Complex left_value = (i == 0) ? left_boundary : field[i - 1];
        const Complex right_value = (i + 1 == number_of_points) ? right_boundary : field[i + 1];
        laplace[i] = (Complex(2.0, 0.0) * field[i] - left_value - right_value) * inverse_dx_sq;
    }

    return laplace;
}

inline StdVec<Complex> computeLineElectricField(const StdVec<Complex> &ax, const StdVec<Complex> &phi, Real dx,
                                                Real angular_frequency, Complex left_phi_boundary,
                                                Complex right_phi_boundary)
{
    const Complex imaginary_unit(0.0, 1.0);
    const StdVec<Complex> grad_phi = computeLineGradient(phi, dx, left_phi_boundary, right_phi_boundary);
    StdVec<Complex> electric_field(ax.size(), Complex(0.0, 0.0));

    for (size_t i = 0; i != ax.size(); ++i)
    {
        electric_field[i] = -imaginary_unit * angular_frequency * ax[i] - grad_phi[i];
    }

    return electric_field;
}

inline void applyLineGaugeTransform(StdVec<Complex> &ax, StdVec<Complex> &phi, const StdVec<Complex> &chi, Real dx,
                                    Real angular_frequency, Complex left_chi_boundary, Complex right_chi_boundary)
{
    const Complex imaginary_unit(0.0, 1.0);
    const StdVec<Complex> grad_chi = computeLineGradient(chi, dx, left_chi_boundary, right_chi_boundary);

    for (size_t i = 0; i != ax.size(); ++i)
    {
        ax[i] -= grad_chi[i];
        phi[i] += imaginary_unit * angular_frequency * chi[i];
    }
}

inline Real computeComplexFieldL2Norm(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (const Complex &value : field)
    {
        squared_sum += std::norm(value);
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

inline Real computeComplexFieldL2Difference(const StdVec<Complex> &lhs, const StdVec<Complex> &rhs)
{
    if (lhs.size() != rhs.size() || lhs.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (size_t i = 0; i != lhs.size(); ++i)
    {
        squared_sum += std::norm(lhs[i] - rhs[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(lhs.size()));
}

inline LineGaugeProjectionDiagnostics evaluateLineGaugeProjectionDiagnostics(
    const StdVec<Complex> &exact_chi, const StdVec<Complex> &solved_chi, const StdVec<Complex> &projection_before,
    const StdVec<Complex> &projection_after, const StdVec<Complex> &gradient_div_a_before,
    const StdVec<Complex> &gradient_div_a_after, const StdVec<Complex> &electric_field_before,
    const StdVec<Complex> &electric_field_after)
{
    LineGaugeProjectionDiagnostics diagnostics;
    diagnostics.projection_before_l2_ = computeComplexFieldL2Norm(projection_before);
    diagnostics.projection_after_l2_ = computeComplexFieldL2Norm(projection_after);
    diagnostics.gradient_div_a_before_l2_ = computeComplexFieldL2Norm(gradient_div_a_before);
    diagnostics.gradient_div_a_after_l2_ = computeComplexFieldL2Norm(gradient_div_a_after);
    diagnostics.electric_field_before_l2_ = computeComplexFieldL2Norm(electric_field_before);
    diagnostics.electric_field_after_l2_ = computeComplexFieldL2Norm(electric_field_after);
    diagnostics.electric_field_change_l2_ =
        computeComplexFieldL2Difference(electric_field_after, electric_field_before);
    diagnostics.chi_l2_error_ = computeComplexFieldL2Difference(solved_chi, exact_chi);

    Real max_error = 0.0;
    for (size_t i = 0; i != solved_chi.size() && i != exact_chi.size(); ++i)
    {
        max_error = std::max(max_error, std::abs(solved_chi[i] - exact_chi[i]));
    }
    diagnostics.chi_max_error_ = max_error;

    return diagnostics;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_GAUGE_PROJECTION_HPP
