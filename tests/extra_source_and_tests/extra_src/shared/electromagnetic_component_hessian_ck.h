#ifndef ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_H
#define ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_H

#include "base_general_dynamics.h"
#include <array>

namespace SPH
{
namespace electromagnetics
{
template <size_t Dimension = Vecd::RowsAtCompileTime>
using ComponentVariableNames = std::array<std::string, Dimension>;

/**
 * @brief Mirror a vector field into already-registered scalar component fields.
 * This is intended to reuse scalar CK gradient / Hessian operators for each component.
 */
class CopyVectorFieldComponentsToScalarVariables : public LocalDynamics
{
  public:
    explicit CopyVectorFieldComponentsToScalarVariables(
        SPHBody &sph_body,
        const std::string &vector_field_name,
        const ComponentVariableNames<> &component_variable_names);
    virtual ~CopyVectorFieldComponentsToScalarVariables() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    static constexpr size_t Components = Vecd::RowsAtCompileTime;
    Vecd *vector_field_;
    std::array<Real *, Components> component_variables_;
};

/**
 * @brief Reconstruct curl(nu curl(A)) from scalar component Hessians using
 *        curl curl(A) = grad(div A) - laplacian(A), assuming locally constant nu.
 * This is a validation / replacement path for the current CurlNuB discretization.
 */
class ReconstructCurlNuBFromScalarComponentHessians : public LocalDynamics
{
  public:
    explicit ReconstructCurlNuBFromScalarComponentHessians(
        SPHBody &sph_body,
        const ComponentVariableNames<> &component_variable_names,
        Real second_order_operator_scaling,
        const std::string &output_name = "CurlNuBFromScalarComponentHessians",
        const std::string &magnetic_reluctivity_name = "MagneticReluctivity");
    virtual ~ReconstructCurlNuBFromScalarComponentHessians() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    static constexpr size_t Components = Vecd::RowsAtCompileTime;
    Real second_order_operator_scaling_;
    Real *magnetic_reluctivity_;
    Vecd *curl_nu_b_;
    std::array<VecMatd *, Components> component_hessians_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_COMPONENT_HESSIAN_CK_H
