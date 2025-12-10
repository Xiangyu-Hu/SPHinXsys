#include "level_set_transformation.h"

#include "mesh_iterators.hpp"

#ifndef LEVEL_SET_TRANSFORMATION_HPP
#define LEVEL_SET_TRANSFORMATION_HPP

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateLevelSetGradient::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : data_spacing_(encloser.index_handler_.DataSpacing()),
      phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &index)
                  {
                      phi_gradient_[package_index](index) =
                          regularizedCentralDifference(
                              phi_, cell_neighborhood_[package_index], index,
                              [](Real dp, Real dm)
                              {
                                  return 0.5 * (dp + dm); // no upwinding
                              }) /
                          data_spacing_;
                  });
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateKernelIntegrals::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : phi_(encloser.mv_phi_.DelegatedData(ex_policy)),
      phi_gradient_(encloser.mv_phi_gradient_.DelegatedData(ex_policy)),
      kernel_weight_(encloser.mv_kernel_weight_.DelegatedData(ex_policy)),
      kernel_gradient_(encloser.mv_kernel_gradient_.DelegatedData(ex_policy)),
      kernel_second_gradient_(encloser.mv_kernel_second_gradient_.DelegatedData(ex_policy)),
      kernel_(ex_policy, encloser.neighbor_method_),
      data_spacing_(encloser.index_handler_.DataSpacing()),
      data_cell_volume_(math::pow(data_spacing_, Dimensions)),
      cell_neighborhood_(encloser.dv_cell_neighborhood_.DelegatedData(ex_policy)),
      cutoff_radius_(encloser.neighbor_method_.CutOffRadius()),
      bounding_box_(BoundingBoxi(Arrayi::Constant(
          static_cast<int>(std::ceil((cutoff_radius_ - Eps) / data_spacing_))))) {}
//=================================================================================================//
inline void UpdateKernelIntegrals::UpdateKernel::update(const UnsignedInt &package_index)
{
    assignByDataIndex(
        kernel_weight_[package_index], [&](const Arrayi &data_index) -> Real
        { return computeKernelIntegral(package_index, data_index); });
    assignByDataIndex(
        kernel_gradient_[package_index], [&](const Arrayi &data_index) -> Vecd
        { return computeKernelGradientIntegral(package_index, data_index); });
    assignByDataIndex(
        kernel_second_gradient_[package_index], [&](const Arrayi &data_index) -> Matd
        { return computeKernelSecondGradientIntegral(package_index, data_index); });
}
//=================================================================================================//
inline Real UpdateKernelIntegrals::UpdateKernel::
    CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing)
{
    Real squared_norm_inv = 1.0 / (phi_gradient.squaredNorm() + TinyReal);
    Real volume_fraction(0);
    for (UnsignedInt i = 0; i != Dimensions; ++i)
    {
        volume_fraction += phi_gradient[i] * phi_gradient[i] * squared_norm_inv *
                           Heaviside(phi / (ABS(phi_gradient[i]) + TinyReal), 0.5 * data_spacing);
    }
    return volume_fraction;
}
//=================================================================================================//
template <typename DataType, typename FunctionByDataIndex>
DataType UpdateKernelIntegrals::UpdateKernel::
    computeIntegral(Real phi, const UnsignedInt &package_index, const Arrayi &grid_index,
                    const DataType &initial_value, const FunctionByDataIndex &function_by_index)
{
    DataType integral = initial_value;
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each(
            bounding_box_.lower_, bounding_box_.upper_ + Arrayi::Ones(),
            [&](const Arrayi &search_index)
            {
                DataPackagePair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + search_index);
                Real phi_neighbor = phi_[neighbor_meta.first](neighbor_meta.second);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first](neighbor_meta.second);
                    Vecd displacement = -search_index.cast<Real>().matrix() * data_spacing_;
                    if (displacement.norm() < cutoff_radius_)
                        integral += function_by_index(displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return integral * data_cell_volume_;
}
//=============================================================================================//
inline Real UpdateKernelIntegrals::UpdateKernel::
    computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    return phi > cutoff_radius_
               ? 1.0
               : computeIntegral(phi, package_index, data_index, 0.0,
                                 [&](const Vecd &displacement) -> Real
                                 { return kernel_.W(displacement); });
}
//=============================================================================================//
inline Vecd UpdateKernelIntegrals::UpdateKernel::
    computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    Vecd integral = Vecd::Zero();
    return computeIntegral(phi, package_index, data_index, integral,
                           [&](const Vecd &displacement) -> Vecd
                           { return kernel_.dW(displacement) * displacement /
                                    (displacement.norm() + TinyReal); });
}
//=============================================================================================//
inline Matd UpdateKernelIntegrals::UpdateKernel::
    computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &data_index)
{
    Real phi = phi_[package_index](data_index);
    Matd integral = Matd::Zero();
    return computeIntegral(phi, package_index, data_index, integral,
                           [&](const Vecd &displacement) -> Matd
                           { return kernel_.d2W(displacement) *
                                    displacement * displacement.transpose() /
                                    (displacement.squaredNorm() + TinyReal); });
}
//=================================================================================================//
} // namespace SPH
#endif // LEVEL_SET_TRANSFORMATION_HPP