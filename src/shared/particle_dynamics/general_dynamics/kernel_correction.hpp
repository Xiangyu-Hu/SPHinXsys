#ifndef KERNEL_CORRECTION_HPP
#define KERNEL_CORRECTION_HPP

#include "kernel_correction.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
LinearGradientCorrectionMatrix<DataDelegationType>::
    LinearGradientCorrectionMatrix(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      B_(this->particles_->template registerStateVariable<Matd>(
          "LinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
KernelGradientCorrection<DataDelegationType>::
    KernelGradientCorrection(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation){};
//=================================================================================================//
template <class DataDelegationType>
template <class PairAverageType>
void KernelGradientCorrection<DataDelegationType>::
    correctKernelGradient(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i)
{
    for (size_t n = 0; n != neighborhood.current_size_; ++n)
    {
        size_t index_j = neighborhood.j_[n];
        Vecd displacement = neighborhood.r_ij_[n] * neighborhood.e_ij_[n];

        Vecd corrected_direction = average_correction_matrix(index_i, index_j) * neighborhood.e_ij_[n];
        Real direction_norm = corrected_direction.norm();
        neighborhood.dW_ij_[n] *= direction_norm;
        neighborhood.e_ij_[n] = corrected_direction / (direction_norm + Eps);
        neighborhood.r_ij_[n] = displacement.dot(neighborhood.e_ij_[n]);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_CORRECTION_HPP
