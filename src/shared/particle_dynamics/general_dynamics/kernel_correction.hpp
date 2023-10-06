#ifndef KERNEL_CORRECTION_HPP
#define KERNEL_CORRECTION_HPP

#include "kernel_correction.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType>
KernelCorrectionMatrix<DataDelegationType>::
    KernelCorrectionMatrix(BaseRelationType &base_relation, , Real alpha)
    : LocalDynamics(base_relation.getSPHBody()),
      DataDelegationType(base_relation),
      alpha_(alpha), B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class KernelCorrectionMatrixType>
KernelGradientCorrection<DataDelegationType>::
    KernelGradientCorrection(KernelCorrectionMatrixType &kernel_correction)
    : LocalDynamics(kernel_correction.getSPHBody()),
      DataDelegationType(kernel_correction.getRelation()),
      average_correction_matrix_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")){};
//=================================================================================================//
template <class DataDelegationType>
template <class PairAverageType>
void KernelGradientCorrection<DataDelegationType>::
    average_correction_matrix(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i)
{
    for (size_t n = 0; n != neighborhood.current_size_; ++n)
    {
        size_t index_j = neighborhood.j_[n];
        Vecd displacement = neighborhood.r_ij_[n] * neighborhood.e_ij_[n];

        Vecd corrected_direction = average_correction_matrix(index_i, index_j) * neighborhood.e_ij_[n];
        Real direction_norm = corrected_direction.norm();
        neighborhood.dW_ijV_j_[n] *= direction_norm;
        neighborhood.e_ij_[n] = corrected_direction / (direction_norm + Eps);
        neighborhood.r_ij_[n] = displacement.dot(neighborhood.e_ij_[n]);
    }
};
} // namespace SPH
#endif // KERNEL_CORRECTION_HPP
