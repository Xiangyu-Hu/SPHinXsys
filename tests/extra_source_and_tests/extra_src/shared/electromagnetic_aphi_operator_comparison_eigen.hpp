#ifndef ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_HPP
#define ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_HPP

#include "electromagnetic_aphi_operator_comparison_eigen.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace
{
inline SpMatC make_repeated_diagonal_sparse_complex(const StdVec<Real> &entries,
                                                    int repeats)
{
    std::vector<TripletC> triplets;
    triplets.reserve(entries.size() * static_cast<size_t>(repeats));
    int base_size = static_cast<int>(entries.size());
    for (int repeat = 0; repeat != repeats; ++repeat)
    {
        int offset = repeat * base_size;
        for (size_t i = 0; i != entries.size(); ++i)
        {
            triplets.emplace_back(offset + static_cast<int>(i),
                                  offset + static_cast<int>(i),
                                  Complex(entries[i], 0.0));
        }
    }
    SpMatC diagonal(base_size * repeats, base_size * repeats);
    diagonal.setFromTriplets(triplets.begin(), triplets.end());
    diagonal.makeCompressed();
    return diagonal;
}
} // namespace

//=================================================================================================//
DiscreteEMOperatorComparator::
    DiscreteEMOperatorComparator(BaseInnerRelation &inner_relation,
                                 const DiscreteOperatorComparisonParameters &parameters)
    : inner_relation_(inner_relation),
      parameters_(parameters),
      number_of_particles_(inner_relation.getSPHBody().getBaseParticles().TotalRealParticles()),
      reference_smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      volumetric_measure_(inner_relation.getSPHBody().getBaseParticles().getVariableDataByName<Real>("VolumetricMeasure")),
      grad_x_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      grad_y_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      grad_z_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      laplace_nu_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      laplace_operator_(static_cast<int>(3 * number_of_particles_), static_cast<int>(3 * number_of_particles_)),
      weak_curlcurl_operator_(static_cast<int>(3 * number_of_particles_), static_cast<int>(3 * number_of_particles_)),
      divergence_matrix_(static_cast<int>(number_of_particles_), static_cast<int>(3 * number_of_particles_)),
      curl_matrix_(static_cast<int>(3 * number_of_particles_), static_cast<int>(3 * number_of_particles_))
{
}
//=================================================================================================//
Real DiscreteEMOperatorComparator::harmonicMean(Real value_i, Real value_j)
{
    return 2.0 * value_i * value_j / (value_i + value_j + TinyReal);
}
//=================================================================================================//
Real DiscreteEMOperatorComparator::
    pairwiseDiffusionWeight(const Neighborhood &inner_neighborhood,
                            size_t neighbor_index) const
{
    Real dW_ijV_j = inner_neighborhood.dW_ij_[neighbor_index] *
                    volumetric_measure_[inner_neighborhood.j_[neighbor_index]];
    Real r_ij = inner_neighborhood.r_ij_[neighbor_index];
    return -2.0 * dW_ijV_j /
           (r_ij + parameters_.pair_weight_regularization * reference_smoothing_length_ + TinyReal);
}
//=================================================================================================//
void DiscreteEMOperatorComparator::
    buildOperators(const StdVec<Real> &magnetic_reluctivity)
{
    std::vector<TripletC> grad_x_triplets;
    std::vector<TripletC> grad_y_triplets;
    std::vector<TripletC> grad_z_triplets;
    std::vector<TripletC> laplace_nu_triplets;

    grad_x_triplets.reserve(number_of_particles_ * 16);
    grad_y_triplets.reserve(number_of_particles_ * 16);
    grad_z_triplets.reserve(number_of_particles_ * 16);
    laplace_nu_triplets.reserve(number_of_particles_ * 16);

    for (size_t index_i = 0; index_i != number_of_particles_; ++index_i)
    {
        Neighborhood &inner_neighborhood = inner_relation_.inner_configuration_[index_i];
        Real grad_x_diagonal = 0.0;
        Real grad_y_diagonal = 0.0;
        Real grad_z_diagonal = 0.0;
        Real laplace_diagonal = 0.0;

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (inner_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }

            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] *
                               volumetric_measure_[index_j] * inner_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }

            grad_x_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                         Complex(gradW_ijV_j[0], 0.0));
            grad_y_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                         Complex(gradW_ijV_j[1], 0.0));
            grad_z_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                         Complex(gradW_ijV_j[2], 0.0));
            grad_x_diagonal -= gradW_ijV_j[0];
            grad_y_diagonal -= gradW_ijV_j[1];
            grad_z_diagonal -= gradW_ijV_j[2];

            Real diffusion_weight = pairwiseDiffusionWeight(inner_neighborhood, n);
            Real nu_ij = harmonicMean(magnetic_reluctivity[index_i], magnetic_reluctivity[index_j]);
            laplace_nu_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                             Complex(-nu_ij * diffusion_weight, 0.0));
            laplace_diagonal += nu_ij * diffusion_weight;
        }

        grad_x_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_x_diagonal, 0.0));
        grad_y_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_y_diagonal, 0.0));
        grad_z_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_z_diagonal, 0.0));
        laplace_nu_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                         Complex(laplace_diagonal, 0.0));
    }

    grad_x_.setFromTriplets(grad_x_triplets.begin(), grad_x_triplets.end());
    grad_y_.setFromTriplets(grad_y_triplets.begin(), grad_y_triplets.end());
    grad_z_.setFromTriplets(grad_z_triplets.begin(), grad_z_triplets.end());
    laplace_nu_.setFromTriplets(laplace_nu_triplets.begin(), laplace_nu_triplets.end());
    grad_x_.makeCompressed();
    grad_y_.makeCompressed();
    grad_z_.makeCompressed();
    laplace_nu_.makeCompressed();

    std::vector<TripletC> divergence_triplets;
    divergence_triplets.reserve(number_of_particles_ * 12);
    for (int outer = 0; outer != grad_x_.outerSize(); ++outer)
    {
        for (SpMatC::InnerIterator it(grad_x_, outer); it; ++it)
        {
            divergence_triplets.emplace_back(it.row(), idxAx(static_cast<size_t>(it.col())), it.value());
        }
    }
    for (int outer = 0; outer != grad_y_.outerSize(); ++outer)
    {
        for (SpMatC::InnerIterator it(grad_y_, outer); it; ++it)
        {
            divergence_triplets.emplace_back(it.row(), idxAy(static_cast<size_t>(it.col())), it.value());
        }
    }
    for (int outer = 0; outer != grad_z_.outerSize(); ++outer)
    {
        for (SpMatC::InnerIterator it(grad_z_, outer); it; ++it)
        {
            divergence_triplets.emplace_back(it.row(), idxAz(static_cast<size_t>(it.col())), it.value());
        }
    }
    divergence_matrix_.setFromTriplets(divergence_triplets.begin(), divergence_triplets.end());
    divergence_matrix_.makeCompressed();

    std::vector<TripletC> curl_triplets;
    curl_triplets.reserve(number_of_particles_ * 18);
    auto append_curl_block = [&](const SpMatC &block, int row_offset, int col_offset, Complex scale)
    {
        for (int outer = 0; outer != block.outerSize(); ++outer)
        {
            for (SpMatC::InnerIterator it(block, outer); it; ++it)
            {
                curl_triplets.emplace_back(row_offset + it.row(),
                                           col_offset + it.col(),
                                           scale * it.value());
            }
        }
    };

    append_curl_block(grad_y_, 0, idxAz(0), Complex(1.0, 0.0));
    append_curl_block(grad_z_, 0, idxAy(0), Complex(-1.0, 0.0));
    append_curl_block(grad_z_, static_cast<int>(number_of_particles_), idxAx(0), Complex(1.0, 0.0));
    append_curl_block(grad_x_, static_cast<int>(number_of_particles_), idxAz(0), Complex(-1.0, 0.0));
    append_curl_block(grad_x_, static_cast<int>(2 * number_of_particles_), idxAy(0), Complex(1.0, 0.0));
    append_curl_block(grad_y_, static_cast<int>(2 * number_of_particles_), idxAx(0), Complex(-1.0, 0.0));
    curl_matrix_.setFromTriplets(curl_triplets.begin(), curl_triplets.end());
    curl_matrix_.makeCompressed();

    std::vector<TripletC> laplace_triplets;
    laplace_triplets.reserve(static_cast<size_t>(laplace_nu_.nonZeros() * 3));
    auto append_diagonal_block = [&](int offset)
    {
        for (int outer = 0; outer != laplace_nu_.outerSize(); ++outer)
        {
            for (SpMatC::InnerIterator it(laplace_nu_, outer); it; ++it)
            {
                laplace_triplets.emplace_back(offset + it.row(),
                                              offset + it.col(),
                                              it.value());
            }
        }
    };
    append_diagonal_block(0);
    append_diagonal_block(static_cast<int>(number_of_particles_));
    append_diagonal_block(static_cast<int>(2 * number_of_particles_));
    laplace_operator_.setFromTriplets(laplace_triplets.begin(), laplace_triplets.end());
    laplace_operator_.makeCompressed();

    StdVec<Real> nu_volume(number_of_particles_);
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        nu_volume[i] = magnetic_reluctivity[i] * volumetric_measure_[i];
    }
    SpMatC curl_mass = make_repeated_diagonal_sparse_complex(nu_volume, Dimensions);
    weak_curlcurl_operator_ = curl_matrix_.transpose() * curl_mass * curl_matrix_;

    if (parameters_.enable_gauge_penalty && parameters_.gauge_penalty > TinyReal)
    {
        SpMatC div_mass = make_repeated_diagonal_sparse_complex(nu_volume, 1);
        weak_curlcurl_operator_ +=
            Complex(parameters_.gauge_penalty, 0.0) *
            (divergence_matrix_.transpose() * div_mass * divergence_matrix_);
        laplace_operator_ +=
            Complex(parameters_.gauge_penalty, 0.0) *
            (divergence_matrix_.transpose() * div_mass * divergence_matrix_);
    }
}
//=================================================================================================//
VecC DiscreteEMOperatorComparator::
    gatherUnknownVector(const StdVec<Complex> &ax,
                        const StdVec<Complex> &ay,
                        const StdVec<Complex> &az) const
{
    VecC unknown_vector(static_cast<int>(3 * number_of_particles_));
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        unknown_vector[idxAx(i)] = ax[i];
        unknown_vector[idxAy(i)] = ay[i];
        unknown_vector[idxAz(i)] = az[i];
    }
    return unknown_vector;
}
//=================================================================================================//
DiscreteModeEnergySummary DiscreteEMOperatorComparator::
    evaluateMode(const StdVec<Complex> &ax,
                 const StdVec<Complex> &ay,
                 const StdVec<Complex> &az) const
{
    DiscreteModeEnergySummary summary;
    VecC unknown_vector = gatherUnknownVector(ax, ay, az);
    VecC curl_vector = curl_matrix_ * unknown_vector;
    VecC div_vector = divergence_matrix_ * unknown_vector;

    summary.laplace_energy = std::real(unknown_vector.dot(laplace_operator_ * unknown_vector));
    summary.weak_curlcurl_energy = std::real(unknown_vector.dot(weak_curlcurl_operator_ * unknown_vector));
    summary.divergence_l2 = div_vector.norm();
    summary.curl_l2 = curl_vector.norm();

    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        Real a_abs = std::sqrt(std::norm(ax[i]) + std::norm(ay[i]) + std::norm(az[i]));
        summary.max_abs_a = SMAX(summary.max_abs_a, a_abs);
    }
    return summary;
}
//=================================================================================================//

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_HPP
