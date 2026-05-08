#ifndef ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_HPP
#define ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_HPP

#include "electromagnetic_aphi_laplace_eigen.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace
{
inline SpMatC make_diagonal_sparse_complex(const StdVec<Real> &entries)
{
    std::vector<TripletC> triplets;
    triplets.reserve(entries.size());
    for (size_t i = 0; i != entries.size(); ++i)
    {
        triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), Complex(entries[i], 0.0));
    }
    SpMatC diagonal(static_cast<int>(entries.size()), static_cast<int>(entries.size()));
    diagonal.setFromTriplets(triplets.begin(), triplets.end());
    diagonal.makeCompressed();
    return diagonal;
}

inline void overwrite_row_with_identity(SpMatC &matrix,
                                        VecC &right_hand_side,
                                        int row_index,
                                        const Complex &value)
{
    for (int column = 0; column != matrix.outerSize(); ++column)
    {
        for (SpMatC::InnerIterator iterator(matrix, column); iterator; ++iterator)
        {
            if (iterator.row() == row_index)
            {
                iterator.valueRef() = Complex(0.0, 0.0);
            }
        }
    }
    matrix.coeffRef(row_index, row_index) = Complex(1.0, 0.0);
    right_hand_side[row_index] = value;
}
} // namespace

//=================================================================================================//
LaplaceStructuredAPhiEigenSolver::
    LaplaceStructuredAPhiEigenSolver(BaseInnerRelation &inner_relation,
                                     const LaplaceStructuredAPhiParameters &parameters)
    : LaplaceStructuredAPhiEigenSolver(
          LaplaceStructuredAPhiDiscreteView{
              inner_relation.getSPHBody().getBaseParticles().TotalRealParticles(),
              inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength(),
              inner_relation.getSPHBody().getBaseParticles().getVariableDataByName<Vecd>("Position"),
              inner_relation.getSPHBody().getBaseParticles().getVariableDataByName<Real>("VolumetricMeasure"),
              nullptr,
              &inner_relation.inner_configuration_,
              nullptr},
          parameters)
{
}
//=================================================================================================//
LaplaceStructuredAPhiEigenSolver::
    LaplaceStructuredAPhiEigenSolver(const LaplaceStructuredAPhiDiscreteView &discrete_view,
                                     const LaplaceStructuredAPhiParameters &parameters)
    : parameters_(parameters),
      number_of_particles_(discrete_view.number_of_particles),
      reference_smoothing_length_(discrete_view.reference_smoothing_length),
      positions_(discrete_view.positions),
      volumetric_measure_(discrete_view.volumetric_measure),
      smoothing_length_ratio_(discrete_view.smoothing_length_ratio),
      particle_configuration_(discrete_view.particle_configuration),
      neighbor_is_contact_(discrete_view.neighbor_is_contact),
      grad_x_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      grad_y_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      grad_z_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      laplace_sigma_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      laplace_nu_(static_cast<int>(number_of_particles_), static_cast<int>(number_of_particles_)),
      system_matrix_(static_cast<int>(4 * number_of_particles_), static_cast<int>(4 * number_of_particles_)),
      right_hand_side_(static_cast<int>(4 * number_of_particles_))
{
    right_hand_side_.setZero();
    used_phi_block_scale_ = 1.0;
    used_diagonal_equilibration_iterations_ = 0;
    last_solver_backend_ = "not_run";
    last_solver_iterations_ = 0;
    last_solver_estimated_error_ = 0.0;
}
//=================================================================================================//
Real LaplaceStructuredAPhiEigenSolver::harmonicMean(Real value_i, Real value_j)
{
    return 2.0 * value_i * value_j / (value_i + value_j + TinyReal);
}
//=================================================================================================//
Real LaplaceStructuredAPhiEigenSolver::complexAbs(const Complex &value)
{
    return std::abs(value);
}
//=================================================================================================//
Real LaplaceStructuredAPhiEigenSolver::
    choosePhiBlockScale(Real angular_frequency, Real requested_scale)
{
    if (requested_scale > TinyReal)
    {
        return requested_scale;
    }
    return std::sqrt(SMAX(static_cast<Real>(1.0), std::abs(angular_frequency)));
}
//=================================================================================================//
Real LaplaceStructuredAPhiEigenSolver::
    pairwiseDiffusionWeight(const Neighborhood &inner_neighborhood,
                            size_t index_i,
                            size_t neighbor_index,
                            Real reference_smoothing_length) const
{
    Real dW_ijV_j = inner_neighborhood.dW_ij_[neighbor_index] *
                    volumetric_measure_[inner_neighborhood.j_[neighbor_index]];
    Real r_ij = inner_neighborhood.r_ij_[neighbor_index];
    size_t index_j = inner_neighborhood.j_[neighbor_index];
    Real smoothing_length_i = reference_smoothing_length;
    Real smoothing_length_j = reference_smoothing_length;
    if (smoothing_length_ratio_ != nullptr)
    {
        smoothing_length_i *= smoothing_length_ratio_[index_i];
        smoothing_length_j *= smoothing_length_ratio_[index_j];
    }
    Real pair_smoothing_length = SMAX(smoothing_length_i, smoothing_length_j);
    return -2.0 * dW_ijV_j /
           (r_ij + parameters_.pair_weight_regularization * pair_smoothing_length + TinyReal);
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::
    buildOperatorMatrices(const StdVec<Real> &electrical_conductivity,
                          const StdVec<Real> &magnetic_reluctivity)
{
    last_electrical_conductivity_ = electrical_conductivity;
    last_magnetic_reluctivity_ = magnetic_reluctivity;

    std::vector<TripletC> grad_x_triplets;
    std::vector<TripletC> grad_y_triplets;
    std::vector<TripletC> grad_z_triplets;
    std::vector<TripletC> laplace_sigma_triplets;
    std::vector<TripletC> laplace_nu_triplets;

    grad_x_triplets.reserve(number_of_particles_ * 16);
    grad_y_triplets.reserve(number_of_particles_ * 16);
    grad_z_triplets.reserve(number_of_particles_ * 16);
    laplace_sigma_triplets.reserve(number_of_particles_ * 16);
    laplace_nu_triplets.reserve(number_of_particles_ * 16);

    for (size_t index_i = 0; index_i != number_of_particles_; ++index_i)
    {
        const Neighborhood &inner_neighborhood = (*particle_configuration_)[index_i];

        Real grad_x_diagonal = 0.0;
        Real grad_y_diagonal = 0.0;
        Real grad_z_diagonal = 0.0;
        Real laplace_sigma_diagonal = 0.0;
        Real laplace_nu_diagonal = 0.0;

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (inner_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            bool is_contact_neighbor =
                neighbor_is_contact_ != nullptr &&
                index_i < neighbor_is_contact_->size() &&
                n < (*neighbor_is_contact_)[index_i].size() &&
                (*neighbor_is_contact_)[index_i][n] != 0;
            Real gradient_scale = is_contact_neighbor ? parameters_.contact_gradient_scale : Real(1.0);
            Real diffusion_scale = is_contact_neighbor ? parameters_.contact_diffusion_scale : Real(1.0);

            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] *
                               volumetric_measure_[index_j] * inner_neighborhood.e_ij_[n];
            gradW_ijV_j *= gradient_scale;
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

            Real diffusion_weight = pairwiseDiffusionWeight(inner_neighborhood, index_i, n, reference_smoothing_length_);
            diffusion_weight *= diffusion_scale;
            Real sigma_ij = harmonicMean(electrical_conductivity[index_i], electrical_conductivity[index_j]);
            Real nu_ij = harmonicMean(magnetic_reluctivity[index_i], magnetic_reluctivity[index_j]);
            Real laplace_sigma_off_diagonal = -sigma_ij * diffusion_weight;
            Real laplace_nu_off_diagonal = -nu_ij * diffusion_weight;

            laplace_sigma_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                                Complex(laplace_sigma_off_diagonal, 0.0));
            laplace_nu_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_j),
                                             Complex(laplace_nu_off_diagonal, 0.0));
            laplace_sigma_diagonal += sigma_ij * diffusion_weight;
            laplace_nu_diagonal += nu_ij * diffusion_weight;
        }

        grad_x_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_x_diagonal, 0.0));
        grad_y_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_y_diagonal, 0.0));
        grad_z_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                     Complex(grad_z_diagonal, 0.0));
        laplace_sigma_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                            Complex(laplace_sigma_diagonal, 0.0));
        laplace_nu_triplets.emplace_back(static_cast<int>(index_i), static_cast<int>(index_i),
                                         Complex(laplace_nu_diagonal, 0.0));
    }

    grad_x_.setFromTriplets(grad_x_triplets.begin(), grad_x_triplets.end());
    grad_y_.setFromTriplets(grad_y_triplets.begin(), grad_y_triplets.end());
    grad_z_.setFromTriplets(grad_z_triplets.begin(), grad_z_triplets.end());
    laplace_sigma_.setFromTriplets(laplace_sigma_triplets.begin(), laplace_sigma_triplets.end());
    laplace_nu_.setFromTriplets(laplace_nu_triplets.begin(), laplace_nu_triplets.end());

    grad_x_.makeCompressed();
    grad_y_.makeCompressed();
    grad_z_.makeCompressed();
    laplace_sigma_.makeCompressed();
    laplace_nu_.makeCompressed();
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::buildSystemMatrix()
{
    Real geom_length_to_m = SMAX(parameters_.geom_length_to_m, TinyReal);
    Real grad_scale = static_cast<Real>(1.0) / geom_length_to_m;
    Real laplace_scale = grad_scale * grad_scale;
    SpMatC mass_sigma = make_diagonal_sparse_complex(last_electrical_conductivity_);

    SpMatC block_ax_phi = Complex(grad_scale, 0.0) * (mass_sigma * grad_x_);
    SpMatC block_ay_phi = Complex(grad_scale, 0.0) * (mass_sigma * grad_y_);
    SpMatC block_az_phi = Complex(grad_scale, 0.0) * (mass_sigma * grad_z_);

    SpMatC block_phi_ax =
        Complex(0.0, parameters_.angular_frequency * grad_scale) * (grad_x_.transpose() * mass_sigma);
    SpMatC block_phi_ay =
        Complex(0.0, parameters_.angular_frequency * grad_scale) * (grad_y_.transpose() * mass_sigma);
    SpMatC block_phi_az =
        Complex(0.0, parameters_.angular_frequency * grad_scale) * (grad_z_.transpose() * mass_sigma);

    SpMatC block_phi_phi = Complex(laplace_scale, 0.0) * laplace_sigma_;

    if (parameters_.enable_gauge_penalty && parameters_.gauge_penalty > TinyReal)
    {
        SpMatC magnetic_diagonal = make_diagonal_sparse_complex(last_magnetic_reluctivity_);
        Complex penalty_scale(parameters_.gauge_penalty, 0.0);
        penalty_scale *= Complex(laplace_scale, 0.0);
        SpMatC gauge_xx = penalty_scale * (grad_x_.transpose() * magnetic_diagonal * grad_x_);
        SpMatC gauge_xy = penalty_scale * (grad_x_.transpose() * magnetic_diagonal * grad_y_);
        SpMatC gauge_xz = penalty_scale * (grad_x_.transpose() * magnetic_diagonal * grad_z_);
        SpMatC gauge_yx = penalty_scale * (grad_y_.transpose() * magnetic_diagonal * grad_x_);
        SpMatC gauge_yy = penalty_scale * (grad_y_.transpose() * magnetic_diagonal * grad_y_);
        SpMatC gauge_yz = penalty_scale * (grad_y_.transpose() * magnetic_diagonal * grad_z_);
        SpMatC gauge_zx = penalty_scale * (grad_z_.transpose() * magnetic_diagonal * grad_x_);
        SpMatC gauge_zy = penalty_scale * (grad_z_.transpose() * magnetic_diagonal * grad_y_);
        SpMatC gauge_zz = penalty_scale * (grad_z_.transpose() * magnetic_diagonal * grad_z_);

        std::vector<TripletC> triplets;
        triplets.reserve(static_cast<size_t>(laplace_nu_.nonZeros() * 16));

        auto append_block = [&](const SpMatC &block, int row_offset, int col_offset)
        {
            for (int outer = 0; outer != block.outerSize(); ++outer)
            {
                for (SpMatC::InnerIterator iterator(block, outer); iterator; ++iterator)
                {
                    triplets.emplace_back(row_offset + iterator.row(),
                                          col_offset + iterator.col(),
                                          iterator.value());
                }
            }
        };

        SpMatC diagonal_block =
            Complex(laplace_scale, 0.0) * laplace_nu_ +
            Complex(0.0, parameters_.angular_frequency) * mass_sigma;
        append_block(diagonal_block + gauge_xx, 0, 0);
        append_block(diagonal_block + gauge_yy, static_cast<int>(number_of_particles_),
                     static_cast<int>(number_of_particles_));
        append_block(diagonal_block + gauge_zz, static_cast<int>(2 * number_of_particles_),
                     static_cast<int>(2 * number_of_particles_));
        append_block(gauge_xy, 0, static_cast<int>(number_of_particles_));
        append_block(gauge_xz, 0, static_cast<int>(2 * number_of_particles_));
        append_block(gauge_yx, static_cast<int>(number_of_particles_), 0);
        append_block(gauge_yz, static_cast<int>(number_of_particles_),
                     static_cast<int>(2 * number_of_particles_));
        append_block(gauge_zx, static_cast<int>(2 * number_of_particles_), 0);
        append_block(gauge_zy, static_cast<int>(2 * number_of_particles_),
                     static_cast<int>(number_of_particles_));
        append_block(block_ax_phi, 0, static_cast<int>(3 * number_of_particles_));
        append_block(block_ay_phi, static_cast<int>(number_of_particles_),
                     static_cast<int>(3 * number_of_particles_));
        append_block(block_az_phi, static_cast<int>(2 * number_of_particles_),
                     static_cast<int>(3 * number_of_particles_));
        append_block(block_phi_ax, static_cast<int>(3 * number_of_particles_), 0);
        append_block(block_phi_ay, static_cast<int>(3 * number_of_particles_),
                     static_cast<int>(number_of_particles_));
        append_block(block_phi_az, static_cast<int>(3 * number_of_particles_),
                     static_cast<int>(2 * number_of_particles_));
        append_block(block_phi_phi, static_cast<int>(3 * number_of_particles_),
                     static_cast<int>(3 * number_of_particles_));

        system_matrix_.setZero();
        system_matrix_.setFromTriplets(triplets.begin(), triplets.end());
        system_matrix_.makeCompressed();
        return;
    }

    std::vector<TripletC> triplets;
    triplets.reserve(static_cast<size_t>(laplace_nu_.nonZeros() * 10));

    auto append_block = [&](const SpMatC &block, int row_offset, int col_offset)
    {
        for (int outer = 0; outer != block.outerSize(); ++outer)
        {
            for (SpMatC::InnerIterator iterator(block, outer); iterator; ++iterator)
            {
                triplets.emplace_back(row_offset + iterator.row(),
                                      col_offset + iterator.col(),
                                      iterator.value());
            }
        }
    };

    SpMatC diagonal_block =
        Complex(laplace_scale, 0.0) * laplace_nu_ +
        Complex(0.0, parameters_.angular_frequency) * mass_sigma;
    append_block(diagonal_block, 0, 0);
    append_block(diagonal_block, static_cast<int>(number_of_particles_),
                 static_cast<int>(number_of_particles_));
    append_block(diagonal_block, static_cast<int>(2 * number_of_particles_),
                 static_cast<int>(2 * number_of_particles_));
    append_block(block_ax_phi, 0, static_cast<int>(3 * number_of_particles_));
    append_block(block_ay_phi, static_cast<int>(number_of_particles_),
                 static_cast<int>(3 * number_of_particles_));
    append_block(block_az_phi, static_cast<int>(2 * number_of_particles_),
                 static_cast<int>(3 * number_of_particles_));
    append_block(block_phi_ax, static_cast<int>(3 * number_of_particles_), 0);
    append_block(block_phi_ay, static_cast<int>(3 * number_of_particles_),
                 static_cast<int>(number_of_particles_));
    append_block(block_phi_az, static_cast<int>(3 * number_of_particles_),
                 static_cast<int>(2 * number_of_particles_));
    append_block(block_phi_phi, static_cast<int>(3 * number_of_particles_),
                 static_cast<int>(3 * number_of_particles_));

    system_matrix_.setZero();
    system_matrix_.setFromTriplets(triplets.begin(), triplets.end());
    system_matrix_.makeCompressed();
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::initializeBlockScaling()
{
    left_block_scaling_.assign(4 * number_of_particles_, 1.0);
    right_block_scaling_.assign(4 * number_of_particles_, 1.0);
    used_phi_block_scale_ = 1.0;
    used_diagonal_equilibration_iterations_ = 0;

    if (!parameters_.enable_block_scaling)
    {
        return;
    }

    used_phi_block_scale_ =
        choosePhiBlockScale(parameters_.angular_frequency, parameters_.phi_block_scale);
    Real safe_phi_scale = SMAX(used_phi_block_scale_, TinyReal);
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        left_block_scaling_[static_cast<size_t>(idxPhi(i))] = 1.0 / safe_phi_scale;
        right_block_scaling_[static_cast<size_t>(idxPhi(i))] = safe_phi_scale;
    }
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::applyBlockScalingToSystem()
{
    if (!parameters_.enable_block_scaling)
    {
        return;
    }

    for (int outer = 0; outer != system_matrix_.outerSize(); ++outer)
    {
        for (SpMatC::InnerIterator iterator(system_matrix_, outer); iterator; ++iterator)
        {
            size_t row_index = static_cast<size_t>(iterator.row());
            size_t col_index = static_cast<size_t>(iterator.col());
            iterator.valueRef() *=
                left_block_scaling_[row_index] * right_block_scaling_[col_index];
        }
    }

    for (int i = 0; i != right_hand_side_.size(); ++i)
    {
        right_hand_side_[i] *= left_block_scaling_[static_cast<size_t>(i)];
    }
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::applyDiagonalEquilibrationToSystem()
{
    if (!parameters_.enable_diagonal_equilibration)
    {
        return;
    }

    int iterations = SMAX(0, parameters_.diagonal_equilibration_iterations);
    if (iterations == 0)
    {
        return;
    }

    const Real scaling_floor = static_cast<Real>(1.0e-12);
    const Real scaling_ceiling = static_cast<Real>(1.0e12);
    const int system_size = system_matrix_.rows();

    for (int iteration = 0; iteration != iterations; ++iteration)
    {
        StdVec<Real> row_max(static_cast<size_t>(system_size), 0.0);
        StdVec<Real> col_max(static_cast<size_t>(system_size), 0.0);

        for (int outer = 0; outer != system_matrix_.outerSize(); ++outer)
        {
            for (SpMatC::InnerIterator iterator(system_matrix_, outer); iterator; ++iterator)
            {
                Real abs_value = complexAbs(iterator.value());
                size_t row_index = static_cast<size_t>(iterator.row());
                size_t col_index = static_cast<size_t>(iterator.col());
                row_max[row_index] = SMAX(row_max[row_index], abs_value);
                col_max[col_index] = SMAX(col_max[col_index], abs_value);
            }
        }

        StdVec<Real> row_scale(static_cast<size_t>(system_size), 1.0);
        StdVec<Real> col_scale(static_cast<size_t>(system_size), 1.0);
        for (int i = 0; i != system_size; ++i)
        {
            Real safe_row_max = SMAX(row_max[static_cast<size_t>(i)], scaling_floor);
            Real safe_col_max = SMAX(col_max[static_cast<size_t>(i)], scaling_floor);
            Real row_scale_candidate =
                static_cast<Real>(1.0) / std::sqrt(safe_row_max);
            Real col_scale_candidate =
                static_cast<Real>(1.0) / std::sqrt(safe_col_max);
            row_scale[static_cast<size_t>(i)] =
                SMIN(scaling_ceiling, SMAX(scaling_floor, row_scale_candidate));
            col_scale[static_cast<size_t>(i)] =
                SMIN(scaling_ceiling, SMAX(scaling_floor, col_scale_candidate));
        }

        for (int outer = 0; outer != system_matrix_.outerSize(); ++outer)
        {
            for (SpMatC::InnerIterator iterator(system_matrix_, outer); iterator; ++iterator)
            {
                size_t row_index = static_cast<size_t>(iterator.row());
                size_t col_index = static_cast<size_t>(iterator.col());
                iterator.valueRef() *= row_scale[row_index] * col_scale[col_index];
            }
        }

        for (int i = 0; i != right_hand_side_.size(); ++i)
        {
            right_hand_side_[i] *= row_scale[static_cast<size_t>(i)];
            left_block_scaling_[static_cast<size_t>(i)] *= row_scale[static_cast<size_t>(i)];
            right_block_scaling_[static_cast<size_t>(i)] *= col_scale[static_cast<size_t>(i)];
        }
    }

    used_diagonal_equilibration_iterations_ = iterations;
}
//=================================================================================================//
VecC LaplaceStructuredAPhiEigenSolver::
    assembleRightHandSideFromExactSolution(const StdVec<Real> &electrical_conductivity,
                                           const StdVec<Real> &magnetic_reluctivity,
                                           const LaplaceStructuredAPhiFields &exact_fields)
{
    buildOperatorMatrices(electrical_conductivity, magnetic_reluctivity);
    buildSystemMatrix();
    // Build the physical right-hand side from the unscaled system first.
    // Any optional block scaling is applied later, consistently, inside
    // assembleSystem() to both the matrix and the physical RHS.
    return system_matrix_ * gatherUnknownVector(exact_fields);
}
//=================================================================================================//
Complex LaplaceStructuredAPhiEigenSolver::
    scaleBoundaryUnknown(int index, const Complex &physical_value) const
{
    return physical_value / right_block_scaling_[static_cast<size_t>(index)];
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::
    applyBoundaryCondition(const LaplaceStructuredAPhiBoundaryCondition &boundary_condition)
{
    if (!boundary_condition.is_dirichlet.empty())
    {
        for (size_t i = 0; i != number_of_particles_; ++i)
        {
            if (!boundary_condition.is_dirichlet[i])
            {
                continue;
            }
            overwrite_row_with_identity(system_matrix_, right_hand_side_, idxAx(i),
                                        scaleBoundaryUnknown(idxAx(i), boundary_condition.ax[i]));
            overwrite_row_with_identity(system_matrix_, right_hand_side_, idxAy(i),
                                        scaleBoundaryUnknown(idxAy(i), boundary_condition.ay[i]));
            overwrite_row_with_identity(system_matrix_, right_hand_side_, idxAz(i),
                                        scaleBoundaryUnknown(idxAz(i), boundary_condition.az[i]));
            overwrite_row_with_identity(system_matrix_, right_hand_side_, idxPhi(i),
                                        scaleBoundaryUnknown(idxPhi(i), boundary_condition.phi[i]));
        }
    }

    if (parameters_.fix_phi_reference && parameters_.phi_reference_index < number_of_particles_)
    {
        overwrite_row_with_identity(system_matrix_, right_hand_side_,
                                    idxPhi(parameters_.phi_reference_index),
                                    scaleBoundaryUnknown(idxPhi(parameters_.phi_reference_index),
                                                         Complex(parameters_.phi_reference_value, 0.0)));
    }
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::
    assembleSystem(const StdVec<Real> &electrical_conductivity,
                   const StdVec<Real> &magnetic_reluctivity,
                   const VecC &rhs_unconstrained,
                   const LaplaceStructuredAPhiBoundaryCondition &boundary_condition)
{
    buildOperatorMatrices(electrical_conductivity, magnetic_reluctivity);
    buildSystemMatrix();
    rhs_unconstrained_physical_ = rhs_unconstrained;
    right_hand_side_ = rhs_unconstrained;
    initializeBlockScaling();
    applyBlockScalingToSystem();
    applyDiagonalEquilibrationToSystem();
    applyBoundaryCondition(boundary_condition);
}
//=================================================================================================//
VecC LaplaceStructuredAPhiEigenSolver::
    gatherUnknownVector(const LaplaceStructuredAPhiFields &fields) const
{
    VecC unknown_vector(static_cast<int>(4 * number_of_particles_));
    unknown_vector.setZero();
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        unknown_vector[idxAx(i)] = fields.ax[i];
        unknown_vector[idxAy(i)] = fields.ay[i];
        unknown_vector[idxAz(i)] = fields.az[i];
        unknown_vector[idxPhi(i)] = fields.phi[i];
    }
    return unknown_vector;
}
//=================================================================================================//
VecC LaplaceStructuredAPhiEigenSolver::
    physicalToScaledUnknownVector(const VecC &physical_unknown_vector) const
{
    VecC scaled_unknown_vector = physical_unknown_vector;
    for (int i = 0; i != scaled_unknown_vector.size(); ++i)
    {
        scaled_unknown_vector[i] /= right_block_scaling_[static_cast<size_t>(i)];
    }
    return scaled_unknown_vector;
}
//=================================================================================================//
VecC LaplaceStructuredAPhiEigenSolver::
    scaledToPhysicalUnknownVector(const VecC &scaled_unknown_vector) const
{
    VecC physical_unknown_vector = scaled_unknown_vector;
    for (int i = 0; i != physical_unknown_vector.size(); ++i)
    {
        physical_unknown_vector[i] *= right_block_scaling_[static_cast<size_t>(i)];
    }
    return physical_unknown_vector;
}
//=================================================================================================//
VecC LaplaceStructuredAPhiEigenSolver::
    gatherScaledUnknownVector(const LaplaceStructuredAPhiFields &fields) const
{
    return physicalToScaledUnknownVector(gatherUnknownVector(fields));
}
//=================================================================================================//
void LaplaceStructuredAPhiEigenSolver::
    scatterUnknownVector(const VecC &unknown_vector,
                         LaplaceStructuredAPhiFields &fields) const
{
    fields.ax.resize(number_of_particles_);
    fields.ay.resize(number_of_particles_);
    fields.az.resize(number_of_particles_);
    fields.phi.resize(number_of_particles_);
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        fields.ax[i] = unknown_vector[idxAx(i)];
        fields.ay[i] = unknown_vector[idxAy(i)];
        fields.az[i] = unknown_vector[idxAz(i)];
        fields.phi[i] = unknown_vector[idxPhi(i)];
    }
}
//=================================================================================================//
bool LaplaceStructuredAPhiEigenSolver::
    solve(LaplaceStructuredAPhiFields &solution_fields,
          std::string &solver_message)
{
    std::string backend = parameters_.solver_backend;
    std::transform(backend.begin(), backend.end(), backend.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });

    last_solver_backend_ = backend;
    last_solver_iterations_ = 0;
    last_solver_estimated_error_ = 0.0;

    if (backend == "bicgstab")
    {
        Eigen::BiCGSTAB<SpMatC, Eigen::IncompleteLUT<Complex>> solver;
        solver.setMaxIterations(parameters_.iterative_max_iterations);
        solver.setTolerance(parameters_.iterative_tolerance);
        solver.compute(system_matrix_);
        if (solver.info() != Eigen::Success)
        {
            solver_message = "BiCGSTAB setup failed.";
            return false;
        }

        VecC solution_vector = solver.solve(right_hand_side_);
        last_solver_iterations_ = solver.iterations();
        last_solver_estimated_error_ = static_cast<Real>(solver.error());
        if (solver.info() != Eigen::Success)
        {
            solver_message = "BiCGSTAB solve failed.";
            return false;
        }

        scatterUnknownVector(scaledToPhysicalUnknownVector(solution_vector), solution_fields);
        solver_message = "BiCGSTAB solve succeeded.";
        return true;
    }

    Eigen::SparseLU<SpMatC> solver;
    solver.analyzePattern(system_matrix_);
    solver.factorize(system_matrix_);
    if (solver.info() != Eigen::Success)
    {
        solver_message = "SparseLU factorization failed.";
        return false;
    }

    VecC solution_vector = solver.solve(right_hand_side_);
    if (solver.info() != Eigen::Success)
    {
        solver_message = "SparseLU solve failed.";
        return false;
    }

    scatterUnknownVector(scaledToPhysicalUnknownVector(solution_vector), solution_fields);
    solver_message = "SparseLU solve succeeded.";
    return true;
}
//=================================================================================================//
bool LaplaceStructuredAPhiEigenSolver::
    solveDiscreteCompatiblePhi(const StdVec<Real> &electrical_conductivity,
                               const StdVec<Real> &magnetic_reluctivity,
                               const LaplaceStructuredAPhiBoundaryCondition &boundary_condition,
                               LaplaceStructuredAPhiFields &fields,
                               std::string &solver_message)
{
    buildOperatorMatrices(electrical_conductivity, magnetic_reluctivity);
    Real geom_length_to_m = SMAX(parameters_.geom_length_to_m, TinyReal);
    Real grad_scale = static_cast<Real>(1.0) / geom_length_to_m;
    Real laplace_scale = grad_scale * grad_scale;

    VecC ax(static_cast<int>(number_of_particles_));
    VecC ay(static_cast<int>(number_of_particles_));
    VecC az(static_cast<int>(number_of_particles_));
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        ax[static_cast<int>(i)] = fields.ax[i];
        ay[static_cast<int>(i)] = fields.ay[i];
        az[static_cast<int>(i)] = fields.az[i];
    }

    SpMatC mass_sigma = make_diagonal_sparse_complex(electrical_conductivity);
    VecC right_hand_side =
        -Complex(0.0, parameters_.angular_frequency * grad_scale) *
        (grad_x_.transpose() * (mass_sigma * ax) +
         grad_y_.transpose() * (mass_sigma * ay) +
         grad_z_.transpose() * (mass_sigma * az));

    SpMatC phi_matrix = Complex(laplace_scale, 0.0) * laplace_sigma_;
    if (!boundary_condition.is_dirichlet.empty())
    {
        for (size_t i = 0; i != number_of_particles_; ++i)
        {
            if (!boundary_condition.is_dirichlet[i])
            {
                continue;
            }
            overwrite_row_with_identity(phi_matrix, right_hand_side, static_cast<int>(i),
                                        boundary_condition.phi[i]);
        }
    }

    if (parameters_.fix_phi_reference && parameters_.phi_reference_index < number_of_particles_)
    {
        overwrite_row_with_identity(phi_matrix, right_hand_side,
                                    static_cast<int>(parameters_.phi_reference_index),
                                    Complex(parameters_.phi_reference_value, 0.0));
    }

    Eigen::SparseLU<SpMatC> solver;
    solver.analyzePattern(phi_matrix);
    solver.factorize(phi_matrix);
    if (solver.info() != Eigen::Success)
    {
        solver_message = "Compatible phi SparseLU factorization failed.";
        return false;
    }

    VecC phi_solution = solver.solve(right_hand_side);
    if (solver.info() != Eigen::Success)
    {
        solver_message = "Compatible phi SparseLU solve failed.";
        return false;
    }

    fields.phi.resize(number_of_particles_);
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        fields.phi[i] = phi_solution[static_cast<int>(i)];
    }
    solver_message = "Compatible phi SparseLU solve succeeded.";
    return true;
}
//=================================================================================================//
LaplaceStructuredAPhiDiagnostics LaplaceStructuredAPhiEigenSolver::
    computeDiagnostics(const StdVec<Real> &electrical_conductivity,
                       const StdVec<Real> &magnetic_reluctivity,
                       LaplaceStructuredAPhiFields &fields) const
{
    (void)magnetic_reluctivity;
    Real geom_length_to_m = SMAX(parameters_.geom_length_to_m, TinyReal);
    Real grad_scale = static_cast<Real>(1.0) / geom_length_to_m;
    Real laplace_scale = grad_scale * grad_scale;
    LaplaceStructuredAPhiDiagnostics diagnostics;
    VecC unknown_vector = gatherScaledUnknownVector(fields);
    VecC residual = system_matrix_ * unknown_vector - right_hand_side_;
    diagnostics.linear_residual_norm = residual.norm();

    VecC ax(static_cast<int>(number_of_particles_));
    VecC ay(static_cast<int>(number_of_particles_));
    VecC az(static_cast<int>(number_of_particles_));
    VecC phi(static_cast<int>(number_of_particles_));
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        ax[static_cast<int>(i)] = fields.ax[i];
        ay[static_cast<int>(i)] = fields.ay[i];
        az[static_cast<int>(i)] = fields.az[i];
        phi[static_cast<int>(i)] = fields.phi[i];
    }

    VecC grad_phi_x = Complex(grad_scale, 0.0) * (grad_x_ * phi);
    VecC grad_phi_y = Complex(grad_scale, 0.0) * (grad_y_ * phi);
    VecC grad_phi_z = Complex(grad_scale, 0.0) * (grad_z_ * phi);

    VecC div_a = Complex(grad_scale, 0.0) * (grad_x_ * ax + grad_y_ * ay + grad_z_ * az);
    diagnostics.divergence_a_l2 = div_a.norm();

    SpMatC mass_sigma = make_diagonal_sparse_complex(electrical_conductivity);
    VecC phi_diffusion_term = Complex(laplace_scale, 0.0) * (laplace_sigma_ * phi);
    VecC sigma_a_coupling_term =
        Complex(0.0, parameters_.angular_frequency * grad_scale) *
        (grad_x_.transpose() * (mass_sigma * ax) +
         grad_y_.transpose() * (mass_sigma * ay) +
         grad_z_.transpose() * (mass_sigma * az));
    VecC phi_lhs = phi_diffusion_term + sigma_a_coupling_term;
    diagnostics.phi_diffusion_term_l2 = phi_diffusion_term.norm();
    diagnostics.sigma_a_coupling_term_l2 = sigma_a_coupling_term.norm();

    VecC phi_row_residual(static_cast<int>(number_of_particles_));
    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        phi_row_residual[static_cast<int>(i)] = residual[idxPhi(i)];
    }
    diagnostics.phi_equation_residual_l2 = phi_row_residual.norm();

    VecC rhs_phi = VecC::Zero(static_cast<int>(number_of_particles_));
    if (rhs_unconstrained_physical_.size() == static_cast<int>(4 * number_of_particles_))
    {
        for (size_t i = 0; i != number_of_particles_; ++i)
        {
            rhs_phi[static_cast<int>(i)] = rhs_unconstrained_physical_[idxPhi(i)];
        }
    }
    diagnostics.rhs_phi_l2 = rhs_phi.norm();

    VecC physical_continuity_residual = phi_lhs - rhs_phi;
    diagnostics.physical_current_continuity_l2 = physical_continuity_residual.norm();
    diagnostics.current_continuity_l2 = diagnostics.physical_current_continuity_l2;
    Real continuity_scale = phi_lhs.norm() + rhs_phi.norm() + TinyReal;
    diagnostics.relative_current_continuity_l2 =
        diagnostics.physical_current_continuity_l2 / continuity_scale;

    fields.electric_field.resize(number_of_particles_);
    fields.magnetic_flux_density.resize(number_of_particles_);
    fields.current_density.resize(number_of_particles_);
    fields.joule_heat.resize(number_of_particles_);

    VecC curl_x = Complex(grad_scale, 0.0) * (grad_y_ * az - grad_z_ * ay);
    VecC curl_y = Complex(grad_scale, 0.0) * (grad_z_ * ax - grad_x_ * az);
    VecC curl_z = Complex(grad_scale, 0.0) * (grad_x_ * ay - grad_y_ * ax);

    for (size_t i = 0; i != number_of_particles_; ++i)
    {
        Vec3c electric_field;
        electric_field[0] = -Complex(0.0, parameters_.angular_frequency) * fields.ax[i] - grad_phi_x[static_cast<int>(i)];
        electric_field[1] = -Complex(0.0, parameters_.angular_frequency) * fields.ay[i] - grad_phi_y[static_cast<int>(i)];
        electric_field[2] = -Complex(0.0, parameters_.angular_frequency) * fields.az[i] - grad_phi_z[static_cast<int>(i)];
        fields.electric_field[i] = electric_field;

        Vec3c magnetic_flux_density;
        magnetic_flux_density[0] = curl_x[static_cast<int>(i)];
        magnetic_flux_density[1] = curl_y[static_cast<int>(i)];
        magnetic_flux_density[2] = curl_z[static_cast<int>(i)];
        fields.magnetic_flux_density[i] = magnetic_flux_density;

        Vec3c current_density = electrical_conductivity[i] * electric_field;
        fields.current_density[i] = current_density;

        Real electric_norm = std::sqrt(std::norm(electric_field[0]) +
                                       std::norm(electric_field[1]) +
                                       std::norm(electric_field[2]));
        Real current_norm = std::sqrt(std::norm(current_density[0]) +
                                      std::norm(current_density[1]) +
                                      std::norm(current_density[2]));
        Real a_norm = std::sqrt(std::norm(fields.ax[i]) +
                                std::norm(fields.ay[i]) +
                                std::norm(fields.az[i]));

        diagnostics.max_abs_a = SMAX(diagnostics.max_abs_a, a_norm);
        diagnostics.max_abs_e = SMAX(diagnostics.max_abs_e, electric_norm);
        diagnostics.max_abs_j = SMAX(diagnostics.max_abs_j, current_norm);

        Real joule_heat = 0.5 * electrical_conductivity[i] *
                          (std::norm(electric_field[0]) +
                           std::norm(electric_field[1]) +
                           std::norm(electric_field[2]));
        fields.joule_heat[i] = joule_heat;
        diagnostics.total_joule_power += volumetric_measure_[i] * joule_heat;
        diagnostics.gauge_penalty_energy += volumetric_measure_[i] * std::norm(div_a[static_cast<int>(i)]);
        diagnostics.magnetic_energy += volumetric_measure_[i] *
                                       (std::norm(curl_x[static_cast<int>(i)]) +
                                        std::norm(curl_y[static_cast<int>(i)]) +
                                        std::norm(curl_z[static_cast<int>(i)]));
    }

    diagnostics.gauge_penalty_energy = std::sqrt(diagnostics.gauge_penalty_energy);
    diagnostics.magnetic_energy *= 0.5;
    return diagnostics;
}
//=================================================================================================//

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_HPP
