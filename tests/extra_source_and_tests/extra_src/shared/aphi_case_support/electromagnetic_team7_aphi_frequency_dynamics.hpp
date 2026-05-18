#ifndef ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_HPP
#define ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_HPP

#include "aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.h"
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace
{
inline Real ComputeLocalMagneticJacobiDiagonal(const Neighborhood &inner_neighborhood,
                                               size_t index_i,
                                               Real *Vol,
                                               Real *magnetic_reluctivity,
                                               Real magnetic_diagonal_scaling)
{
    Real diagonal = 0.0;
    Real nu_i = magnetic_reluctivity[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j =
            inner_neighborhood.dW_ij_[n] * Vol[index_j] * inner_neighborhood.e_ij_[n];
        Real grad_norm_sq = gradW_ijV_j.squaredNorm();
        if (!std::isfinite(grad_norm_sq))
        {
            continue;
        }
        Real nu_ij = 0.5 * (nu_i + magnetic_reluctivity[index_j]);
        diagonal += nu_ij * grad_norm_sq;
    }
    return magnetic_diagonal_scaling * diagonal;
}

inline Real ComputeContactMagneticJacobiDiagonal(const StdVec<ParticleConfiguration *> &contact_configuration,
                                                 size_t index_i,
                                                 const StdVec<Real *> &contact_vol,
                                                 Real nu_i,
                                                 const StdVec<Real *> &contact_magnetic_reluctivity,
                                                 Real magnetic_diagonal_scaling)
{
    Real diagonal = 0.0;
    for (size_t k = 0; k != contact_configuration.size(); ++k)
    {
        Real *contact_vol_k = contact_vol[k];
        Real *contact_nu_k = contact_magnetic_reluctivity[k];
        Neighborhood &contact_neighborhood = (*contact_configuration[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j =
                contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Real grad_norm_sq = gradW_ijV_j.squaredNorm();
            if (!std::isfinite(grad_norm_sq))
            {
                continue;
            }
            Real nu_ij = 0.5 * (nu_i + contact_nu_k[index_j]);
            diagonal += nu_ij * grad_norm_sq;
        }
    }
    return magnetic_diagonal_scaling * diagonal;
}

inline Real ComputeLocalMagneticGradientRowSum(const Neighborhood &inner_neighborhood,
                                               size_t index_i,
                                               Real *Vol,
                                               Real *magnetic_reluctivity)
{
    Real row_sum = 0.0;
    Real nu_i = magnetic_reluctivity[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j =
            inner_neighborhood.dW_ij_[n] * Vol[index_j] * inner_neighborhood.e_ij_[n];
        Real grad_norm_sq = gradW_ijV_j.squaredNorm();
        if (!std::isfinite(grad_norm_sq))
        {
            continue;
        }
        Real nu_ij = 0.5 * (nu_i + magnetic_reluctivity[index_j]);
        row_sum += sqrt(SMAX(static_cast<Real>(0.0), nu_ij)) * sqrt(grad_norm_sq);
    }
    return row_sum;
}

inline Real ComputeContactMagneticGradientRowSum(const StdVec<ParticleConfiguration *> &contact_configuration,
                                                 size_t index_i,
                                                 const StdVec<Real *> &contact_vol,
                                                 Real nu_i,
                                                 const StdVec<Real *> &contact_magnetic_reluctivity)
{
    Real row_sum = 0.0;
    for (size_t k = 0; k != contact_configuration.size(); ++k)
    {
        Real *contact_vol_k = contact_vol[k];
        Real *contact_nu_k = contact_magnetic_reluctivity[k];
        Neighborhood &contact_neighborhood = (*contact_configuration[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j =
                contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Real grad_norm_sq = gradW_ijV_j.squaredNorm();
            if (!std::isfinite(grad_norm_sq))
            {
                continue;
            }
            Real nu_ij = 0.5 * (nu_i + contact_nu_k[index_j]);
            row_sum += sqrt(SMAX(static_cast<Real>(0.0), nu_ij)) * sqrt(grad_norm_sq);
        }
    }
    return row_sum;
}

inline Real ComputeConservativeMagneticDiagonal(const Neighborhood &inner_neighborhood,
                                                size_t index_i,
                                                Real *Vol,
                                                Real *magnetic_reluctivity,
                                                Real magnetic_diagonal_scaling)
{
    Real local_diagonal =
        ComputeLocalMagneticJacobiDiagonal(inner_neighborhood, index_i, Vol,
                                          magnetic_reluctivity, magnetic_diagonal_scaling);
    Real local_row_sum =
        ComputeLocalMagneticGradientRowSum(inner_neighborhood, index_i, Vol,
                                          magnetic_reluctivity);
    Real conservative_diagonal =
        magnetic_diagonal_scaling * local_row_sum * local_row_sum;
    return SMAX(local_diagonal, conservative_diagonal);
}

inline Real ComputeBalancedMagneticDiagonal(const Neighborhood &inner_neighborhood,
                                            size_t index_i,
                                            Real *Vol,
                                            Real *magnetic_reluctivity,
                                            Real magnetic_diagonal_scaling)
{
    Real jacobi_diagonal =
        ComputeLocalMagneticJacobiDiagonal(inner_neighborhood, index_i, Vol,
                                           magnetic_reluctivity, magnetic_diagonal_scaling);
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol,
                                            magnetic_reluctivity, magnetic_diagonal_scaling);
    Real safe_jacobi = SMAX(jacobi_diagonal, TinyReal);
    Real safe_conservative = SMAX(conservative_diagonal, safe_jacobi);
    Real balanced_diagonal = sqrt(safe_jacobi * safe_conservative);
    return SMAX(jacobi_diagonal, balanced_diagonal);
}

inline Real BlendMagneticDiagonal(Real conservative_diagonal,
                                  Real balanced_diagonal,
                                  Real balanced_weight)
{
    Real clamped_weight =
        SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), balanced_weight));
    return conservative_diagonal +
           clamped_weight * (balanced_diagonal - conservative_diagonal);
}

inline Real ComputeConservativeMagneticDiagonal(const Neighborhood &inner_neighborhood,
                                                size_t index_i,
                                                Real *Vol,
                                                Real *magnetic_reluctivity,
                                                const StdVec<ParticleConfiguration *> &contact_configuration,
                                                const StdVec<Real *> &contact_vol,
                                                const StdVec<Real *> &contact_magnetic_reluctivity,
                                                Real magnetic_diagonal_scaling)
{
    Real nu_i = magnetic_reluctivity[index_i];
    Real local_diagonal =
        ComputeLocalMagneticJacobiDiagonal(inner_neighborhood, index_i, Vol,
                                          magnetic_reluctivity, magnetic_diagonal_scaling);
    Real contact_diagonal =
        ComputeContactMagneticJacobiDiagonal(contact_configuration, index_i, contact_vol,
                                            nu_i, contact_magnetic_reluctivity,
                                            magnetic_diagonal_scaling);
    Real local_row_sum =
        ComputeLocalMagneticGradientRowSum(inner_neighborhood, index_i, Vol,
                                          magnetic_reluctivity);
    Real contact_row_sum =
        ComputeContactMagneticGradientRowSum(contact_configuration, index_i, contact_vol,
                                            nu_i, contact_magnetic_reluctivity);
    Real conservative_diagonal =
        magnetic_diagonal_scaling *
        (local_row_sum + contact_row_sum) * (local_row_sum + contact_row_sum);
    return SMAX(local_diagonal + contact_diagonal, conservative_diagonal);
}

inline Real ComputeJacobiMagneticDiagonal(const Neighborhood &inner_neighborhood,
                                          size_t index_i,
                                          Real *Vol,
                                          Real *magnetic_reluctivity,
                                          const StdVec<ParticleConfiguration *> &contact_configuration,
                                          const StdVec<Real *> &contact_vol,
                                          const StdVec<Real *> &contact_magnetic_reluctivity,
                                          Real magnetic_diagonal_scaling)
{
    Real nu_i = magnetic_reluctivity[index_i];
    Real local_diagonal =
        ComputeLocalMagneticJacobiDiagonal(inner_neighborhood, index_i, Vol,
                                           magnetic_reluctivity, magnetic_diagonal_scaling);
    Real contact_diagonal =
        ComputeContactMagneticJacobiDiagonal(contact_configuration, index_i, contact_vol,
                                             nu_i, contact_magnetic_reluctivity,
                                             magnetic_diagonal_scaling);
    return local_diagonal + contact_diagonal;
}

inline Real ComputeBalancedMagneticDiagonal(const Neighborhood &inner_neighborhood,
                                            size_t index_i,
                                            Real *Vol,
                                            Real *magnetic_reluctivity,
                                            const StdVec<ParticleConfiguration *> &contact_configuration,
                                            const StdVec<Real *> &contact_vol,
                                            const StdVec<Real *> &contact_magnetic_reluctivity,
                                            Real magnetic_diagonal_scaling)
{
    Real nu_i = magnetic_reluctivity[index_i];
    Real local_diagonal =
        ComputeLocalMagneticJacobiDiagonal(inner_neighborhood, index_i, Vol,
                                           magnetic_reluctivity, magnetic_diagonal_scaling);
    Real contact_diagonal =
        ComputeContactMagneticJacobiDiagonal(contact_configuration, index_i, contact_vol,
                                             nu_i, contact_magnetic_reluctivity,
                                             magnetic_diagonal_scaling);
    Real jacobi_diagonal = local_diagonal + contact_diagonal;
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol,
                                            magnetic_reluctivity,
                                            contact_configuration, contact_vol,
                                            contact_magnetic_reluctivity,
                                            magnetic_diagonal_scaling);
    Real safe_jacobi = SMAX(jacobi_diagonal, TinyReal);
    Real safe_conservative = SMAX(conservative_diagonal, safe_jacobi);
    Real balanced_diagonal = sqrt(safe_jacobi * safe_conservative);
    return SMAX(jacobi_diagonal, balanced_diagonal);
}

inline Real ApplyContactDiagonalRatioCap(Real diagonal_value,
                                         Real jacobi_diagonal,
                                         Real ratio_cap)
{
    if (!(ratio_cap > static_cast<Real>(1.0)))
    {
        return diagonal_value;
    }
    Real safe_jacobi = SMAX(jacobi_diagonal, TinyReal);
    Real capped_upper = ratio_cap * safe_jacobi;
    return SMAX(safe_jacobi, SMIN(diagonal_value, capped_upper));
}

inline Real ComputeAdaptiveContactDiagonalRatioCap(Real ratio_cap,
                                                   Real conservative_diagonal,
                                                   Real jacobi_diagonal,
                                                   Real adaptive_strength)
{
    if (!(ratio_cap > static_cast<Real>(1.0)) ||
        !(adaptive_strength > static_cast<Real>(0.0)))
    {
        return ratio_cap;
    }
    Real safe_jacobi = SMAX(jacobi_diagonal, TinyReal);
    Real safe_conservative = SMAX(conservative_diagonal, safe_jacobi);
    Real conservative_over_jacobi = safe_conservative / safe_jacobi;
    Real ratio_growth =
        sqrt(SMAX(static_cast<Real>(1.0), conservative_over_jacobi)) -
        static_cast<Real>(1.0);
    Real tighten_factor = static_cast<Real>(1.0) +
                          adaptive_strength * ratio_growth;
    Real adaptive_cap = ratio_cap / (tighten_factor + TinyReal);
    return SMAX(static_cast<Real>(2.0),
                SMIN(ratio_cap, adaptive_cap));
}

inline Vecd ClampVectorByNorm(const Vecd &value, Real max_norm)
{
    Real norm_sq = value.squaredNorm();
    if (!std::isfinite(norm_sq))
    {
        return ZeroData<Vecd>::value;
    }
    Real norm = sqrt(norm_sq);
    if (norm > max_norm)
    {
        return value * (max_norm / (norm + TinyReal));
    }
    return value;
}

inline Vecd ComputeAxisAlignedBoxBoundaryNormal(const Vecd &position,
                                                const Vecd &box_center,
                                                const Vecd &box_halfsize)
{
    Vecd relative = position - box_center;
    size_t normal_axis = 0;
    Real min_face_distance = MaxReal;
    for (size_t d = 0; d != Dimensions; ++d)
    {
        Real half_extent = box_halfsize[d];
        if (!(half_extent > TinyReal))
        {
            continue;
        }
        Real face_distance = fabs(half_extent - fabs(relative[d]));
        if (face_distance < min_face_distance)
        {
            min_face_distance = face_distance;
            normal_axis = d;
        }
    }

    Vecd normal = ZeroData<Vecd>::value;
    Real sign = relative[normal_axis] >= 0.0 ? 1.0 : -1.0;
    normal[normal_axis] = sign;
    return normal;
}

inline Real ComputeNormalizedPseudoTimeStep(Real dt, Real reference_pseudo_time_step)
{
    Real safe_reference_dt = reference_pseudo_time_step;
    if (!(safe_reference_dt > TinyReal))
    {
        safe_reference_dt = 1.0;
    }
    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        return 0.0;
    }
    return dt / safe_reference_dt;
}

inline Matd ComputeCurlCurlPairSelfBlock(const Vecd &gradW_ijV_j, Real magnetic_reluctivity_j)
{
    Real grad_norm_sq = gradW_ijV_j.squaredNorm();
    if (!(grad_norm_sq > TinyReal) || !std::isfinite(grad_norm_sq))
    {
        return ZeroData<Matd>::value;
    }
    return magnetic_reluctivity_j *
           (grad_norm_sq * Matd::Identity() - gradW_ijV_j * gradW_ijV_j.transpose());
}

inline Matd ComputeLocalMagneticSelfBlock(const Neighborhood &inner_neighborhood,
                                          size_t index_i,
                                          Real *Vol,
                                          Real *magnetic_reluctivity,
                                          Real magnetic_diagonal_scaling)
{
    Matd self_block = ZeroData<Matd>::value;
    Vecd accumulated_gradient = ZeroData<Vecd>::value;
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j =
            inner_neighborhood.dW_ij_[n] * Vol[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        accumulated_gradient += gradW_ijV_j;
        self_block += ComputeCurlCurlPairSelfBlock(gradW_ijV_j, magnetic_reluctivity[index_j]);
    }
    self_block -= ComputeCurlCurlPairSelfBlock(accumulated_gradient, magnetic_reluctivity[index_i]);
    return magnetic_diagonal_scaling * self_block;
}

inline Matd ComputeLocalMagneticSelfBlock(const Neighborhood &inner_neighborhood,
                                          size_t index_i,
                                          Real *Vol,
                                          Real *magnetic_reluctivity,
                                          const StdVec<ParticleConfiguration *> &contact_configuration,
                                          const StdVec<Real *> &contact_vol,
                                          const StdVec<Real *> &contact_magnetic_reluctivity,
                                          Real magnetic_diagonal_scaling)
{
    Matd self_block = ZeroData<Matd>::value;
    Vecd accumulated_gradient = ZeroData<Vecd>::value;
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j =
            inner_neighborhood.dW_ij_[n] * Vol[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        accumulated_gradient += gradW_ijV_j;
        self_block += ComputeCurlCurlPairSelfBlock(gradW_ijV_j, magnetic_reluctivity[index_j]);
    }

    for (size_t k = 0; k != contact_configuration.size(); ++k)
    {
        Real *contact_vol_k = contact_vol[k];
        Real *contact_nu_k = contact_magnetic_reluctivity[k];
        Neighborhood &contact_neighborhood = (*contact_configuration[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j =
                contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j] *
                contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            accumulated_gradient += gradW_ijV_j;
            self_block += ComputeCurlCurlPairSelfBlock(gradW_ijV_j, contact_nu_k[index_j]);
        }
    }

    self_block -= ComputeCurlCurlPairSelfBlock(accumulated_gradient, magnetic_reluctivity[index_i]);
    return magnetic_diagonal_scaling * self_block;
}

inline Vecd SolveMagneticSelfBlockOrFallback(const Matd &self_block,
                                             const Vecd &residual,
                                             Real fallback_diagonal)
{
    Vecd fallback_delta = residual / (fallback_diagonal + TinyReal);
    Matd regularized_block =
        self_block + (TinyReal * (fabs(fallback_diagonal) + 1.0)) * Matd::Identity();
    Eigen::LDLT<Matd> ldlt(regularized_block);
    if (ldlt.info() != Eigen::Success)
    {
        return fallback_delta;
    }

    Vecd block_delta = ldlt.solve(residual);
    if (!std::isfinite(block_delta.squaredNorm()))
    {
        return fallback_delta;
    }

    Real block_norm_sq = block_delta.squaredNorm();
    Real fallback_norm_sq = fallback_delta.squaredNorm();
    if (!(block_norm_sq > TinyReal))
    {
        return block_delta;
    }
    if (!(fallback_norm_sq > TinyReal) || !std::isfinite(fallback_norm_sq))
    {
        return block_delta;
    }
    if (block_norm_sq <= fallback_norm_sq + TinyReal)
    {
        return block_delta;
    }
    return block_delta * sqrt(fallback_norm_sq / (block_norm_sq + TinyReal));
}

using CoupledBlockMat = Eigen::Matrix<Real, 2 * Dimensions, 2 * Dimensions>;
using CoupledBlockVec = Eigen::Matrix<Real, 2 * Dimensions, 1>;

inline void SolveCoupledMagneticSelfBlockOrFallback(const Matd &self_block,
                                                    const Vecd &residual_real,
                                                    const Vecd &residual_imag,
                                                    Real omega_sigma,
                                                    Real fallback_operator_diagonal,
                                                    Vecd &delta_real,
                                                    Vecd &delta_imag)
{
    Real fallback_denom =
        fallback_operator_diagonal * fallback_operator_diagonal +
        omega_sigma * omega_sigma + TinyReal;
    Vecd fallback_delta_real =
        (fallback_operator_diagonal * residual_real + omega_sigma * residual_imag) / fallback_denom;
    Vecd fallback_delta_imag =
        (fallback_operator_diagonal * residual_imag - omega_sigma * residual_real) / fallback_denom;

    CoupledBlockMat coupled_block = CoupledBlockMat::Zero();
    coupled_block.template block<Dimensions, Dimensions>(0, 0) = self_block;
    coupled_block.template block<Dimensions, Dimensions>(Dimensions, Dimensions) = self_block;
    coupled_block.template block<Dimensions, Dimensions>(0, Dimensions) =
        -omega_sigma * Matd::Identity();
    coupled_block.template block<Dimensions, Dimensions>(Dimensions, 0) =
        omega_sigma * Matd::Identity();

    Real regularization_scale =
        TinyReal * (fabs(fallback_operator_diagonal) + fabs(omega_sigma) + 1.0);
    coupled_block += regularization_scale * CoupledBlockMat::Identity();

    CoupledBlockVec rhs = CoupledBlockVec::Zero();
    rhs.template head<Dimensions>() = residual_real;
    rhs.template tail<Dimensions>() = residual_imag;

    Eigen::LDLT<CoupledBlockMat> ldlt(coupled_block);
    if (ldlt.info() != Eigen::Success)
    {
        delta_real = fallback_delta_real;
        delta_imag = fallback_delta_imag;
        return;
    }

    CoupledBlockVec block_delta = ldlt.solve(rhs);
    if (!std::isfinite(block_delta.squaredNorm()))
    {
        delta_real = fallback_delta_real;
        delta_imag = fallback_delta_imag;
        return;
    }

    Vecd block_delta_real = block_delta.template head<Dimensions>();
    Vecd block_delta_imag = block_delta.template tail<Dimensions>();
    Real block_norm_sq =
        block_delta_real.squaredNorm() + block_delta_imag.squaredNorm();
    Real fallback_norm_sq =
        fallback_delta_real.squaredNorm() + fallback_delta_imag.squaredNorm();
    if (!(block_norm_sq > TinyReal))
    {
        delta_real = block_delta_real;
        delta_imag = block_delta_imag;
        return;
    }
    if (!(fallback_norm_sq > TinyReal) || !std::isfinite(fallback_norm_sq))
    {
        delta_real = block_delta_real;
        delta_imag = block_delta_imag;
        return;
    }
    if (block_norm_sq <= fallback_norm_sq + TinyReal)
    {
        delta_real = block_delta_real;
        delta_imag = block_delta_imag;
        return;
    }

    Real scale = sqrt(fallback_norm_sq / (block_norm_sq + TinyReal));
    delta_real = block_delta_real * scale;
    delta_imag = block_delta_imag * scale;
}
} // namespace
//=================================================================================================//
InitializeAphiFrequencyElectromagneticVariables::
    InitializeAphiFrequencyElectromagneticVariables(SPHBody &sph_body,
                                                    Real default_conductivity,
                                                    Real default_rho_cp,
                                                    Real default_magnetic_reluctivity)
    : LocalDynamics(sph_body),
      default_conductivity_(default_conductivity),
      default_rho_cp_(default_rho_cp),
      default_magnetic_reluctivity_(default_magnetic_reluctivity),
      vector_potential_real_(particles_->registerStateVariableData<Vecd>("VectorPotentialReal")),
      vector_potential_imag_(particles_->registerStateVariableData<Vecd>("VectorPotentialImag")),
      vector_potential_change_rate_real_(particles_->registerStateVariableData<Vecd>("VectorPotentialChangeRateReal")),
      vector_potential_change_rate_imag_(particles_->registerStateVariableData<Vecd>("VectorPotentialChangeRateImag")),
      electric_potential_gradient_real_(particles_->registerStateVariableData<Vecd>("ElectricPotentialGradientReal")),
      electric_potential_gradient_imag_(particles_->registerStateVariableData<Vecd>("ElectricPotentialGradientImag")),
      source_current_density_real_(particles_->registerStateVariableData<Vecd>("SourceCurrentDensityReal")),
      source_current_density_imag_(particles_->registerStateVariableData<Vecd>("SourceCurrentDensityImag")),
      curl_nu_b_real_(particles_->registerStateVariableData<Vecd>("CurlNuBReal")),
      curl_nu_b_imag_(particles_->registerStateVariableData<Vecd>("CurlNuBImag")),
      vector_potential_curl_real_(particles_->registerStateVariableData<AngularVecd>("VectorPotentialCurlReal")),
      vector_potential_curl_imag_(particles_->registerStateVariableData<AngularVecd>("VectorPotentialCurlImag")),
      electric_potential_real_(particles_->registerStateVariableData<Real>("ElectricPotentialReal")),
      electric_potential_imag_(particles_->registerStateVariableData<Real>("ElectricPotentialImag")),
      electric_potential_source_real_(particles_->registerStateVariableData<Real>("ElectricPotentialSourceReal")),
      electric_potential_source_imag_(particles_->registerStateVariableData<Real>("ElectricPotentialSourceImag")),
      electric_potential_change_rate_real_(particles_->registerStateVariableData<Real>("ElectricPotentialChangeRateReal")),
      electric_potential_change_rate_imag_(particles_->registerStateVariableData<Real>("ElectricPotentialChangeRateImag")),
      electric_field_real_(particles_->registerStateVariableData<Vecd>("ElectricFieldReal")),
      electric_field_imag_(particles_->registerStateVariableData<Vecd>("ElectricFieldImag")),
      current_density_real_(particles_->registerStateVariableData<Vecd>("CurrentDensityReal")),
      current_density_imag_(particles_->registerStateVariableData<Vecd>("CurrentDensityImag")),
      joule_heat_source_(particles_->registerStateVariableData<Real>("JouleHeatSource")),
      temperature_change_rate_by_joule_(particles_->registerStateVariableData<Real>("TemperatureChangeRateByJoule")),
      electrical_conductivity_(particles_->registerStateVariableData<Real>("ElectricalConductivity", default_conductivity_)),
      rho_cp_(particles_->registerStateVariableData<Real>("RhoCp", default_rho_cp_)),
      magnetic_reluctivity_(particles_->registerStateVariableData<Real>("MagneticReluctivity", default_magnetic_reluctivity_))
{
    particles_->addVariableToWrite<Vecd>("VectorPotentialReal");
    particles_->addVariableToWrite<Vecd>("VectorPotentialImag");
    particles_->addVariableToWrite<Vecd>("VectorPotentialChangeRateReal");
    particles_->addVariableToWrite<Vecd>("VectorPotentialChangeRateImag");
    particles_->addVariableToWrite<Vecd>("SourceCurrentDensityReal");
    particles_->addVariableToWrite<Vecd>("SourceCurrentDensityImag");
    particles_->addVariableToWrite<Real>("ElectricPotentialReal");
    particles_->addVariableToWrite<Real>("ElectricPotentialImag");
    particles_->addVariableToWrite<Vecd>("ElectricPotentialGradientReal");
    particles_->addVariableToWrite<Vecd>("ElectricPotentialGradientImag");
    particles_->addVariableToWrite<Vecd>("CurlNuBReal");
    particles_->addVariableToWrite<Vecd>("CurlNuBImag");
    particles_->addVariableToWrite<AngularVecd>("VectorPotentialCurlReal");
    particles_->addVariableToWrite<AngularVecd>("VectorPotentialCurlImag");
    particles_->addVariableToWrite<Vecd>("ElectricFieldReal");
    particles_->addVariableToWrite<Vecd>("ElectricFieldImag");
    particles_->addVariableToWrite<Vecd>("CurrentDensityReal");
    particles_->addVariableToWrite<Vecd>("CurrentDensityImag");
    particles_->addVariableToWrite<Real>("JouleHeatSource");
    particles_->addVariableToWrite<Real>("TemperatureChangeRateByJoule");
    particles_->addVariableToWrite<Real>("ElectricalConductivity");
    particles_->addVariableToWrite<Real>("RhoCp");
    particles_->addVariableToWrite<Real>("MagneticReluctivity");
}
//=================================================================================================//
void InitializeAphiFrequencyElectromagneticVariables::update(size_t index_i, Real dt)
{
    (void)dt;
    if (electrical_conductivity_[index_i] <= TinyReal)
    {
        electrical_conductivity_[index_i] = default_conductivity_;
    }
    if (rho_cp_[index_i] <= TinyReal)
    {
        rho_cp_[index_i] = default_rho_cp_;
    }
    if (magnetic_reluctivity_[index_i] <= TinyReal)
    {
        magnetic_reluctivity_[index_i] = default_magnetic_reluctivity_;
    }

    vector_potential_real_[index_i] = ZeroData<Vecd>::value;
    vector_potential_imag_[index_i] = ZeroData<Vecd>::value;
    vector_potential_change_rate_real_[index_i] = ZeroData<Vecd>::value;
    vector_potential_change_rate_imag_[index_i] = ZeroData<Vecd>::value;
    electric_potential_gradient_real_[index_i] = ZeroData<Vecd>::value;
    electric_potential_gradient_imag_[index_i] = ZeroData<Vecd>::value;
    source_current_density_real_[index_i] = ZeroData<Vecd>::value;
    source_current_density_imag_[index_i] = ZeroData<Vecd>::value;
    curl_nu_b_real_[index_i] = ZeroData<Vecd>::value;
    curl_nu_b_imag_[index_i] = ZeroData<Vecd>::value;
    vector_potential_curl_real_[index_i] = ZeroData<AngularVecd>::value;
    vector_potential_curl_imag_[index_i] = ZeroData<AngularVecd>::value;
    electric_potential_real_[index_i] = 0.0;
    electric_potential_imag_[index_i] = 0.0;
    electric_potential_source_real_[index_i] = 0.0;
    electric_potential_source_imag_[index_i] = 0.0;
    electric_potential_change_rate_real_[index_i] = 0.0;
    electric_potential_change_rate_imag_[index_i] = 0.0;
    electric_field_real_[index_i] = ZeroData<Vecd>::value;
    electric_field_imag_[index_i] = ZeroData<Vecd>::value;
    current_density_real_[index_i] = ZeroData<Vecd>::value;
    current_density_imag_[index_i] = ZeroData<Vecd>::value;
    joule_heat_source_[index_i] = 0.0;
    temperature_change_rate_by_joule_[index_i] = 0.0;
}
//=================================================================================================//
PrescribedComplexSourceCurrentDensity::
    PrescribedComplexSourceCurrentDensity(SPHBody &sph_body,
                                          const Vecd &source_current_density_real,
                                          const Vecd &source_current_density_imag)
    : LocalDynamics(sph_body),
      source_current_density_real_value_(source_current_density_real),
      source_current_density_imag_value_(source_current_density_imag),
      source_current_density_real_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityReal")),
      source_current_density_imag_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityImag")) {}
//=================================================================================================//
void PrescribedComplexSourceCurrentDensity::update(size_t index_i, Real dt)
{
    (void)dt;
    source_current_density_real_[index_i] = source_current_density_real_value_;
    source_current_density_imag_[index_i] = source_current_density_imag_value_;
}
//=================================================================================================//
ConstrainScalarFieldByName::
    ConstrainScalarFieldByName(BodyPartByParticle &body_part,
                               const std::string &field_name,
                               Real reference_value)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      reference_value_(reference_value),
      scalar_field_(particles_->getVariableDataByName<Real>(field_name)) {}
//=================================================================================================//
void ConstrainScalarFieldByName::update(size_t index_i, Real dt)
{
    (void)dt;
    scalar_field_[index_i] = reference_value_;
}
//=================================================================================================//
ConstrainVectorFieldByName::
    ConstrainVectorFieldByName(BodyPartByParticle &body_part,
                               const std::string &field_name,
                               const Vecd &reference_value)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      reference_value_(reference_value),
      vector_field_(particles_->getVariableDataByName<Vecd>(field_name)) {}
//=================================================================================================//
void ConstrainVectorFieldByName::update(size_t index_i, Real dt)
{
    (void)dt;
    vector_field_[index_i] = reference_value_;
}
//=================================================================================================//
ConstrainVectorFieldToAxisAlignedBoxNormalByName::
    ConstrainVectorFieldToAxisAlignedBoxNormalByName(BodyPartByParticle &body_part,
                                                     const std::string &field_name,
                                                     const Vecd &box_center,
                                                     const Vecd &box_halfsize)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      box_center_(box_center),
      box_halfsize_(box_halfsize),
      positions_(particles_->getVariableDataByName<Vecd>("Position")),
      vector_field_(particles_->getVariableDataByName<Vecd>(field_name)) {}
//=================================================================================================//
void ConstrainVectorFieldToAxisAlignedBoxNormalByName::update(size_t index_i, Real dt)
{
    (void)dt;
    Vecd normal = ComputeAxisAlignedBoxBoundaryNormal(positions_[index_i], box_center_, box_halfsize_);
    Vecd field_value = vector_field_[index_i];
    vector_field_[index_i] = field_value.dot(normal) * normal;
}
//=================================================================================================//
ElectricPotentialSourceFromVectorFieldInner::
    ElectricPotentialSourceFromVectorFieldInner(BaseInnerRelation &inner_relation,
                                                const std::string &vector_field_name,
                                                const std::string &source_name,
                                                Real vector_field_scaling,
                                                Real divergence_scaling)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      vector_field_scaling_(vector_field_scaling),
      divergence_scaling_(divergence_scaling),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_source_(particles_->getVariableDataByName<Real>(source_name)),
      vector_field_(particles_->getVariableDataByName<Vecd>(vector_field_name)) {}
//=================================================================================================//
void ElectricPotentialSourceFromVectorFieldInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real source_i = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * electrical_conductivity_[index_j] /
                        (electrical_conductivity_[index_i] + electrical_conductivity_[index_j] + TinyReal);
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        Vecd vector_diff = vector_field_scaling_ * (vector_field_[index_i] - vector_field_[index_j]);
        source_i -= divergence_scaling_ * sigma_ij * vector_diff.dot(gradW_ijV_j);
    }
    electric_potential_source_[index_i] = source_i;
}
//=================================================================================================//
ElectricPotentialSourceFromVectorFieldContact::
    ElectricPotentialSourceFromVectorFieldContact(BaseContactRelation &contact_relation,
                                                  const std::string &vector_field_name,
                                                  const std::string &source_name,
                                                  Real vector_field_scaling,
                                                  Real divergence_scaling)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      vector_field_scaling_(vector_field_scaling),
      divergence_scaling_(divergence_scaling),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_source_(particles_->getVariableDataByName<Real>(source_name)),
      vector_field_(particles_->getVariableDataByName<Vecd>(vector_field_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electrical_conductivity_.push_back(contact_particles->getVariableDataByName<Real>("ElectricalConductivity"));
        contact_vector_field_.push_back(contact_particles->getVariableDataByName<Vecd>(vector_field_name));
    }
}
//=================================================================================================//
void ElectricPotentialSourceFromVectorFieldContact::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real source_contact = 0.0;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_sigma_k = contact_electrical_conductivity_[k];
        Vecd *contact_vector_k = contact_vector_field_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * contact_sigma_k[index_j] /
                            (electrical_conductivity_[index_i] + contact_sigma_k[index_j] + TinyReal);
            Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] *
                               contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            Vecd vector_diff = vector_field_scaling_ * (vector_field_[index_i] - contact_vector_k[index_j]);
            source_contact -= divergence_scaling_ * sigma_ij * vector_diff.dot(gradW_ijV_j);
        }
    }
    electric_potential_source_[index_i] += source_contact;
}
//=================================================================================================//
ScalarRelaxationInnerByName::
    ScalarRelaxationInnerByName(BaseInnerRelation &inner_relation,
                                const std::string &potential_name,
                                const std::string &source_name,
                                const std::string &change_rate_name,
                                Real laplacian_scaling)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      laplacian_scaling_(laplacian_scaling),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      scalar_potential_(particles_->getVariableDataByName<Real>(potential_name)),
      scalar_source_(particles_->getVariableDataByName<Real>(source_name)),
      scalar_change_rate_(particles_->getVariableDataByName<Real>(change_rate_name)) {}
//=================================================================================================//
void ScalarRelaxationInnerByName::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real change_rate = scalar_source_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * electrical_conductivity_[index_j] /
                        (electrical_conductivity_[index_i] + electrical_conductivity_[index_j] + TinyReal);
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
        Real scalar_diff = scalar_potential_[index_j] - scalar_potential_[index_i];
        change_rate += laplacian_scaling_ * sigma_ij * scalar_diff * surface_area_ij;
    }
    scalar_change_rate_[index_i] = change_rate;
}
//=================================================================================================//
ScalarRelaxationContactByName::
    ScalarRelaxationContactByName(BaseContactRelation &contact_relation,
                                  const std::string &potential_name,
                                  const std::string &change_rate_name,
                                  Real laplacian_scaling)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      smoothing_length_(contact_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      laplacian_scaling_(laplacian_scaling),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      scalar_change_rate_(particles_->getVariableDataByName<Real>(change_rate_name)),
      scalar_potential_(particles_->getVariableDataByName<Real>(potential_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electrical_conductivity_.push_back(contact_particles->getVariableDataByName<Real>("ElectricalConductivity"));
        contact_scalar_potential_.push_back(contact_particles->getVariableDataByName<Real>(potential_name));
    }
}
//=================================================================================================//
void ScalarRelaxationContactByName::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real contact_change_rate = 0.0;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_sigma_k = contact_electrical_conductivity_[k];
        Real *contact_scalar_k = contact_scalar_potential_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * contact_sigma_k[index_j] /
                            (electrical_conductivity_[index_i] + contact_sigma_k[index_j] + TinyReal);
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
            Real scalar_diff = contact_scalar_k[index_j] - scalar_potential_[index_i];
            contact_change_rate += laplacian_scaling_ * sigma_ij * scalar_diff * surface_area_ij;
        }
    }
    scalar_change_rate_[index_i] += contact_change_rate;
}
//=================================================================================================//
ScalarRelaxationComplexByName::
    ScalarRelaxationComplexByName(BaseInnerRelation &inner_relation,
                                  BaseContactRelation &contact_relation,
                                  const std::string &potential_name,
                                  const std::string &source_name,
                                  const std::string &change_rate_name,
                                  Real laplacian_scaling,
                                  bool use_contact)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      laplacian_scaling_(laplacian_scaling),
      use_contact_(use_contact),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      scalar_potential_(particles_->getVariableDataByName<Real>(potential_name)),
      scalar_source_(particles_->getVariableDataByName<Real>(source_name)),
      scalar_change_rate_(particles_->getVariableDataByName<Real>(change_rate_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electrical_conductivity_.push_back(contact_particles->getVariableDataByName<Real>("ElectricalConductivity"));
        contact_scalar_potential_.push_back(contact_particles->getVariableDataByName<Real>(potential_name));
    }
}
//=================================================================================================//
void ScalarRelaxationComplexByName::interaction(size_t index_i, Real dt)
{
    Real residual = scalar_source_[index_i];
    Real diagonal = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * electrical_conductivity_[index_j] /
                        (electrical_conductivity_[index_i] + electrical_conductivity_[index_j] + TinyReal);
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
        Real operator_coeff = laplacian_scaling_ * sigma_ij * surface_area_ij;
        residual += operator_coeff * (scalar_potential_[index_j] - scalar_potential_[index_i]);
        diagonal += operator_coeff;
    }

    if (use_contact_)
    {
        for (size_t k = 0; k != contact_configuration_.size(); ++k)
        {
            Real *contact_vol_k = contact_vol_[k];
            Real *contact_sigma_k = contact_electrical_conductivity_[k];
            Real *contact_scalar_k = contact_scalar_potential_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * contact_sigma_k[index_j] /
                                (electrical_conductivity_[index_i] + contact_sigma_k[index_j] + TinyReal);
                Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j];
                Real r_ij = contact_neighborhood.r_ij_[n];
                Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
                Real operator_coeff = laplacian_scaling_ * sigma_ij * surface_area_ij;
                residual += operator_coeff * (contact_scalar_k[index_j] - scalar_potential_[index_i]);
                diagonal += operator_coeff;
            }
        }
    }

    Real safe_dt = (dt > TinyReal && std::isfinite(dt)) ? dt : 1.0;
    scalar_change_rate_[index_i] = residual / (diagonal * safe_dt + TinyReal);
}
//=================================================================================================//
UpdateScalarByRelaxationRateByName::
    UpdateScalarByRelaxationRateByName(SPHBody &sph_body,
                                       const std::string &scalar_name,
                                       const std::string &change_rate_name,
                                       Real relaxation_scaling,
                                       Real max_abs_value)
    : LocalDynamics(sph_body),
      relaxation_scaling_(relaxation_scaling),
      max_abs_value_(max_abs_value),
      scalar_field_(particles_->getVariableDataByName<Real>(scalar_name)),
      scalar_change_rate_(particles_->getVariableDataByName<Real>(change_rate_name)) {}
//=================================================================================================//
void UpdateScalarByRelaxationRateByName::update(size_t index_i, Real dt)
{
    Real delta = relaxation_scaling_ * dt * scalar_change_rate_[index_i];
    if (!std::isfinite(delta))
    {
        delta = 0.0;
    }
    scalar_field_[index_i] += delta;
    if (!std::isfinite(scalar_field_[index_i]))
    {
        scalar_field_[index_i] = 0.0;
    }
    scalar_field_[index_i] = SMIN(max_abs_value_, SMAX(-max_abs_value_, scalar_field_[index_i]));
}
//=================================================================================================//
ScalarGradientInnerByName::
    ScalarGradientInnerByName(BaseInnerRelation &inner_relation,
                              const std::string &scalar_name,
                              const std::string &gradient_name,
                              Real gradient_scaling)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      gradient_scaling_(gradient_scaling),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      scalar_field_(particles_->getVariableDataByName<Real>(scalar_name)),
      scalar_gradient_(particles_->getVariableDataByName<Vecd>(gradient_name)) {}
//=================================================================================================//
void ScalarGradientInnerByName::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Vecd gradient = ZeroData<Vecd>::value;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        gradient -= gradient_scaling_ * (scalar_field_[index_i] - scalar_field_[index_j]) * gradW_ijV_j;
    }
    scalar_gradient_[index_i] = gradient;
}
//=================================================================================================//
ScalarGradientContactByName::
    ScalarGradientContactByName(BaseContactRelation &contact_relation,
                                const std::string &scalar_name,
                                const std::string &gradient_name,
                                Real gradient_scaling)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      gradient_scaling_(gradient_scaling),
      scalar_field_(particles_->getVariableDataByName<Real>(scalar_name)),
      scalar_gradient_(particles_->getVariableDataByName<Vecd>(gradient_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_scalar_field_.push_back(contact_particles->getVariableDataByName<Real>(scalar_name));
    }
}
//=================================================================================================//
void ScalarGradientContactByName::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Vecd gradient_contact = ZeroData<Vecd>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_scalar_k = contact_scalar_field_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] *
                               contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            gradient_contact -= gradient_scaling_ * (scalar_field_[index_i] - contact_scalar_k[index_j]) * gradW_ijV_j;
        }
    }
    scalar_gradient_[index_i] += gradient_contact;
}
//=================================================================================================//
VectorPotentialFrequencyEquationInner::
    VectorPotentialFrequencyEquationInner(BaseInnerRelation &inner_relation,
                                          Real angular_frequency,
                                          Real omega_coupling_sign,
                                          const std::string &vector_potential_name,
                                          const std::string &coupled_vector_potential_name,
                                          const std::string &source_current_density_name,
                                          const std::string &electric_potential_gradient_name,
                                          const std::string &curl_nu_b_name,
                                          const std::string &change_rate_name,
                                          Real sigma_relaxation_scaling,
                                          Real sigma_relaxation_floor,
                                          Real magnetic_diagonal_scaling,
                                          Real reference_pseudo_time_step,
                                          Real relaxation_scaling,
                                          Real max_change_rate)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      angular_frequency_(angular_frequency),
      omega_coupling_sign_(omega_coupling_sign),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      coupled_vector_potential_(particles_->getVariableDataByName<Vecd>(coupled_vector_potential_name)),
      source_current_density_(particles_->getVariableDataByName<Vecd>(source_current_density_name)),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>(change_rate_name)) {}
//=================================================================================================//
void VectorPotentialFrequencyEquationInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real magnetic_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_, magnetic_diagonal_scaling_);
    Vecd rhs = source_current_density_[index_i] - curl_nu_b_[index_i] -
               sigma_i * electric_potential_gradient_[index_i] +
               omega_coupling_sign_ * sigma_i * angular_frequency_ * coupled_vector_potential_[index_i];
    Real operator_diagonal =
        magnetic_diagonal + sigma_relaxation_scaling_ * sigma_i * angular_frequency_ +
        sigma_relaxation_floor_;
    vector_potential_change_rate_[index_i] =
        ClampVectorByNorm(rhs / (operator_diagonal + TinyReal), max_change_rate_);
}
//=================================================================================================//
void VectorPotentialFrequencyEquationInner::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_[index_i] += pseudo_step_scale *
                                  vector_potential_change_rate_[index_i];
    if (!std::isfinite(vector_potential_[index_i].squaredNorm()))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyCoupledEquationInner::
    VectorPotentialFrequencyCoupledEquationInner(BaseInnerRelation &inner_relation,
                                                 Real angular_frequency,
                                                 const std::string &vector_potential_real_name,
                                                 const std::string &vector_potential_imag_name,
                                                 const std::string &source_current_density_real_name,
                                                 const std::string &source_current_density_imag_name,
                                                 const std::string &electric_potential_gradient_real_name,
                                                 const std::string &electric_potential_gradient_imag_name,
                                                 const std::string &curl_nu_b_real_name,
                                                 const std::string &curl_nu_b_imag_name,
                                                 const std::string &change_rate_real_name,
                                                 const std::string &change_rate_imag_name,
                                                 Real sigma_relaxation_scaling,
                                                 Real sigma_relaxation_floor,
                                                 Real magnetic_diagonal_scaling,
                                                 Real reference_pseudo_time_step,
                                                 Real relaxation_scaling,
                                                 Real max_change_rate)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      angular_frequency_(angular_frequency),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_real_(particles_->getVariableDataByName<Vecd>(vector_potential_real_name)),
      vector_potential_imag_(particles_->getVariableDataByName<Vecd>(vector_potential_imag_name)),
      source_current_density_real_(particles_->getVariableDataByName<Vecd>(source_current_density_real_name)),
      source_current_density_imag_(particles_->getVariableDataByName<Vecd>(source_current_density_imag_name)),
      electric_potential_gradient_real_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_real_name)),
      electric_potential_gradient_imag_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_imag_name)),
      curl_nu_b_real_(particles_->getVariableDataByName<Vecd>(curl_nu_b_real_name)),
      curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>(curl_nu_b_imag_name)),
      vector_potential_change_rate_real_(particles_->getVariableDataByName<Vecd>(change_rate_real_name)),
      vector_potential_change_rate_imag_(particles_->getVariableDataByName<Vecd>(change_rate_imag_name)) {}
//=================================================================================================//
void VectorPotentialFrequencyCoupledEquationInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real magnetic_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_, magnetic_diagonal_scaling_);
    Real operator_diagonal =
        magnetic_diagonal + sigma_relaxation_floor_;
    Real omega_sigma = sigma_relaxation_scaling_ * sigma_i * angular_frequency_;

    // Coupled update uses delta-corrections, so the residual must include the
    // current real/imag omega-coupling terms from the full phasor equation.
    Vecd residual_real = source_current_density_real_[index_i] - curl_nu_b_real_[index_i] -
                         sigma_i * electric_potential_gradient_real_[index_i] +
                         sigma_i * angular_frequency_ * vector_potential_imag_[index_i];
    Vecd residual_imag = source_current_density_imag_[index_i] - curl_nu_b_imag_[index_i] -
                         sigma_i * electric_potential_gradient_imag_[index_i] -
                         sigma_i * angular_frequency_ * vector_potential_real_[index_i];
    Real coupled_denom =
        operator_diagonal * operator_diagonal + omega_sigma * omega_sigma + TinyReal;
    Vecd delta_real =
        (operator_diagonal * residual_real + omega_sigma * residual_imag) / coupled_denom;
    Vecd delta_imag =
        (operator_diagonal * residual_imag - omega_sigma * residual_real) / coupled_denom;
    Real pair_norm =
        sqrt(delta_real.squaredNorm() + delta_imag.squaredNorm());
    if (!std::isfinite(pair_norm))
    {
        vector_potential_change_rate_real_[index_i] = ZeroData<Vecd>::value;
        vector_potential_change_rate_imag_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    if (pair_norm > max_change_rate_)
    {
        Real scale = max_change_rate_ / (pair_norm + TinyReal);
        delta_real *= scale;
        delta_imag *= scale;
    }
    vector_potential_change_rate_real_[index_i] = delta_real;
    vector_potential_change_rate_imag_[index_i] = delta_imag;
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledEquationInner::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_real_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_real_[index_i];
    vector_potential_imag_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_imag_[index_i];

    if (!std::isfinite(vector_potential_real_[index_i].squaredNorm()))
    {
        vector_potential_real_[index_i] = ZeroData<Vecd>::value;
    }
    if (!std::isfinite(vector_potential_imag_[index_i].squaredNorm()))
    {
        vector_potential_imag_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyEquationComplex::
    VectorPotentialFrequencyEquationComplex(BaseInnerRelation &inner_relation,
                                            BaseContactRelation &contact_relation,
                                            Real angular_frequency,
                                            Real omega_coupling_sign,
                                            const std::string &vector_potential_name,
                                            const std::string &coupled_vector_potential_name,
                                            const std::string &source_current_density_name,
                                            const std::string &electric_potential_gradient_name,
                                            const std::string &curl_nu_b_name,
                                            const std::string &change_rate_name,
                                            Real sigma_relaxation_scaling,
                                            Real sigma_relaxation_floor,
                                            Real magnetic_diagonal_scaling,
                                            Real reference_pseudo_time_step,
                                            Real relaxation_scaling,
                                            Real max_change_rate,
                                            Real balanced_magnetic_diagonal_weight,
                                            Real contact_diagonal_ratio_cap,
                                            Real adaptive_contact_diagonal_cap_strength)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      angular_frequency_(angular_frequency),
      omega_coupling_sign_(omega_coupling_sign),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      balanced_magnetic_diagonal_weight_(balanced_magnetic_diagonal_weight),
      contact_diagonal_ratio_cap_(contact_diagonal_ratio_cap),
      adaptive_contact_diagonal_cap_strength_(adaptive_contact_diagonal_cap_strength),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      coupled_vector_potential_(particles_->getVariableDataByName<Vecd>(coupled_vector_potential_name)),
      source_current_density_(particles_->getVariableDataByName<Vecd>(source_current_density_name)),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>(change_rate_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
    }
}
//=================================================================================================//
void VectorPotentialFrequencyEquationComplex::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_,
                                            contact_configuration_, contact_vol_,
                                            contact_magnetic_reluctivity_,
                                            magnetic_diagonal_scaling_);
    Real jacobi_diagonal =
        ComputeJacobiMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Real balanced_diagonal =
        ComputeBalancedMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                        magnetic_reluctivity_,
                                        contact_configuration_, contact_vol_,
                                        contact_magnetic_reluctivity_,
                                        magnetic_diagonal_scaling_);
    Real magnetic_diagonal =
        BlendMagneticDiagonal(conservative_diagonal, balanced_diagonal,
                             balanced_magnetic_diagonal_weight_);
    Real adaptive_ratio_cap =
        ComputeAdaptiveContactDiagonalRatioCap(
            contact_diagonal_ratio_cap_, conservative_diagonal,
            jacobi_diagonal, adaptive_contact_diagonal_cap_strength_);
    magnetic_diagonal =
        ApplyContactDiagonalRatioCap(magnetic_diagonal, jacobi_diagonal,
                                     adaptive_ratio_cap);
    Vecd rhs = source_current_density_[index_i] - curl_nu_b_[index_i] -
               sigma_i * electric_potential_gradient_[index_i] +
               omega_coupling_sign_ * sigma_i * angular_frequency_ * coupled_vector_potential_[index_i];
    Real operator_diagonal =
        magnetic_diagonal + sigma_relaxation_scaling_ * sigma_i * angular_frequency_ +
        sigma_relaxation_floor_;
    vector_potential_change_rate_[index_i] =
        ClampVectorByNorm(rhs / (operator_diagonal + TinyReal), max_change_rate_);
}
//=================================================================================================//
void VectorPotentialFrequencyEquationComplex::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_[index_i] += pseudo_step_scale *
                                  vector_potential_change_rate_[index_i];
    if (!std::isfinite(vector_potential_[index_i].squaredNorm()))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyCoupledEquationComplex::
    VectorPotentialFrequencyCoupledEquationComplex(BaseInnerRelation &inner_relation,
                                                   BaseContactRelation &contact_relation,
                                                   Real angular_frequency,
                                                   const std::string &vector_potential_real_name,
                                                   const std::string &vector_potential_imag_name,
                                                   const std::string &source_current_density_real_name,
                                                   const std::string &source_current_density_imag_name,
                                                   const std::string &electric_potential_gradient_real_name,
                                                   const std::string &electric_potential_gradient_imag_name,
                                                   const std::string &curl_nu_b_real_name,
                                                   const std::string &curl_nu_b_imag_name,
                                                   const std::string &change_rate_real_name,
                                                   const std::string &change_rate_imag_name,
                                                   Real sigma_relaxation_scaling,
                                                   Real sigma_relaxation_floor,
                                                   Real magnetic_diagonal_scaling,
                                                   Real reference_pseudo_time_step,
                                                   Real relaxation_scaling,
                                                   Real max_change_rate,
                                                   Real balanced_magnetic_diagonal_weight,
                                                   Real contact_diagonal_ratio_cap,
                                                   Real adaptive_contact_diagonal_cap_strength)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      angular_frequency_(angular_frequency),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      balanced_magnetic_diagonal_weight_(balanced_magnetic_diagonal_weight),
      contact_diagonal_ratio_cap_(contact_diagonal_ratio_cap),
      adaptive_contact_diagonal_cap_strength_(adaptive_contact_diagonal_cap_strength),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_real_(particles_->getVariableDataByName<Vecd>(vector_potential_real_name)),
      vector_potential_imag_(particles_->getVariableDataByName<Vecd>(vector_potential_imag_name)),
      source_current_density_real_(particles_->getVariableDataByName<Vecd>(source_current_density_real_name)),
      source_current_density_imag_(particles_->getVariableDataByName<Vecd>(source_current_density_imag_name)),
      electric_potential_gradient_real_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_real_name)),
      electric_potential_gradient_imag_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_imag_name)),
      curl_nu_b_real_(particles_->getVariableDataByName<Vecd>(curl_nu_b_real_name)),
      curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>(curl_nu_b_imag_name)),
      vector_potential_change_rate_real_(particles_->getVariableDataByName<Vecd>(change_rate_real_name)),
      vector_potential_change_rate_imag_(particles_->getVariableDataByName<Vecd>(change_rate_imag_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
    }
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledEquationComplex::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_,
                                            contact_configuration_, contact_vol_,
                                            contact_magnetic_reluctivity_,
                                            magnetic_diagonal_scaling_);
    Real jacobi_diagonal =
        ComputeJacobiMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Real balanced_diagonal =
        ComputeBalancedMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                        magnetic_reluctivity_,
                                        contact_configuration_, contact_vol_,
                                        contact_magnetic_reluctivity_,
                                        magnetic_diagonal_scaling_);
    Real magnetic_diagonal =
        BlendMagneticDiagonal(conservative_diagonal, balanced_diagonal,
                             balanced_magnetic_diagonal_weight_);
    Real adaptive_ratio_cap =
        ComputeAdaptiveContactDiagonalRatioCap(
            contact_diagonal_ratio_cap_, conservative_diagonal,
            jacobi_diagonal, adaptive_contact_diagonal_cap_strength_);
    magnetic_diagonal =
        ApplyContactDiagonalRatioCap(magnetic_diagonal, jacobi_diagonal,
                                     adaptive_ratio_cap);
    Real operator_diagonal =
        magnetic_diagonal + sigma_relaxation_floor_;
    Real omega_sigma = sigma_relaxation_scaling_ * sigma_i * angular_frequency_;

    Vecd residual_real = source_current_density_real_[index_i] - curl_nu_b_real_[index_i] -
                         sigma_i * electric_potential_gradient_real_[index_i] +
                         sigma_i * angular_frequency_ * vector_potential_imag_[index_i];
    Vecd residual_imag = source_current_density_imag_[index_i] - curl_nu_b_imag_[index_i] -
                         sigma_i * electric_potential_gradient_imag_[index_i] -
                         sigma_i * angular_frequency_ * vector_potential_real_[index_i];
    Real coupled_denom =
        operator_diagonal * operator_diagonal + omega_sigma * omega_sigma + TinyReal;
    Vecd delta_real =
        (operator_diagonal * residual_real + omega_sigma * residual_imag) / coupled_denom;
    Vecd delta_imag =
        (operator_diagonal * residual_imag - omega_sigma * residual_real) / coupled_denom;
    Real pair_norm =
        sqrt(delta_real.squaredNorm() + delta_imag.squaredNorm());
    if (!std::isfinite(pair_norm))
    {
        vector_potential_change_rate_real_[index_i] = ZeroData<Vecd>::value;
        vector_potential_change_rate_imag_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    if (pair_norm > max_change_rate_)
    {
        Real scale = max_change_rate_ / (pair_norm + TinyReal);
        delta_real *= scale;
        delta_imag *= scale;
    }
    vector_potential_change_rate_real_[index_i] = delta_real;
    vector_potential_change_rate_imag_[index_i] = delta_imag;
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledEquationComplex::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_real_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_real_[index_i];
    vector_potential_imag_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_imag_[index_i];

    if (!std::isfinite(vector_potential_real_[index_i].squaredNorm()))
    {
        vector_potential_real_[index_i] = ZeroData<Vecd>::value;
    }
    if (!std::isfinite(vector_potential_imag_[index_i].squaredNorm()))
    {
        vector_potential_imag_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyCoupledBlockEquationComplex::
    VectorPotentialFrequencyCoupledBlockEquationComplex(BaseInnerRelation &inner_relation,
                                                        BaseContactRelation &contact_relation,
                                                        Real angular_frequency,
                                                        const std::string &vector_potential_real_name,
                                                        const std::string &vector_potential_imag_name,
                                                        const std::string &source_current_density_real_name,
                                                        const std::string &source_current_density_imag_name,
                                                        const std::string &electric_potential_gradient_real_name,
                                                        const std::string &electric_potential_gradient_imag_name,
                                                        const std::string &curl_nu_b_real_name,
                                                        const std::string &curl_nu_b_imag_name,
                                                        const std::string &change_rate_real_name,
                                                        const std::string &change_rate_imag_name,
                                                        Real sigma_relaxation_scaling,
                                                        Real sigma_relaxation_floor,
                                                        Real magnetic_diagonal_scaling,
                                                        Real reference_pseudo_time_step,
                                                        Real relaxation_scaling,
                                                        Real max_change_rate,
                                                        Real balanced_magnetic_diagonal_weight,
                                                        Real contact_diagonal_ratio_cap,
                                                        Real adaptive_contact_diagonal_cap_strength)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      angular_frequency_(angular_frequency),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      balanced_magnetic_diagonal_weight_(balanced_magnetic_diagonal_weight),
      contact_diagonal_ratio_cap_(contact_diagonal_ratio_cap),
      adaptive_contact_diagonal_cap_strength_(adaptive_contact_diagonal_cap_strength),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_real_(particles_->getVariableDataByName<Vecd>(vector_potential_real_name)),
      vector_potential_imag_(particles_->getVariableDataByName<Vecd>(vector_potential_imag_name)),
      source_current_density_real_(particles_->getVariableDataByName<Vecd>(source_current_density_real_name)),
      source_current_density_imag_(particles_->getVariableDataByName<Vecd>(source_current_density_imag_name)),
      electric_potential_gradient_real_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_real_name)),
      electric_potential_gradient_imag_(particles_->getVariableDataByName<Vecd>(electric_potential_gradient_imag_name)),
      curl_nu_b_real_(particles_->getVariableDataByName<Vecd>(curl_nu_b_real_name)),
      curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>(curl_nu_b_imag_name)),
      vector_potential_change_rate_real_(particles_->getVariableDataByName<Vecd>(change_rate_real_name)),
      vector_potential_change_rate_imag_(particles_->getVariableDataByName<Vecd>(change_rate_imag_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
    }
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledBlockEquationComplex::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_,
                                            contact_configuration_, contact_vol_,
                                            contact_magnetic_reluctivity_,
                                            magnetic_diagonal_scaling_);
    Real jacobi_diagonal =
        ComputeJacobiMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Real balanced_diagonal =
        ComputeBalancedMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                        magnetic_reluctivity_,
                                        contact_configuration_, contact_vol_,
                                        contact_magnetic_reluctivity_,
                                        magnetic_diagonal_scaling_);
    Real conservative_magnetic_diagonal =
        BlendMagneticDiagonal(conservative_diagonal, balanced_diagonal,
                             balanced_magnetic_diagonal_weight_);
    Real adaptive_ratio_cap =
        ComputeAdaptiveContactDiagonalRatioCap(
            contact_diagonal_ratio_cap_, conservative_diagonal,
            jacobi_diagonal, adaptive_contact_diagonal_cap_strength_);
    conservative_magnetic_diagonal =
        ApplyContactDiagonalRatioCap(conservative_magnetic_diagonal, jacobi_diagonal,
                                     adaptive_ratio_cap);
    Matd local_self_block =
        ComputeLocalMagneticSelfBlock(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Matd operator_block =
        local_self_block + sigma_relaxation_floor_ * Matd::Identity();
    Real fallback_operator_diagonal =
        conservative_magnetic_diagonal + sigma_relaxation_floor_;
    Real omega_sigma = sigma_relaxation_scaling_ * sigma_i * angular_frequency_;

    Vecd residual_real = source_current_density_real_[index_i] - curl_nu_b_real_[index_i] -
                         sigma_i * electric_potential_gradient_real_[index_i] +
                         sigma_i * angular_frequency_ * vector_potential_imag_[index_i];
    Vecd residual_imag = source_current_density_imag_[index_i] - curl_nu_b_imag_[index_i] -
                         sigma_i * electric_potential_gradient_imag_[index_i] -
                         sigma_i * angular_frequency_ * vector_potential_real_[index_i];

    Vecd delta_real = ZeroData<Vecd>::value;
    Vecd delta_imag = ZeroData<Vecd>::value;
    SolveCoupledMagneticSelfBlockOrFallback(operator_block,
                                            residual_real,
                                            residual_imag,
                                            omega_sigma,
                                            fallback_operator_diagonal,
                                            delta_real,
                                            delta_imag);
    Real pair_norm =
        sqrt(delta_real.squaredNorm() + delta_imag.squaredNorm());
    if (!std::isfinite(pair_norm))
    {
        vector_potential_change_rate_real_[index_i] = ZeroData<Vecd>::value;
        vector_potential_change_rate_imag_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    if (pair_norm > max_change_rate_)
    {
        Real scale = max_change_rate_ / (pair_norm + TinyReal);
        delta_real *= scale;
        delta_imag *= scale;
    }
    vector_potential_change_rate_real_[index_i] = delta_real;
    vector_potential_change_rate_imag_[index_i] = delta_imag;
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledBlockEquationComplex::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_real_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_real_[index_i];
    vector_potential_imag_[index_i] +=
        pseudo_step_scale * vector_potential_change_rate_imag_[index_i];

    if (!std::isfinite(vector_potential_real_[index_i].squaredNorm()))
    {
        vector_potential_real_[index_i] = ZeroData<Vecd>::value;
    }
    if (!std::isfinite(vector_potential_imag_[index_i].squaredNorm()))
    {
        vector_potential_imag_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyMagneticOnlyEquationComplex::
    VectorPotentialFrequencyMagneticOnlyEquationComplex(BaseInnerRelation &inner_relation,
                                                        BaseContactRelation &contact_relation,
                                                        const std::string &vector_potential_name,
                                                        const std::string &source_current_density_name,
                                                        const std::string &curl_nu_b_name,
                                                        const std::string &change_rate_name,
                                                        Real magnetic_diagonal_scaling,
                                                        Real reference_pseudo_time_step,
                                                        Real relaxation_scaling,
                                                        Real max_change_rate,
                                                        Real balanced_magnetic_diagonal_weight,
                                                        Real contact_diagonal_ratio_cap,
                                                        Real adaptive_contact_diagonal_cap_strength)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      balanced_magnetic_diagonal_weight_(balanced_magnetic_diagonal_weight),
      contact_diagonal_ratio_cap_(contact_diagonal_ratio_cap),
      adaptive_contact_diagonal_cap_strength_(adaptive_contact_diagonal_cap_strength),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      source_current_density_(particles_->getVariableDataByName<Vecd>(source_current_density_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>(change_rate_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
    }
}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyEquationComplex::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_,
                                            contact_configuration_, contact_vol_,
                                            contact_magnetic_reluctivity_,
                                            magnetic_diagonal_scaling_);
    Real jacobi_diagonal =
        ComputeJacobiMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Real balanced_diagonal =
        ComputeBalancedMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                        magnetic_reluctivity_,
                                        contact_configuration_, contact_vol_,
                                        contact_magnetic_reluctivity_,
                                        magnetic_diagonal_scaling_);
    Real magnetic_diagonal =
        BlendMagneticDiagonal(conservative_diagonal, balanced_diagonal,
                             balanced_magnetic_diagonal_weight_);
    Real adaptive_ratio_cap =
        ComputeAdaptiveContactDiagonalRatioCap(
            contact_diagonal_ratio_cap_, conservative_diagonal,
            jacobi_diagonal, adaptive_contact_diagonal_cap_strength_);
    magnetic_diagonal =
        ApplyContactDiagonalRatioCap(magnetic_diagonal, jacobi_diagonal,
                                     adaptive_ratio_cap);
    Vecd residual = source_current_density_[index_i] - curl_nu_b_[index_i];
    vector_potential_change_rate_[index_i] =
        ClampVectorByNorm(residual / (magnetic_diagonal + TinyReal), max_change_rate_);
}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyEquationComplex::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_[index_i] += pseudo_step_scale *
                                  vector_potential_change_rate_[index_i];
    if (!std::isfinite(vector_potential_[index_i].squaredNorm()))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyMagneticOnlyBlockEquationComplex::
    VectorPotentialFrequencyMagneticOnlyBlockEquationComplex(BaseInnerRelation &inner_relation,
                                                             BaseContactRelation &contact_relation,
                                                             const std::string &vector_potential_name,
                                                             const std::string &source_current_density_name,
                                                             const std::string &curl_nu_b_name,
                                                             const std::string &change_rate_name,
                                                             Real magnetic_diagonal_scaling,
                                                             Real reference_pseudo_time_step,
                                                             Real relaxation_scaling,
                                                             Real max_change_rate,
                                                             Real balanced_magnetic_diagonal_weight,
                                                             Real contact_diagonal_ratio_cap,
                                                             Real adaptive_contact_diagonal_cap_strength)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      DataDelegateContact(contact_relation),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      balanced_magnetic_diagonal_weight_(balanced_magnetic_diagonal_weight),
      contact_diagonal_ratio_cap_(contact_diagonal_ratio_cap),
      adaptive_contact_diagonal_cap_strength_(adaptive_contact_diagonal_cap_strength),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      source_current_density_(particles_->getVariableDataByName<Vecd>(source_current_density_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>(change_rate_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
    }
}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyBlockEquationComplex::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_only_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_,
                                            contact_configuration_, contact_vol_,
                                            contact_magnetic_reluctivity_,
                                            magnetic_diagonal_scaling_);
    Real jacobi_diagonal =
        ComputeJacobiMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Real balanced_diagonal =
        ComputeBalancedMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                        magnetic_reluctivity_,
                                        contact_configuration_, contact_vol_,
                                        contact_magnetic_reluctivity_,
                                        magnetic_diagonal_scaling_);
    Real conservative_diagonal =
        BlendMagneticDiagonal(conservative_only_diagonal, balanced_diagonal,
                             balanced_magnetic_diagonal_weight_);
    Real adaptive_ratio_cap =
        ComputeAdaptiveContactDiagonalRatioCap(
            contact_diagonal_ratio_cap_, conservative_only_diagonal,
            jacobi_diagonal, adaptive_contact_diagonal_cap_strength_);
    conservative_diagonal =
        ApplyContactDiagonalRatioCap(conservative_diagonal, jacobi_diagonal,
                                     adaptive_ratio_cap);
    Matd local_self_block =
        ComputeLocalMagneticSelfBlock(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_,
                                      contact_configuration_, contact_vol_,
                                      contact_magnetic_reluctivity_,
                                      magnetic_diagonal_scaling_);
    Vecd residual = source_current_density_[index_i] - curl_nu_b_[index_i];
    Vecd delta =
        SolveMagneticSelfBlockOrFallback(local_self_block, residual, conservative_diagonal);
    vector_potential_change_rate_[index_i] =
        ClampVectorByNorm(delta, max_change_rate_);
}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyBlockEquationComplex::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_[index_i] += pseudo_step_scale *
                                  vector_potential_change_rate_[index_i];
    if (!std::isfinite(vector_potential_[index_i].squaredNorm()))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialFrequencyMagneticOnlyEquationInner::
    VectorPotentialFrequencyMagneticOnlyEquationInner(BaseInnerRelation &inner_relation,
                                                      const std::string &vector_potential_name,
                                                      const std::string &source_current_density_name,
                                                      const std::string &curl_nu_b_name,
                                                      const std::string &change_rate_name,
                                                      Real magnetic_diagonal_scaling,
                                                      Real reference_pseudo_time_step,
                                                      Real relaxation_scaling,
                                                      Real max_change_rate)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      reference_pseudo_time_step_(reference_pseudo_time_step),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      source_current_density_(particles_->getVariableDataByName<Vecd>(source_current_density_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>(change_rate_name)) {}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyEquationInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real magnetic_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_, magnetic_diagonal_scaling_);
    Vecd residual = source_current_density_[index_i] - curl_nu_b_[index_i];
    vector_potential_change_rate_[index_i] =
        ClampVectorByNorm(residual / (magnetic_diagonal + TinyReal), max_change_rate_);
}
//=================================================================================================//
void VectorPotentialFrequencyMagneticOnlyEquationInner::update(size_t index_i, Real dt)
{
    Real pseudo_step_scale =
        relaxation_scaling_ *
        ComputeNormalizedPseudoTimeStep(dt, reference_pseudo_time_step_);
    vector_potential_[index_i] += pseudo_step_scale *
                                  vector_potential_change_rate_[index_i];
    if (!std::isfinite(vector_potential_[index_i].squaredNorm()))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
FrequencyAEquationResidualDiagnostic::
    FrequencyAEquationResidualDiagnostic(SPHBody &sph_body,
                                         Real angular_frequency,
                                         Real omega_coupling_sign,
                                         const std::string &a_name,
                                         const std::string &coupled_a_name,
                                         const std::string &source_name,
                                         const std::string &grad_phi_name,
                                         const std::string &curl_nu_b_name,
                                         const std::string &residual_name,
                                         const std::string &relative_residual_name)
    : LocalDynamics(sph_body),
      angular_frequency_(angular_frequency),
      omega_coupling_sign_(omega_coupling_sign),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      a_(particles_->getVariableDataByName<Vecd>(a_name)),
      coupled_a_(particles_->getVariableDataByName<Vecd>(coupled_a_name)),
      source_(particles_->getVariableDataByName<Vecd>(source_name)),
      grad_phi_(particles_->getVariableDataByName<Vecd>(grad_phi_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(curl_nu_b_name)),
      residual_(particles_->registerStateVariableData<Vecd>(residual_name)),
      relative_residual_(particles_->registerStateVariableData<Real>(relative_residual_name))
{
    particles_->addVariableToWrite<Vecd>(residual_name);
    particles_->addVariableToWrite<Real>(relative_residual_name);
}
//=================================================================================================//
void FrequencyAEquationResidualDiagnostic::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = electrical_conductivity_[index_i];
    const Vecd omega_term = omega_coupling_sign_ * sigma_i * angular_frequency_ * coupled_a_[index_i];

    residual_[index_i] = source_[index_i] - curl_nu_b_[index_i] - sigma_i * grad_phi_[index_i] + omega_term;

    const Real numerator = residual_[index_i].norm();
    const Real denominator =
        source_[index_i].norm() + curl_nu_b_[index_i].norm() +
        (sigma_i * grad_phi_[index_i]).norm() + omega_term.norm() + TinyReal;
    relative_residual_[index_i] = numerator / denominator;
}
//=================================================================================================//
FrequencyElectricFieldCurrentAndJouleHeatInner::
    FrequencyElectricFieldCurrentAndJouleHeatInner(BaseInnerRelation &inner_relation,
                                                   Real angular_frequency,
                                                   const std::string &a_real_name,
                                                   const std::string &a_imag_name,
                                                   const std::string &grad_phi_real_name,
                                                   const std::string &grad_phi_imag_name,
                                                   const std::string &e_real_name,
                                                   const std::string &e_imag_name,
                                                   const std::string &j_real_name,
                                                   const std::string &j_imag_name,
                                                   const std::string &joule_name,
                                                   const std::string &temperature_change_rate_name)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      angular_frequency_(angular_frequency),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      rho_cp_(particles_->getVariableDataByName<Real>("RhoCp")),
      a_real_(particles_->getVariableDataByName<Vecd>(a_real_name)),
      a_imag_(particles_->getVariableDataByName<Vecd>(a_imag_name)),
      grad_phi_real_(particles_->getVariableDataByName<Vecd>(grad_phi_real_name)),
      grad_phi_imag_(particles_->getVariableDataByName<Vecd>(grad_phi_imag_name)),
      electric_field_real_(particles_->getVariableDataByName<Vecd>(e_real_name)),
      electric_field_imag_(particles_->getVariableDataByName<Vecd>(e_imag_name)),
      current_density_real_(particles_->getVariableDataByName<Vecd>(j_real_name)),
      current_density_imag_(particles_->getVariableDataByName<Vecd>(j_imag_name)),
      joule_heat_source_(particles_->getVariableDataByName<Real>(joule_name)),
      temperature_change_rate_by_joule_(particles_->getVariableDataByName<Real>(temperature_change_rate_name)) {}
//=================================================================================================//
void FrequencyElectricFieldCurrentAndJouleHeatInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Vecd e_real = angular_frequency_ * a_imag_[index_i] - grad_phi_real_[index_i];
    Vecd e_imag = -angular_frequency_ * a_real_[index_i] - grad_phi_imag_[index_i];
    if (!std::isfinite(e_real.squaredNorm()) || !std::isfinite(e_imag.squaredNorm()))
    {
        electric_field_real_[index_i] = ZeroData<Vecd>::value;
        electric_field_imag_[index_i] = ZeroData<Vecd>::value;
        current_density_real_[index_i] = ZeroData<Vecd>::value;
        current_density_imag_[index_i] = ZeroData<Vecd>::value;
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }

    electric_field_real_[index_i] = e_real;
    electric_field_imag_[index_i] = e_imag;

    Vecd j_real = electrical_conductivity_[index_i] * e_real;
    Vecd j_imag = electrical_conductivity_[index_i] * e_imag;
    if (!std::isfinite(j_real.squaredNorm()) || !std::isfinite(j_imag.squaredNorm()))
    {
        current_density_real_[index_i] = ZeroData<Vecd>::value;
        current_density_imag_[index_i] = ZeroData<Vecd>::value;
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }

    current_density_real_[index_i] = j_real;
    current_density_imag_[index_i] = j_imag;

    Real joule_heat = 0.5 * (j_real.dot(e_real) + j_imag.dot(e_imag));
    if (!std::isfinite(joule_heat))
    {
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }

    joule_heat_source_[index_i] = joule_heat;
    temperature_change_rate_by_joule_[index_i] = joule_heat / (rho_cp_[index_i] + TinyReal);
}
//=================================================================================================//
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_TEAM7_APHI_FREQUENCY_DYNAMICS_HPP
