/**
 * @file 	general_dynamics_refinement.cpp
 * @author	Yijie Sun and Xiangyu Hu
 */

#include "general_dynamics_refinement.h"
#include "base_particles.hpp"
#include "level_set_shape.h"

namespace SPH
{
//================================================================================================ =//
Vecd ComputeDensityErrorInner::
    getPositionFromDensityError(const StdVec<size_t> &original_indices, const StdVec<Vecd> &initial_new_positions,
                                const StdVec<size_t> &new_indices, Real min_distance, Real max_distance)
{
    StdVec<Vecd> new_positions = initial_new_positions;
    size_t index_center = new_indices[0];
    Vecd pos_final = new_positions[0];
    Vecd pos_iteration = new_positions[0];
    Real residual = 1.0e-4;
    int iteration = 0;
    Real error_sum_min = rho0_ * 10.0;
    Real temp_E = 0.0;

    while (iteration < 100)
    {
        temp_E = ABS(sigma_E_);
        pos_iteration = getPosition(original_indices, new_positions, new_indices);
        Vecd displacement = pos_iteration - particles_->pos_[index_center];
        pos_iteration = particles_->pos_[index_center] + positionLimitation(displacement, min_distance, max_distance);

        if (ABS(ABS(sigma_E_) - temp_E) < residual)
        {
            pos_final = pos_iteration;
            break;
        }
        if (error_sum_min > sigma_E_)
        {
            error_sum_min = ABS(sigma_E_);
            pos_final = pos_iteration;
        }
        iteration += 1;
        new_positions[0] = pos_iteration;
        new_positions[1] = 2.0 * particles_->pos_[index_center] - pos_iteration;
    }
    for (size_t n = 0; n != new_indices.size(); ++n)
        density_error_[new_indices[n]] = error_sum_min;

    return pos_final;
}
//=================================================================================================//
Vecd ComputeDensityErrorInner::positionLimitation(Vecd displacement, Real min_distance, Real max_distance)
{
    Vecd modify_displacement = displacement;
    if (displacement.norm() > max_distance)
        modify_displacement = displacement * max_distance / (displacement.norm() + TinyReal);
    if (displacement.norm() < min_distance)
        modify_displacement = displacement * min_distance / (displacement.norm() + TinyReal);

    return modify_displacement;
}
//=================================================================================================//
Vecd ComputeDensityErrorInner::
    getPosition(const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions, const StdVec<size_t> &new_indices)
{
    E_cof_sigma_ = 0.0;
    sigma_E_ = 0.0;
    E_cof_ = Vecd::Zero();
    densityErrorOfNewGeneratedParticles(new_indices, new_positions);
    densityErrorOfNeighborParticles(new_indices, original_indices, new_positions);

    E_cof_sigma_ += E_cof_.dot(E_cof_);
    Real cof = 1.0 / (E_cof_sigma_ + TinyReal);
    Real k = sigma_E_ * cof;
    Vecd dr = E_cof_ * k;
    Vecd update_position = new_positions[0] + dr;

    return update_position;
}
//=================================================================================================//
void ComputeDensityErrorInner::
    densityErrorOfNewGeneratedParticles(const StdVec<size_t> &new_indices, const StdVec<Vecd> &new_positions)
{
    sign_new_indices_.clear();
    size_t index_rho = new_indices[0];
    Vecd grad_kernel = computeKernelGradient(index_rho);
    for (size_t n = 0; n != new_positions.size(); ++n)
    {
        Vecd displacement = new_positions[n] - particles_->pos_[index_rho];
        Real rho_newIndex = computeNewGeneratedParticleDensity(index_rho, new_positions[n]);
        Real error = particles_->rho_[index_rho] - rho_newIndex + grad_kernel.dot(-displacement);

        sign_new_indices_.push_back(error / (ABS(error) + TinyReal));
        sigma_E_ += sign_new_indices_[n] * error;
    }
    E_cof_ += sign_new_indices_[0] * 2.0 * dW_new_indices_[0] - sign_new_indices_[1] * 2.0 * dW_new_indices_[1];
    E_cof_ += sign_new_indices_[0] * grad_new_indices_[0] - sign_new_indices_[1] * grad_new_indices_[1];
}
//=================================================================================================//
Vecd ComputeDensityErrorInner::computeKernelGradient(size_t index_rho)
{
    Neighborhood &inner_neighborhood = inner_configuration_[index_rho];
    Vecd grad_kernel = Vecd::Zero();
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        grad_kernel += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] *
                       particles_->mass_[inner_neighborhood.j_[n]] / rho0_ /
                       particles_->Vol_[inner_neighborhood.j_[n]];
    }
    return grad_kernel;
}
//=================================================================================================//
Real ComputeDensityErrorInner::
    computeKernelWeightBetweenParticles(Real h_ratio, Vecd displacement, Real Vol_ratio)
{
    Real distance = displacement.norm();
    Kernel *kernel_ptr_ = particle_adaptation_.getKernel();
    Real cutoff_radius = kernel_ptr_->CutOffRadius(h_ratio);
    Real kernel_weight = 0.0;
    if (distance <= cutoff_radius)
        kernel_weight = kernel_ptr_->W(h_ratio, distance, displacement) * Vol_ratio;

    return kernel_weight;
}
//=================================================================================================//
Vecd ComputeDensityErrorInner::
    computeKernelWeightGradientBetweenParticles(Real h_ratio_min, Vecd displacement, Real Vol)
{
    Real distance = displacement.norm();
    Vecd e_ij = displacement / (distance + TinyReal);
    Kernel *kernel_ptr_ = particle_adaptation_.getKernel();
    Real cutoff_radius = kernel_ptr_->CutOffRadius(h_ratio_min);
    Vecd grad_kernel = Vecd::Zero();
    if (distance <= cutoff_radius)
    {
        Real d_weight = kernel_ptr_->dW(h_ratio_min, distance, displacement) * Vol;
        grad_kernel = d_weight * e_ij;
    }
    return grad_kernel;
}
//=================================================================================================//
Real ComputeDensityErrorInner::
    computeNewGeneratedParticleDensity(size_t index_rho, const Vecd &position)
{
    Real Vol_newIndex = particles_->Vol_[index_rho] / 2.0;
    Real h_newIndex = pow(particles_->Vol_[index_rho] / Vol_newIndex, 1.0 / (Real)Dimensions);

    Real W0 = particle_adaptation_.getKernel()->W0(h_newIndex, ZeroVecd);
    Real inv_sigma_i = inv_sigma0_ / particle_adaptation_.NumberDensityScaleFactor(h_newIndex);
    Real sigma_newIndex = W0;

    Vecd displacement = 2.0 * (position - particles_->pos_[index_rho]);
    dW_new_indices_.push_back(computeKernelWeightGradientBetweenParticles(h_newIndex, displacement, Vol_newIndex));
    sigma_newIndex += computeKernelWeightBetweenParticles(h_newIndex, displacement);

    Vecd grad_sigma = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_rho];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Vecd displacement = position - particles_->pos_[inner_neighborhood.j_[n]];
        Real h_ratio_min = SMIN(h_ratio_[inner_neighborhood.j_[n]], h_newIndex);
        Real Vol_j = particles_->mass_[inner_neighborhood.j_[n]] / rho0_;
        Real Vol_ratio = Vol_j / Vol_newIndex;
        sigma_newIndex += computeKernelWeightBetweenParticles(h_newIndex, displacement, Vol_ratio);
        grad_sigma += computeKernelWeightGradientBetweenParticles(h_ratio_min, displacement, Vol_j);
    }
    grad_new_indices_.push_back(grad_sigma);
    sigma_newIndex = sigma_newIndex * rho0_ * inv_sigma_i;
    return sigma_newIndex;
}
//=================================================================================================//
void ComputeDensityErrorInner::
    computeDensityErrorOnNeighborParticles(Neighborhood &neighborhood, size_t index_rho,
                                           const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions)
{
    Real Vol_newIndex = particles_->Vol_[index_rho] / 2.0;
    Real h_newIndex = pow(particles_->Vol_[index_rho] / Vol_newIndex, 1.0 / (Real)Dimensions);

    for (size_t k = 0; k != neighborhood.current_size_; ++k)
    {
        Real h_ratio_j = h_ratio_[neighborhood.j_[k]];
        Vecd pos_j = particles_->pos_[neighborhood.j_[k]];
        Real Vol_j = particles_->mass_[neighborhood.j_[k]] / rho0_;
        Real inv_sigma_j = inv_sigma0_ / particle_adaptation_.NumberDensityScaleFactor(h_newIndex);
        Real sigma_split_j = 0.0;

        for (size_t n = 0; n != original_indices.size(); ++n)
        {
            Vecd displacement = pos_j - particles_->pos_[original_indices[n]];
            Real Vol_ratio = particles_->Vol_[original_indices[n]] / Vol_j;
            sigma_split_j += computeKernelWeightBetweenParticles(h_ratio_j, displacement, Vol_ratio);
        }

        StdVec<Vecd> grad_sigma_j;
        Vecd sigma_co_j = Vecd::Zero();
        for (size_t n = 0; n != new_positions.size(); ++n)
        {
            Real h_ratio_min = SMIN(h_ratio_j, h_newIndex);
            Vecd displacement = pos_j - new_positions[n];
            Real Vol_ratio = Vol_newIndex / Vol_j;
            sigma_split_j -= computeKernelWeightBetweenParticles(h_ratio_j, displacement, Vol_ratio);
            Vecd grad_j = computeKernelWeightGradientBetweenParticles(h_ratio_min, -displacement, Vol_j);
            sigma_co_j -= grad_j * sign_new_indices_[n];
            grad_sigma_j.push_back(computeKernelWeightGradientBetweenParticles(h_ratio_min, displacement, Vol_newIndex));
        }

        sigma_split_j = sigma_split_j * rho0_ * inv_sigma_j;

        Real sign = sigma_split_j / (ABS(sigma_split_j) + TinyReal);
        sigma_co_j += sign * (grad_sigma_j[0] + grad_sigma_j[1]);

        sigma_E_ += sign * sigma_split_j;
        E_cof_ += sign * (grad_sigma_j[1] - grad_sigma_j[0]);
        E_cof_sigma_ += sigma_co_j.dot(sigma_co_j);
    }
}
//=================================================================================================//
void ComputeDensityErrorInner::
    densityErrorOfNeighborParticles(const StdVec<size_t> &new_indices,
                                    const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions)
{
    Neighborhood &neighborhood = inner_configuration_[new_indices[0]];
    computeDensityErrorOnNeighborParticles(neighborhood, new_indices[0], original_indices, new_positions);
}
//================================================================================================ =//
void ComputeDensityErrorInner::initializeDensityError()
{
    for (size_t index_i = 0; index_i != particles_->total_real_particles_; ++index_i)
        density_error_[index_i] = 0.0;
}
//=================================================================================================//
Vecd ComputeDensityErrorWithWall::computeKernelGradient(size_t index_rho)
{
    Vecd grad_kernel = ComputeDensityErrorInner::computeKernelGradient(index_rho);
    Neighborhood &contact_neighborhood = inner_configuration_[index_rho];
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            grad_kernel += contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
    }
    return grad_kernel;
}
//=================================================================================================//
void ComputeDensityErrorWithWall::
    densityErrorOfNeighborParticles(const StdVec<size_t> &new_indices,
                                    const StdVec<size_t> &original_indices, const StdVec<Vecd> &new_positions)
{
    ComputeDensityErrorInner::densityErrorOfNeighborParticles(new_indices, original_indices, new_positions);
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Neighborhood &neighborhood = (*contact_configuration_[k])[new_indices[0]];
        computeDensityErrorOnNeighborParticles(neighborhood, new_indices[0], original_indices, new_positions);
    }
}
//=================================================================================================//
Real ComputeDensityErrorWithWall::
    computeNewGeneratedParticleDensity(size_t index_rho, const Vecd &position)
{
    Real Vol_newIndex = particles_->Vol_[index_rho] / 2.0;
    Real h_newIndex = pow(particles_->Vol_[index_rho] / Vol_newIndex, 1.0 / (Real)Dimensions);

    Real inv_sigma_i = inv_sigma0_ / particle_adaptation_.NumberDensityScaleFactor(h_newIndex);
    ;
    Real sigma_inner = ComputeDensityErrorInner::computeNewGeneratedParticleDensity(index_rho, position);
    Vecd grad_sigma = Vecd::Zero();
    Real sigma_newIndex = 0.0;
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_rho];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Vecd displacement = position - particles_->pos_[contact_neighborhood.j_[n]];
            Real h_ratio_min = SMIN(h_ratio_[contact_neighborhood.j_[n]], h_newIndex);
            Real Vol_j = particles_->mass_[contact_neighborhood.j_[n]] / rho0_;
            Real Vol_ratio = Vol_j / Vol_newIndex;
            sigma_newIndex += computeKernelWeightBetweenParticles(h_newIndex, displacement, Vol_ratio);
            grad_sigma += computeKernelWeightGradientBetweenParticles(h_ratio_min, displacement, Vol_j);
        }
    }
    grad_new_indices_[grad_new_indices_.size() - 1] += grad_sigma;
    sigma_newIndex = sigma_newIndex * rho0_ * inv_sigma_i + sigma_inner;
    return sigma_newIndex;
}
//=================================================================================================//
ParticleRefinementWithPrescribedArea::
    ParticleRefinementWithPrescribedArea(SPHBody &sph_body, Shape &refinement_region)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      refinement_region_bounds_(refinement_region.getBounds()),
      particle_adaptation_(DynamicCast<ParticleSplitAndMerge>(this, *sph_body.sph_adaptation_)),
      inv_rho0_(1.0 / sph_body_.base_material_->ReferenceDensity()),
      Vol_(particles_->Vol_), pos_(particles_->pos_), rho_(particles_->rho_),
      mass_(particles_->mass_),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
bool ParticleRefinementWithPrescribedArea::
    checkLocation(const BoundingBox &refinement_region_bounds, Vecd position, Real volume)
{
    int bound_number = 0;
    for (int axis_direction = 0; axis_direction != Dimensions; ++axis_direction)
    {
        Real particle_spacing = pow(volume, 1.0 / (Real)Dimensions);
        if (position[axis_direction] > (refinement_region_bounds.first_[axis_direction] + particle_spacing) &&
            position[axis_direction] < (refinement_region_bounds.second_[axis_direction] - particle_spacing))
            bound_number += 1;
    }
    return bound_number != Dimensions ? false : true;
}
//=================================================================================================//
ParticleSplitWithPrescribedArea::
    ParticleSplitWithPrescribedArea(SPHBody &sph_body, Shape &refinement_region, size_t body_buffer_width)
    : ParticleRefinementWithPrescribedArea(sph_body, refinement_region)
{
    particles_->addBufferParticles(body_buffer_width);
    sph_body.allocateConfigurationMemoriesForBufferParticles();
}
//=================================================================================================//
bool ParticleSplitWithPrescribedArea::splitCriteria(size_t index_i)
{
    Real non_deformed_volume = mass_[index_i] * inv_rho0_;
    Vecd &position = pos_[index_i];
    bool is_split_allowed = particle_adaptation_.isSplitAllowed(non_deformed_volume);
    bool is_split_inside = checkLocation(refinement_region_bounds_, position, non_deformed_volume);

    return (is_split_allowed && is_split_inside) ? true : false;
}
//=================================================================================================//
void ParticleSplitWithPrescribedArea::splittingModel(size_t index_i, StdVec<size_t> &new_indices)
{
    size_t particle_real_number = particles_->total_real_particles_ + particle_number_change;
    particles_->copyFromAnotherParticle(particle_real_number, index_i);
    new_indices.push_back(index_i);
    new_indices.push_back(particle_real_number);
    Vecd pos_splitting = getSplittingPosition(new_indices);
    updateNewlySplittingParticle(index_i, particle_real_number, pos_splitting);

    particle_number_change += 1;
    split_position_.push_back(2.0 * pos_[index_i] - pos_splitting);
    split_index_.push_back(std::make_pair(index_i, particle_real_number));
}
//=================================================================================================//
void ParticleSplitWithPrescribedArea::setupDynamics(Real dt)
{
    split_index_.clear();
    split_position_.clear();
    particle_number_change = 0;
}
//=================================================================================================//
void ParticleSplitWithPrescribedArea::update(size_t index_i, Real dt)
{
    for (size_t num = 0; num != particle_number_change; ++num)
    {
        if (index_i == split_index_[num].first)
        {
            size_t index_i = split_index_[num].first;
            size_t index_j = split_index_[num].second;
            pos_[index_i] = split_position_[num];
            mass_[index_i] = mass_[index_j];
            Vol_[index_i] = Vol_[index_j];
            h_ratio_[index_i] = h_ratio_[index_j];

            particles_->total_real_particles_ += 1;
            if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
            {
                std::cout << "ParticleSplitWithPrescribedArea: \n"
                          << "Not enough body buffer particles! Exit the code."
                          << "\n";
                exit(0);
            }
        }
    }
}
//=================================================================================================//
Vecd ParticleSplitWithPrescribedArea::getSplittingPosition(const StdVec<size_t> &new_indices)
{
    srand(int(new_indices[0]));
    Real delta_random = 0 + 2.0 * Pi * rand() / RAND_MAX * (2.0 * Pi - 0);
    Vecd pos_ = particles_->pos_[new_indices[0]];
    Real delta = delta_random + Pi;
    Real Vol_split = Vol_[new_indices[0]] / 2.0;
    Real particle_spacing_j = pow(Vol_split, 1.0 / (Real)Dimensions);

    return particle_adaptation_.splittingPattern(pos_, particle_spacing_j, delta);
}
//=================================================================================================//
void ParticleSplitWithPrescribedArea::
    updateNewlySplittingParticle(size_t index_center, size_t index_new, Vecd pos_split)
{
    mass_[index_new] = 0.5 * mass_[index_center];
    Vol_[index_new] = mass_[index_new] * inv_rho0_;
    Real particle_spacing_j = pow(Vol_[index_new], 1.0 / (Real)Dimensions);
    h_ratio_[index_new] = sph_body_.sph_adaptation_->ReferenceSpacing() / particle_spacing_j;
    particles_->pos_[index_new] = pos_split;
}
//=================================================================================================//
Vecd SplitWithMinimumDensityErrorInner::getSplittingPosition(const StdVec<size_t> &new_indices)
{
    StdVec<size_t> original_indices;
    StdVec<Vecd> new_positions;
    size_t index_i = new_indices[0];
    original_indices.push_back(index_i);
    Real particle_spacing = pow(Vol_[index_i] / 2.0, 1.0 / (Real)Dimensions);
    Vecd pos_j = particle_adaptation_.splittingPattern(pos_[index_i], particle_spacing, 0.25 * Pi);
    Vecd pos_i = 2.0 * pos_[index_i] - pos_j;
    new_positions.push_back(pos_i);
    new_positions.push_back(pos_j);

    Vecd position = compute_density_error.getPositionFromDensityError(
        original_indices, new_positions, new_indices, 0.2 * particle_spacing, 0.65 * particle_spacing);
    for (size_t n = 0; n != new_indices.size(); ++n)
        total_split_error_[new_indices[n]] = compute_density_error.density_error_[new_indices[n]];
    return position;
}
//=================================================================================================//
void SplitWithMinimumDensityErrorInner::setupDynamics(Real dt)
{
    ParticleSplitWithPrescribedArea::setupDynamics(dt);
    compute_density_error.initializeDensityError();
}
//=================================================================================================//
void SplitWithMinimumDensityErrorInner::update(size_t index_i, Real dt)
{
    total_split_error_[index_i] = compute_density_error.density_error_[index_i];
    for (size_t num = 0; num != particle_number_change; ++num)
    {
        if (index_i == split_index_[num].first)
        {
            size_t index_j = split_index_[num].second;
            total_split_error_[index_j] = compute_density_error.density_error_[index_j];
        }
    }
    ParticleSplitWithPrescribedArea::update(index_i, dt);
}
//=================================================================================================//
ParticleMergeWithPrescribedArea::
    ParticleMergeWithPrescribedArea(BaseInnerRelation &inner_relation, Shape &refinement_region)
    : ParticleRefinementWithPrescribedArea(inner_relation.getSPHBody(), refinement_region),
      DataDelegateInner<BaseParticles, DataDelegateEmptyBase>(inner_relation),
      all_particle_data_(particles_->getAllParticleData()), vel_n_(particles_->vel_)
{
    particles_->registerVariable(total_merge_error_, "MergeDensityError", Real(0));
}
//=================================================================================================//
void ParticleMergeWithPrescribedArea::setupDynamics(Real dt)
{
    tag_merged_.clear();
    for (size_t index_i = 0; index_i != particles_->total_real_particles_; ++index_i)
    {
        tag_merged_.push_back(false);
    }
}
//=================================================================================================//
bool ParticleMergeWithPrescribedArea::mergeCriteria(size_t index_i, StdVec<size_t> &merge_indices)
{
    Real non_deformed_volume = mass_[index_i] * inv_rho0_;
    bool resolution_check = particle_adaptation_.mergeResolutionCheck(non_deformed_volume);
    Real particle_spacing = pow(non_deformed_volume, 1.0 / (Real)Dimensions);
    Real search_threshold = 1.2;
    Real search_distance = search_threshold * particle_spacing;
    if (resolution_check)
    {
        bool neighbor_check = findMergeParticles(index_i, merge_indices, particle_spacing, search_distance);
        if (neighbor_check)
        {
            Vecd &position = pos_[index_i];
            bool location_check = !checkLocation(refinement_region_bounds_, position, non_deformed_volume);
            if (location_check)
                return true;
        }
    }
    return false;
}
//=================================================================================================//
bool ParticleMergeWithPrescribedArea::
    findMergeParticles(size_t index_i, StdVec<size_t> &merge_indices, Real search_size, Real search_distance)
{
    Vecd &position = pos_[index_i];
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real distance_ = (particles_->pos_[index_j] - position).norm();
        Real Vol_j = mass_[index_j] * inv_rho0_;
        if (!tag_merged_[index_j])
            if (distance_ < search_distance)
            {
                bool resolution_check = particle_adaptation_.mergeResolutionCheck(Vol_j);
                if (ABS(search_size - particle_adaptation_.ReferenceSpacing()) < TinyReal)
                {
                    if (!resolution_check)
                    {
                        merge_indices.push_back(index_j);
                        tag_merged_[index_j] = true;
                        return true;
                    }
                }
                else if (resolution_check)
                {
                    merge_indices.push_back(index_j);
                    tag_merged_[index_j] = true;
                    return true;
                }
            }
    }
    return false;
}
//=================================================================================================//
void ParticleMergeWithPrescribedArea::mergingModel(const StdVec<size_t> &merge_indices)
{
    updateMergedParticleInformation(particles_->total_real_particles_, merge_indices);
    size_t merge_change_number = merge_indices.size() - 1;
    for (size_t k = merge_change_number; k != size_t(-1); --k)
        particles_->copyFromAnotherParticle(merge_indices[k], particles_->total_real_particles_ - (merge_change_number - k));
    particles_->total_real_particles_ -= merge_change_number;
}
//=================================================================================================//
void ParticleMergeWithPrescribedArea::
    updateMergedParticleInformation(size_t merged_index, const StdVec<size_t> &merge_indices)
{
    StdVec<Real> merge_mass_;
    Real total_mass = 0.0;
    for (size_t k = 0; k != merge_indices.size(); ++k)
    {
        merge_mass_.push_back(mass_[merge_indices[k]]);
        total_mass += mass_[merge_indices[k]];
    }
    merge_particle_value_(all_particle_data_, merged_index, merge_indices, merge_mass_);
    mass_[merged_index] = total_mass;
    Vol_[merged_index] = mass_[merged_index] * inv_rho0_;
    Real particle_spacing = pow(Vol_[merged_index], 1.0 / (Real)Dimensions);
    h_ratio_[merged_index] = sph_body_.sph_adaptation_->ReferenceSpacing() / particle_spacing;
}
//=================================================================================================//
void MergeWithMinimumDensityErrorInner::setupDynamics(Real dt)
{
    compute_density_error.tag_split_.clear();
    tag_merged_.clear();
    for (size_t index_i = 0; index_i != particles_->total_real_particles_; ++index_i)
    {
        tag_merged_.push_back(false);
        compute_density_error.tag_split_.push_back(false);
    }
    compute_density_error.initializeDensityError();
}
//=================================================================================================//
bool MergeWithMinimumDensityErrorInner::mergeCriteria(size_t index_i, StdVec<size_t> &merge_indices)
{
    Real non_deformed_volume = mass_[index_i] * inv_rho0_;
    bool resolution_check = particle_adaptation_.mergeResolutionCheck(non_deformed_volume);
    Real particle_spacing_small = pow(non_deformed_volume, 1.0 / (Real)Dimensions);
    Real particle_spacing_large = pow(non_deformed_volume * 2.0, 1.0 / (Real)Dimensions);
    Real search_threshold = 1.2;
    Real search_distance_small = search_threshold * particle_spacing_small;
    Real search_distance_large = search_threshold * particle_spacing_large;
    if (resolution_check)
    {
        bool neighbor_check_large = findMergeParticles(index_i, merge_indices, particle_spacing_large, search_distance_large);
        bool neighbor_check_small = findMergeParticles(index_i, merge_indices, particle_spacing_small, search_distance_small);
        if (neighbor_check_large && neighbor_check_small)
        {
            Vecd &position = pos_[index_i];
            bool location_check = !checkLocation(refinement_region_bounds_, position, non_deformed_volume);
            if (location_check)
                return true;
        }
    }
    return false;
}
//=================================================================================================//
void MergeWithMinimumDensityErrorInner::mergingModel(const StdVec<size_t> &merge_indices)
{
    StdVec<size_t> new_indices; // the first and third particles
    new_indices.push_back(merge_indices[0]);
    new_indices.push_back(merge_indices[merge_indices.size() - 1]);

    size_t temporary_index = particles_->total_real_particles_; // the first buffer particle as temporary
    updateMergedParticleInformation(temporary_index, merge_indices);
    Vecd pos_merging = getMergingPosition(new_indices, merge_indices);
    updateNewlyMergingParticle(temporary_index, new_indices, pos_merging);
    kineticEnergyConservation(merge_indices);

    // the follow for deleting the second particle
    particles_->copyFromAnotherParticle(merge_indices[1], particles_->total_real_particles_ - 1);
    particles_->total_real_particles_ -= 1;
}
//=================================================================================================//
Vecd MergeWithMinimumDensityErrorInner::
    getMergingPosition(const StdVec<size_t> &new_indices, const StdVec<size_t> &merge_indices)
{
    StdVec<Vecd> initial_new_positions;
    size_t index_center = new_indices[0];
    Real particle_spacing = pow(Vol_[index_center] / 2.0, 1.0 / (Real)Dimensions);
    Real distance_min = 0.2 * particle_spacing; // angularMomentumConservation(index_center, merge_indices);
    Real distance_max = 0.65 * particle_spacing;

    Vecd pos_j = pos_[merge_indices[0]]; // pos_[index_center] + distance_min * position / (position.norm() + TinyReal);
    Vecd pos_i = 2.0 * pos_[index_center] - pos_j;
    initial_new_positions.push_back(pos_i);
    initial_new_positions.push_back(pos_j);
    compute_density_error.tag_split_[new_indices[0]] = true;
    compute_density_error.tag_split_[new_indices[1]] = true;

    Vecd position_final = compute_density_error.getPositionFromDensityError(
        merge_indices, initial_new_positions, new_indices, distance_min, distance_max);
    for (size_t n = 0; n != new_indices.size(); ++n)
        total_merge_error_[new_indices[n]] = compute_density_error.density_error_[new_indices[n]];
    return position_final;
}
//=================================================================================================//
void MergeWithMinimumDensityErrorInner::updateNewlyMergingParticle(size_t index_center, const StdVec<size_t> &new_indices, Vecd pos_split)
{
    for (size_t n = 0; n != new_indices.size(); ++n)
    {
        mass_[new_indices[n]] = 0.5 * mass_[index_center];
        Vol_[new_indices[n]] = mass_[new_indices[n]] * inv_rho0_;
        Real particle_spacing_j = pow(Vol_[new_indices[n]], 1.0 / (Real)Dimensions);
        h_ratio_[new_indices[n]] = sph_body_.sph_adaptation_->ReferenceSpacing() / particle_spacing_j;
    }
    pos_[new_indices[0]] = pos_split;
    pos_[new_indices[1]] = 2.0 * pos_[index_center] - pos_[new_indices[0]];
}
//=================================================================================================//
Real MergeWithMinimumDensityErrorInner::angularMomentumConservation(size_t index_center, const StdVec<size_t> &merge_indices)
{
    rotation = 0.0;
    Real mass = 0.0;
    Real vel_square = 0.0;
    for (size_t n = 0; n != merge_indices.size(); ++n)
    {
        Vecd pos = pos_[merge_indices[n]] - pos_[index_center];
        Vecd vel = vel_n_[merge_indices[n]] - vel_n_[index_center];
        rotation += mass_[merge_indices[n]] * (pos[0] * vel[1] - pos[1] * vel[0]);
        mass += mass_[merge_indices[n]];
        vel_square += mass_[merge_indices[n]] / mass_[index_center] * vel_n_[merge_indices[n]].squaredNorm();
    }
    rotation = rotation / mass;
    Real E = vel_square - vel_n_[index_center].squaredNorm();
    Real distance_min = sqrt((rotation) * (rotation) / (ABS(E) + TinyReal));

    return distance_min;
}
//=================================================================================================//
void MergeWithMinimumDensityErrorInner::kineticEnergyConservation(const StdVec<size_t> &merge_indices)
{
    Real E_total = 0.5 * (2.0 * vel_n_[merge_indices[0]].squaredNorm() +
                          vel_n_[merge_indices[1]].squaredNorm() + vel_n_[merge_indices[2]].squaredNorm());
    Vecd linear_m = 0.5 * (2.0 * vel_n_[merge_indices[0]] + vel_n_[merge_indices[1]] + vel_n_[merge_indices[2]]);
    Real angular_m = rotation * 2.0;
    if (ABS(pos_[merge_indices[0]][1]) > TinyReal && ABS(pos_[merge_indices[0]][0]) > TinyReal)
    {
        Real cof1 = pos_[merge_indices[0]][1] / (pos_[merge_indices[0]][0] + TinyReal);
        Real cof2 = (-linear_m[0] * pos_[merge_indices[0]][1] + linear_m[1] * pos_[merge_indices[0]][0] + angular_m) /
                    (2.0 * pos_[merge_indices[0]][0] + TinyReal);
        Real a = 2.0 * (cof1 * cof1 + 1.0);
        Real b = 4.0 * cof1 * cof2 - 2.0 * linear_m[0] - 2.0 * cof1 * linear_m[1];
        Real c = 2.0 * cof2 * cof2 - E_total - 2.0 * cof2 * linear_m[1] + linear_m.squaredNorm();
        if ((b * b - 4.0 * a * c) >= 0.0)
            vel_n_[merge_indices[0]][0] = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a + TinyReal);
        else
        {
            vel_n_[merge_indices[0]][0] = -b / (2.0 * a + TinyReal);
        }
        vel_n_[merge_indices[0]][1] = cof1 * vel_n_[merge_indices[0]][0] + cof2;
    }
    else if (ABS(pos_[merge_indices[0]][1]) <= TinyReal) // y=0.0
    {
        vel_n_[merge_indices[0]][1] = 0.5 * (angular_m / (pos_[merge_indices[0]][0] + TinyReal) + linear_m[1]);
        Real a = 2.0;
        Real b = -2.0 * linear_m[0];
        Real c = -E_total + linear_m.squaredNorm() - 2.0 * vel_n_[merge_indices[0]][1] * linear_m[1] +
                 2.0 * vel_n_[merge_indices[0]][1] * vel_n_[merge_indices[0]][1];
        if ((b * b - 4.0 * a * c) >= 0.0)
            vel_n_[merge_indices[0]][0] = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a + TinyReal);
        else
        {
            vel_n_[merge_indices[0]][0] = -b / (2.0 * a + TinyReal);
        }
    }
    else if (ABS(pos_[merge_indices[0]][0]) <= TinyReal) // x=0.0
    {
        vel_n_[merge_indices[0]][0] = 0.5 * (-angular_m / (pos_[merge_indices[0]][1] + TinyReal) + linear_m[0]);
        Real a = 2.0;
        Real b = -2.0 * linear_m[1];
        Real c = -E_total + linear_m.squaredNorm() - 2.0 * vel_n_[merge_indices[0]][0] * linear_m[0] +
                 2.0 * vel_n_[merge_indices[0]][0] * vel_n_[merge_indices[0]][0];
        if ((b * b - 4.0 * a * c) >= 0.0)
            vel_n_[merge_indices[0]][1] = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a + TinyReal);
        else
        {
            vel_n_[merge_indices[0]][1] = -b / (2.0 * a + TinyReal);
        }
    }
    vel_n_[merge_indices[1]][0] = linear_m[0] - vel_n_[merge_indices[0]][0];
    vel_n_[merge_indices[1]][1] = linear_m[1] - vel_n_[merge_indices[0]][1];
}
//=================================================================================================//
} // namespace SPH
