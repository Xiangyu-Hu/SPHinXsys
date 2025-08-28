#include "eulerian_ghost_boundary.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
GhostCreationInESPH::GhostCreationInESPH(BaseInnerRelation &inner_relation, Ghost<ReserveSizeFactor> &ghost_boundary)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      ghost_boundary_(ghost_boundary), get_inner_neighbor_(&inner_relation.getSPHBody()),
      indicator_(particles_->getVariableDataByName<int>("Indicator")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      ghost_bound_(ghost_boundary.GhostBound())
{
    ghostGenerationAndAddToConfiguration();
};
//=================================================================================================//
void GhostCreationInESPH::ghostGenerationAndAddToConfiguration()
{
    ghost_bound_.second = ghost_bound_.first;
    RealAndGhostParticleData real_and_ghost_necessary_data;

    for (size_t index_i = 0; index_i != particles_->TotalRealParticles(); ++index_i)
    {
        if (indicator_[index_i] == 1)
        {
            Vecd gradient_summation = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
                gradient_summation += dW_ijV_j * inner_neighborhood.e_ij_[n];
            }
            Real ghost_particle_dW_ijV_j = -gradient_summation.norm();
            Vecd ghost_particle_eij = -gradient_summation / ghost_particle_dW_ijV_j;
            Real distance_to_ghost = fabs(sph_body_.getInitialShape().findSignedDistance(pos_[index_i]));
            Vecd displacement_to_ghost = distance_to_ghost * ghost_particle_eij;
            Vecd ghost_position = pos_[index_i] - displacement_to_ghost;
            mutex_create_ghost_particle_.lock();
            size_t ghost_particle_index = ghost_bound_.second;
            ghost_bound_.second++;
            ghost_boundary_.checkWithinGhostSize(ghost_bound_);
            particles_->updateGhostParticle(ghost_particle_index, index_i);
            // Here, we ensure the volume as the same the real particle
            Vol_[ghost_particle_index] = Vol_[index_i];
            pos_[ghost_particle_index] = ghost_position;
            mutex_create_ghost_particle_.unlock();
            Real double_distance_to_ghost = 2.0 * distance_to_ghost;
            Real ghost_particle_dW_ij = ghost_particle_dW_ijV_j / (Vol_[ghost_particle_index] + TinyReal);
            get_inner_neighbor_(inner_neighborhood, double_distance_to_ghost, ghost_particle_dW_ij, ghost_particle_eij, ghost_particle_index);

            // add all necessary data of real particle and corresponding ghost particle
            real_and_ghost_necessary_data.real_index_ = index_i;
            real_and_ghost_necessary_data.ghost_index_ = ghost_particle_index;
            real_and_ghost_necessary_data.e_ij_ghost_ = ghost_particle_eij;
            real_and_ghost_particle_data_.push_back(real_and_ghost_necessary_data);
        }
    }
}
//=================================================================================================//
GhostBoundaryConditionSetupInESPH::
    GhostBoundaryConditionSetupInESPH(BaseInnerRelation &inner_relation, GhostCreationInESPH &ghost_creation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      mom_(particles_->getVariableDataByName<Vecd>("Momentum")),
      ghost_bound_(ghost_creation.ghost_bound_),
      real_and_ghost_particle_data_(ghost_creation.real_and_ghost_particle_data_),
      boundary_type_(particles_->registerStateVariable<int>("BoundaryType")),
      W0_(sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd))
{
    setupBoundaryTypes();
}
//=================================================================================================//
void GhostBoundaryConditionSetupInESPH::resetBoundaryConditions()
{
    for (size_t ghost_number = 0; ghost_number != real_and_ghost_particle_data_.size(); ++ghost_number)
    {
        size_t index_i = real_and_ghost_particle_data_[ghost_number].real_index_;
        size_t ghost_index = real_and_ghost_particle_data_[ghost_number].ghost_index_;
        Vecd e_ig = real_and_ghost_particle_data_[ghost_number].e_ij_ghost_;

        // Dispatch the appropriate boundary condition
        switch (boundary_type_[index_i])
        {
        case 3: // this refer to the different types of wall boundary conditions
            applyNonSlipWallBoundary(ghost_index, index_i);
            applyReflectiveWallBoundary(ghost_index, index_i, e_ig);
            break;
        case 4:
            applyTopBoundary(ghost_index, index_i);
            break;
        case 5:
            applyPressureOutletBC(ghost_index, index_i);
            break;
        case 7:
            applySymmetryBoundary(ghost_index, index_i, e_ig);
            break;
        case 9:
            applyFarFieldBoundary(ghost_index, index_i);
            break;
        case 10:
            applyGivenValueInletFlow(ghost_index);
            applyVelocityInletFlow(ghost_index, index_i);
            break;
        case 36:
            applyOutletBoundary(ghost_index, index_i);
            break;
        }
    }
}
//=================================================================================================//
GhostKernelGradientUpdate::GhostKernelGradientUpdate(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      kernel_gradient_original_summation_(
          particles_->registerStateVariable<Vecd>("KernelGradientOriginalSummation")),
      indicator_(particles_->getVariableDataByName<int>("Indicator")) {};
//=================================================================================================//
void GhostKernelGradientUpdate::interaction(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1)
    {
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vecd gradient_summation = Vecd::Zero();
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
            if (index_j < particles_->TotalRealParticles())
            {
                gradient_summation += dW_ijV_j * e_ij;
            }
        }
        kernel_gradient_original_summation_[index_i] = gradient_summation;
    }
}
//=================================================================================================//
void GhostKernelGradientUpdate::update(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1)
    {
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (index_j >= particles_->TotalRealParticles())
            {
                Vecd ghost_dW_ijV_j_with_e_ij = -kernel_gradient_original_summation_[index_i];
                Real ghost_dW_ijV_j = -fabs(ghost_dW_ijV_j_with_e_ij.norm());
                inner_neighborhood.dW_ij_[n] = ghost_dW_ijV_j / Vol_[index_j];
                inner_neighborhood.e_ij_[n] = ghost_dW_ijV_j_with_e_ij / (ghost_dW_ijV_j + TinyReal);
            }
        }
    }
}
//=================================================================================================//
} // namespace SPH
