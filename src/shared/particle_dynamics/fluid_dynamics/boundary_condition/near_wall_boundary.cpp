#include "near_wall_boundary.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
NearWallDistance::NearWallDistance(BaseContactRelation &wall_contact_relation)
    : LocalDynamics(wall_contact_relation.getSPHBody()), DataDelegateContact(wall_contact_relation),
      spacing_ref_(sph_body_.getSPHAdaptation().ReferenceSpacing()),
      distance_default_(100.0 * spacing_ref_),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_pos_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Position"));
        wall_n_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        wall_phi_.push_back(contact_particles_[k]->getVariableDataByName<Real>("SignedDistance"));
    }
}
//=================================================================================================//
void NearWallDistance::evaluateDistanceAndNormal(size_t index_i, Vecd &distance, Vecd &normal)
{
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *pos_k = wall_pos_[k];
        Vecd *n_k = wall_n_[k];
        Real *phi_k = wall_phi_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd temp = (pos_[index_i] - pos_k[index_j]) + phi_k[index_j] * n_k[index_j];
            if (temp.squaredNorm() < distance.squaredNorm())
            {
                distance = temp;       // more reliable distance
                normal = n_k[index_j]; // more reliable normal
            }
        }
    }
}
//=================================================================================================//
DistanceFromWall::DistanceFromWall(BaseContactRelation &wall_contact_relation)
    : NearWallDistance(wall_contact_relation),
      distance_from_wall_(particles_->registerStateVariable<Vecd>("DistanceFromWall")) {}
//=================================================================================================//
void DistanceFromWall::interaction(size_t index_i, Real dt)
{
    Vecd distance = distance_default_ * Vecd::Ones();
    Vecd normal = Vecd::Ones();
    NearWallDistance::evaluateDistanceAndNormal(index_i, distance, normal);
    // prediction with regularization
    Vecd normal_distance = distance.dot(normal) * normal;
    Real limiter = SMIN(3.0 * (distance - normal_distance).norm() / spacing_ref_, 1.0);
    distance_from_wall_[index_i] = (1.0 - limiter) * normal_distance + limiter * distance;
}
//=================================================================================================//
BoundingFromWall::BoundingFromWall(BaseContactRelation &wall_contact_relation)
    : NearWallDistance(wall_contact_relation),
      distance_min_(0.25 * spacing_ref_) {}
//=================================================================================================//
void BoundingFromWall::interaction(size_t index_i, Real dt)
{
    Vecd distance = distance_default_ * Vecd::Ones();
    Vecd normal = Vecd::Ones();
    NearWallDistance::evaluateDistanceAndNormal(index_i, distance, normal);
    // bounding if the particle cross the wall
    Real projection = distance.dot(normal);
    if (projection < distance_min_)
    {
        pos_[index_i] += 0.5 * spacing_ref_ * normal - distance; // flip near wall distance
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
