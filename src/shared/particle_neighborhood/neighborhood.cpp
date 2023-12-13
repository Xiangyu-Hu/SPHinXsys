/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "neighborhood.h"

#include "all_complex_bodies.h"
#include "base_particle_dynamics.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void Neighborhood::removeANeighbor(size_t neighbor_n)
{
    current_size_--;
    j_[neighbor_n] = j_[current_size_];
    W_ij_[neighbor_n] = W_ij_[current_size_];
    dW_ijV_j_[neighbor_n] = dW_ijV_j_[current_size_];
    r_ij_[neighbor_n] = r_ij_[current_size_];
    e_ij_[neighbor_n] = e_ij_[current_size_];
}
//=================================================================================================//
void NeighborBuilder::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j, const Real &Vol_j)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
    neighborhood.dW_ijV_j_.push_back(kernel_->dW(distance, displacement) * Vol_j);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(kernel_->e(distance, displacement));
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilder::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                         const Vecd &displacement, size_t index_j, const Real &Vol_j)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
    neighborhood.dW_ijV_j_[current_size] = kernel_->dW(distance, displacement) * Vol_j;
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = kernel_->e(distance, displacement);
}
//=================================================================================================//
void NeighborBuilder::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j, const Real &Vol_j,
                                     Real i_h_ratio, Real h_ratio_min)
{
    neighborhood.j_.push_back(index_j);
    Real weight = distance < kernel_->CutOffRadius(i_h_ratio) ? kernel_->W(i_h_ratio, distance, displacement) : 0.0;
    neighborhood.W_ij_.push_back(weight);
    neighborhood.dW_ijV_j_.push_back(kernel_->dW(h_ratio_min, distance, displacement) * Vol_j);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilder::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                         const Vecd &displacement, size_t index_j, const Real &Vol_j,
                                         Real i_h_ratio, Real h_ratio_min)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = distance < kernel_->CutOffRadius(i_h_ratio)
                                           ? kernel_->W(i_h_ratio, distance, displacement)
                                           : 0.0;
    neighborhood.dW_ijV_j_[current_size] = kernel_->dW(h_ratio_min, distance, displacement) * Vol_j;
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
}
//=================================================================================================//
NeighborBuilderInner::NeighborBuilderInner(SPHBody &body) : NeighborBuilder()
{
    kernel_ = body.sph_adaptation_->getKernel();
}
//=================================================================================================//
void NeighborBuilderInner::operator()(Neighborhood &neighborhood,
                                      const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance_metric = displacement.squaredNorm();
    if (kernel_->checkIfWithinCutOffRadius(displacement) && index_i != index_j)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, std::sqrt(distance_metric), displacement, index_j, std::get<2>(list_data_j))
            : initializeNeighbor(neighborhood, std::sqrt(distance_metric), displacement, index_j, std::get<2>(list_data_j));
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderInnerAdaptive::
    NeighborBuilderInnerAdaptive(SPHBody &body)
    : NeighborBuilder(),
      h_ratio_(*body.getBaseParticles().getVariableByName<Real>("SmoothingLengthRatio"))
{
    kernel_ = body.sph_adaptation_->getKernel();
}
//=================================================================================================//
void NeighborBuilderInnerAdaptive::
operator()(Neighborhood &neighborhood, const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance = displacement.norm();
    Real i_h_ratio = h_ratio_[index_i];
    Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[index_j]);
    Real cutoff_radius = kernel_->CutOffRadius(h_ratio_min);
    if (distance < cutoff_radius && index_i != index_j)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j), i_h_ratio, h_ratio_min)
            : initializeNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j), i_h_ratio, h_ratio_min);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderSelfContact::
    NeighborBuilderSelfContact(SPHBody &body)
    : NeighborBuilder(),
      pos0_(*body.getBaseParticles().getVariableByName<Vecd>("InitialPosition"))
{
    kernel_ = body.sph_adaptation_->getKernel();
}
//=================================================================================================//
void NeighborBuilderSelfContact::operator()(Neighborhood &neighborhood,
                                            const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance = displacement.norm();
    Real distance0 = (pos0_[index_i] - pos0_[index_j]).norm();
    if (distance < kernel_->CutOffRadius() && distance0 > kernel_->CutOffRadius())
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j))
            : initializeNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j));
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderContact::
    NeighborBuilderContact(SPHBody &body, SPHBody &contact_body) : NeighborBuilder()
{
    Kernel *source_kernel = body.sph_adaptation_->getKernel();
    Kernel *target_kernel = contact_body.sph_adaptation_->getKernel();
    kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
}
//=================================================================================================//
void NeighborBuilderContact::operator()(Neighborhood &neighborhood,
                                        const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance = displacement.norm();
    if (distance < kernel_->CutOffRadius())
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j))
            : initializeNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j));
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderSurfaceContact::NeighborBuilderSurfaceContact(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilderContact(body, contact_body)
{
    Real source_smoothing_length = body.sph_adaptation_->ReferenceSmoothingLength();
    Real target_smoothing_length = contact_body.sph_adaptation_->ReferenceSmoothingLength();
    kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(0.5 * (source_smoothing_length + target_smoothing_length));
}
//=================================================================================================//
NeighborBuilderContactBodyPart::
    NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part) : NeighborBuilder()
{
    contact_body_part.getSPHBody().getBaseParticles().registerVariable(part_indicator_, "BodyPartByParticleIndicator");
    Kernel *source_kernel = body.sph_adaptation_->getKernel();
    Kernel *target_kernel = contact_body_part.getSPHBody().sph_adaptation_->getKernel();
    kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;

    BodyPartByParticle &contact_body_part_by_particle = DynamicCast<BodyPartByParticle>(this, contact_body_part);
    IndexVector part_particles = contact_body_part_by_particle.body_part_particles_;

    for (size_t i = 0; i != part_particles.size(); ++i)
    {
        part_indicator_[part_particles[i]] = 1;
    }
}
//=================================================================================================//
void NeighborBuilderContactBodyPart::operator()(Neighborhood &neighborhood,
                                                const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance = displacement.norm();
    if (distance < kernel_->CutOffRadius() && part_indicator_[index_j] == 1)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j))
            : initializeNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j));
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
NeighborBuilderContactAdaptive::
    NeighborBuilderContactAdaptive(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilder(), adaptation_(*body.sph_adaptation_), contact_adaptation_(*contact_body.sph_adaptation_),
      relative_h_ref_(adaptation_.ReferenceSmoothingLength() / contact_adaptation_.ReferenceSmoothingLength())
{
    kernel_ = adaptation_.getKernel();
}
//=================================================================================================//
void NeighborBuilderContactAdaptive::operator()(Neighborhood &neighborhood,
                                                const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance_metric = displacement.squaredNorm();
    Real i_h_ratio = adaptation_.SmoothingLengthRatio(index_i);
    Real h_ratio_min = SMIN(i_h_ratio, relative_h_ref_ * contact_adaptation_.SmoothingLengthRatio(index_j));
    if (distance_metric < kernel_->CutOffRadiusSqr(h_ratio_min))
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, std::sqrt(distance_metric),
                             displacement, index_j, std::get<2>(list_data_j), i_h_ratio, h_ratio_min)
            : initializeNeighbor(neighborhood, std::sqrt(distance_metric),
                                 displacement, index_j, std::get<2>(list_data_j), i_h_ratio, h_ratio_min);
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
NeighborBuilderContactShell::NeighborBuilderContactShell(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilder(),
      n_(*contact_body.getBaseParticles().getVariableByName<Vecd>("NormalDirection")),
      thickness_(*contact_body.getBaseParticles().getVariableByName<Real>("Thickness")),
      mean_curvature_(*contact_body.getBaseParticles().getVariableByName<Real>("MeanCurvature")),
      particle_distance_(contact_body.getSPHBodyResolutionRef())
{
    Kernel *source_kernel = body.sph_adaptation_->getKernel();
    Kernel *target_kernel = contact_body.sph_adaptation_->getKernel();
    kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
}
//=================================================================================================//
void NeighborBuilderContactShell::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                 size_t index_j, const Real &W_ij, const Real &dW_ijV_j, const Vecd &e_ij)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(W_ij);
    neighborhood.dW_ijV_j_.push_back(dW_ijV_j);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(e_ij);
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilderContactShell::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                     size_t index_j, const Real &W_ij, const Real &dW_ijV_j, const Vecd &e_ij)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = W_ij;
    neighborhood.dW_ijV_j_[current_size] = dW_ijV_j;
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = e_ij;
}
//=================================================================================================//
void NeighborBuilderContactShell::operator()(Neighborhood &neighborhood,
                                             const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    const Vecd n_j = n_[index_j]; // normal direction of shell particle in cell linked list with index j

    const Vecd pos_j = std::get<1>(list_data_j);
    const Vecd displacement = pos_i - pos_j;
    const Real distance = displacement.norm();

    const Real Vol_j = std::get<2>(list_data_j);

    // correct normal direction, make sure it points from fluid to shell
    const Real direction_corrector = -displacement.dot(n_j) >= 0 ? 1 : -1;
    const Vecd n_j_corrected = direction_corrector * n_j; // sign of r_ij and n_j
    // H is the total mean curvature
    // for 2D, curvature=H
    // for 3D, mean curvature=H/2
    const Real H_j_corrected = direction_corrector * mean_curvature_[index_j] / (Dimensions - 1); // mean curvature with corrected sign

    if (distance < kernel_->CutOffRadius())
    {
        Real W_ijV_j_ttl = kernel_->W(distance, displacement) * Vol_j;
        Real dW_ijV_j_ttl = kernel_->dW(distance, displacement) * Vol_j;
        Vecd dW_ijV_j_e_ij_ttl = dW_ijV_j_ttl * displacement / (distance + TinyReal);

        Real vol_ratio = H_j_corrected >= 0 ? 1 + 1 / (H_j_corrected + TinyReal) / particle_distance_
                                            : 1 + 1 / (H_j_corrected - TinyReal) / particle_distance_; // 1+R/dp
        Real Vol_j_dummy = Vol_j;
        Vecd pos_j_dummy = pos_j + n_j_corrected * particle_distance_;
        Vecd displacement_dummy = pos_i - pos_j_dummy;
        Real distance_dummy = displacement_dummy.norm();

        while (distance_dummy < kernel_->CutOffRadius())
        {
            const Real dVol = vol_ratio >= 0 ? 1 + 1.0 / (vol_ratio + TinyReal) : 1 + 1.0 / (vol_ratio - TinyReal);
            Vol_j_dummy *= std::pow(dVol, Dimensions - 1);
            Real dW_ijV_j = kernel_->dW(distance_dummy, displacement_dummy) * Vol_j_dummy;
            Vecd e_ij = displacement_dummy / (distance_dummy + TinyReal);
            W_ijV_j_ttl += kernel_->W(distance_dummy, displacement_dummy) * Vol_j_dummy;
            dW_ijV_j_ttl += dW_ijV_j;
            dW_ijV_j_e_ij_ttl += dW_ijV_j * e_ij;

            // calculate the position and volume of the next dummy particle
            vol_ratio += 1.0;
            pos_j_dummy += n_j_corrected * particle_distance_;
            displacement_dummy = pos_i - pos_j_dummy;
            distance_dummy = displacement_dummy.norm();
        }

        Vecd e_ij_corrected = dW_ijV_j_e_ij_ttl / dW_ijV_j_ttl;
        Real W_ij_corrected = W_ijV_j_ttl / Vol_j * particle_distance_ / thickness_[index_j];
        Real dW_ijV_j_corrected = dW_ijV_j_ttl * particle_distance_ / thickness_[index_j];

        // create new neighborhood
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected)
            : initializeNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
