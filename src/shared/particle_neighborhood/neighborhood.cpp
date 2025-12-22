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
    dW_ij_[neighbor_n] = dW_ij_[current_size_];
    r_ij_[neighbor_n] = r_ij_[current_size_];
    e_ij_[neighbor_n] = e_ij_[current_size_];
}
//=================================================================================================//
void NeighborBuilder::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
    neighborhood.dW_ij_.push_back(kernel_->dW(distance, displacement));
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(kernel_->e(distance, displacement));
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilder::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                         const Vecd &displacement, size_t index_j)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
    neighborhood.dW_ij_[current_size] = kernel_->dW(distance, displacement);
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = kernel_->e(distance, displacement);
}
//=================================================================================================//
void NeighborBuilder::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j,
                                     Real i_h_ratio, Real h_ratio_min)
{
    neighborhood.j_.push_back(index_j);
    Real weight = distance < kernel_->CutOffRadius(i_h_ratio) ? kernel_->W(i_h_ratio, distance, displacement) : 0.0;
    neighborhood.W_ij_.push_back(weight);
    neighborhood.dW_ij_.push_back(kernel_->dW(h_ratio_min, distance, displacement));
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilder::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                         const Vecd &displacement, size_t index_j,
                                         Real i_h_ratio, Real h_ratio_min)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = distance < kernel_->CutOffRadius(i_h_ratio)
                                           ? kernel_->W(i_h_ratio, distance, displacement)
                                           : 0.0;
    neighborhood.dW_ij_[current_size] = kernel_->dW(h_ratio_min, distance, displacement);
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
}
//=================================================================================================//
Kernel *NeighborBuilder::chooseKernel(SPHBody &body, SPHBody &target_body)
{
    Kernel *kernel = body.getSPHAdaptation().getKernel();
    Kernel *target_kernel = target_body.getSPHAdaptation().getKernel();
    return kernel->SmoothingLength() > target_kernel->SmoothingLength() ? kernel : target_kernel;
}
//=================================================================================================//
NeighborBuilderInner::NeighborBuilderInner(SPHBody &body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()) {}
//=================================================================================================//
void NeighborBuilderInner::operator()(Neighborhood &neighborhood,
                                      const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance_metric = displacement.squaredNorm();
    if (kernel_->checkIfWithinCutOffRadius(displacement) && index_i != index_j)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, std::sqrt(distance_metric), displacement, index_j)
            : initializeNeighbor(neighborhood, std::sqrt(distance_metric), displacement, index_j);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderInnerAdaptive::
    NeighborBuilderInnerAdaptive(SPHBody &body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()),
      h_ratio_(body.getBaseParticles().getVariableDataByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void NeighborBuilderInnerAdaptive::
operator()(Neighborhood &neighborhood, const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    Real i_h_ratio = h_ratio_[index_i];
    Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[index_j]);
    Real cutoff_radius = kernel_->CutOffRadius(h_ratio_min);
    if (distance < cutoff_radius && index_i != index_j)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, i_h_ratio, h_ratio_min)
            : initializeNeighbor(neighborhood, distance, displacement, index_j, i_h_ratio, h_ratio_min);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderSelfContact::
    NeighborBuilderSelfContact(SPHBody &body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()),
      pos0_(body.getBaseParticles().registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")) {}
//=================================================================================================//
void NeighborBuilderSelfContact::operator()(Neighborhood &neighborhood,
                                            const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    Real distance0 = (pos0_[index_i] - pos0_[index_j]).norm();
    if (distance < kernel_->CutOffRadius() && distance0 > kernel_->CutOffRadius())
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderContact::NeighborBuilderContact(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilder(NeighborBuilder::chooseKernel(body, contact_body)) {}
//=================================================================================================//
void NeighborBuilderContact::operator()(Neighborhood &neighborhood,
                                        const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    if (kernel_->checkIfWithinCutOffRadius(displacement))
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderSurfaceContact::NeighborBuilderSurfaceContact(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilderContact(body, contact_body)
{
    Real source_smoothing_length = body.getSPHAdaptation().ReferenceSmoothingLength();
    Real target_smoothing_length = contact_body.getSPHAdaptation().ReferenceSmoothingLength();
    kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(0.5 * (source_smoothing_length + target_smoothing_length));
}
//=================================================================================================//
NeighborBuilderContactBodyPart::NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part)
    : NeighborBuilder(NeighborBuilder::chooseKernel(body, contact_body_part.getSPHBody())),
      part_indicator_(body.getBaseParticles().registerStateVariableData<int>("BodyPartByParticleIndicator"))
{
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
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    if (distance < kernel_->CutOffRadius() && part_indicator_[index_j] == 1)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j);
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
NeighborBuilderContactAdaptive::NeighborBuilderContactAdaptive(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()), adaptation_(body.getSPHAdaptation()),
      contact_adaptation_(contact_body.getSPHAdaptation()),
      relative_h_ref_(adaptation_.ReferenceSmoothingLength() /
                      contact_adaptation_.ReferenceSmoothingLength()) {}
//=================================================================================================//
void NeighborBuilderContactAdaptive::operator()(Neighborhood &neighborhood,
                                                const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance_metric = displacement.squaredNorm();
    Real i_h_ratio = adaptation_.SmoothingLengthRatio(index_i);
    Real h_ratio_min = SMIN(i_h_ratio, relative_h_ref_ * contact_adaptation_.SmoothingLengthRatio(index_j));
    if (distance_metric < kernel_->CutOffRadiusSqr(h_ratio_min))
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, std::sqrt(distance_metric),
                             displacement, index_j, i_h_ratio, h_ratio_min)
            : initializeNeighbor(neighborhood, std::sqrt(distance_metric),
                                 displacement, index_j, i_h_ratio, h_ratio_min);
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
BaseNeighborBuilderContactShell::BaseNeighborBuilderContactShell(SPHBody &shell_body)
    : NeighborBuilder(shell_body.getSPHAdaptation().getKernel()),
      n_(shell_body.getBaseParticles().getVariableDataByName<Vecd>("NormalDirection")),
      thickness_(shell_body.getBaseParticles().getVariableDataByName<Real>("Thickness")),
      k1_ave_(shell_body.getBaseParticles().registerStateVariableData<Real>("Average1stPrincipleCurvature")),
      k2_ave_(shell_body.getBaseParticles().registerStateVariableData<Real>("Average2ndPrincipleCurvature")),
      particle_distance_(shell_body.getSPHBodyResolutionRef()) {}
//=================================================================================================//
void BaseNeighborBuilderContactShell::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                     size_t index_j, const Real &W_ij, const Real &dW_ij, const Vecd &e_ij)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(W_ij);
    neighborhood.dW_ij_.push_back(dW_ij);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(e_ij);
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void BaseNeighborBuilderContactShell::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                         size_t index_j, const Real &W_ij, const Real &dW_ij, const Vecd &e_ij)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = W_ij;
    neighborhood.dW_ij_[current_size] = dW_ij;
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = e_ij;
}
//=================================================================================================//
BaseNeighborBuilderContactFromShell::BaseNeighborBuilderContactFromShell(SPHBody &body, SPHBody &contact_body, bool normal_correction)
    : BaseNeighborBuilderContactShell(contact_body),
      direction_corrector_(normal_correction ? -1 : 1)
{
}
//=================================================================================================//
void BaseNeighborBuilderContactFromShell::update_neighbors(Neighborhood &neighborhood,
                                                           const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;

    const Vecd pos_j = list_data_j.second;
    const Vecd displacement = pos_i - pos_j;
    const Real distance = displacement.norm();

    if (distance < kernel_->CutOffRadius())
    {
        Real W_ijV_j_ttl = kernel_->W(distance, displacement);
        Real dW_ijV_j_ttl = kernel_->dW(distance, displacement);
        Vecd dW_ijV_j_e_ij_ttl = dW_ijV_j_ttl * displacement / (distance + TinyReal);

        // correct normal direction pointing from fluid to shell
        const Vecd n_j = direction_corrector_ * n_[index_j]; // normal direction of shell particle in cell linked list with index j

        Vecd pos_j_dummy = pos_j + n_j * particle_distance_;
        Vecd displacement_dummy = pos_i - pos_j_dummy;
        Real distance_dummy = displacement_dummy.norm();

        // k1 and k2 are principle curvatures of the shell
        // For 2d, k1 is the curvature, k2 = 0
        const Real k1_j = direction_corrector_ * k1_ave_[index_j];
        const Real k2_j = direction_corrector_ * k2_ave_[index_j];

        int counter = 0;
        while (distance_dummy < kernel_->CutOffRadius())
        {
            counter++;
            const Real factor_1 = 1 + counter * k1_j * particle_distance_;
            const Real factor_2 = 1 + counter * k2_j * particle_distance_;
            if (factor_1 <= 0 || factor_2 <= 0)
                break;
            const Real Vol_j_dummy = factor_1 * factor_2;
            Real dW_ijV_j = kernel_->dW(distance_dummy, displacement_dummy) * Vol_j_dummy;
            Vecd e_ij = displacement_dummy / distance_dummy;
            W_ijV_j_ttl += kernel_->W(distance_dummy, displacement_dummy) * Vol_j_dummy;
            dW_ijV_j_ttl += dW_ijV_j;
            dW_ijV_j_e_ij_ttl += dW_ijV_j * e_ij;

            // calculate the position and volume of the next dummy particle
            pos_j_dummy += n_j * particle_distance_;
            displacement_dummy = pos_i - pos_j_dummy;
            distance_dummy = displacement_dummy.norm();
        }
        Vecd e_ij_corrected = dW_ijV_j_e_ij_ttl / dW_ijV_j_ttl;
        Real W_ij_corrected = W_ijV_j_ttl * particle_distance_ / thickness_[index_j]; // from surface area to volume
        Real dW_ijV_j_corrected = dW_ijV_j_ttl * particle_distance_;                  // from surface area to volume

        // create new neighborhood
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected)
            : initializeNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderContactFromShellToFluid::NeighborBuilderContactFromShellToFluid(SPHBody &body, SPHBody &contact_body, bool normal_correction)
    : BaseNeighborBuilderContactFromShell(body, contact_body, normal_correction)
{
    Real fluid_reference_spacing = body.getSPHAdaptation().ReferenceSpacing();
    Real shell_reference_spacing = contact_body.getSPHAdaptation().ReferenceSpacing();
    if (fluid_reference_spacing < shell_reference_spacing)
        throw std::runtime_error("NeighborBuilderContactToShell: fluid spacing should be larger or equal than shell spacing...");
    kernel_ = body.getSPHAdaptation().getKernel();
}
//=================================================================================================//
NeighborBuilderContactFromFluidToShell::NeighborBuilderContactFromFluidToShell(SPHBody &body, SPHBody &contact_body, bool normal_correction)
    : BaseNeighborBuilderContactShell(body),
      direction_corrector_(normal_correction ? -1 : 1)
{
    kernel_ = NeighborBuilder::chooseKernel(body, contact_body);
}
//=================================================================================================//
void NeighborBuilderContactFromFluidToShell::operator()(Neighborhood &neighborhood,
                                                        const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;

    const Vecd pos_j = list_data_j.second;
    const Vecd displacement = pos_i - pos_j;
    const Real distance = displacement.norm();

    if (distance < kernel_->CutOffRadius())
    {
        // W_ij is not used in force from fluid classes, no need to modify
        const Real W_ij = kernel_->W(distance, displacement);
        Real dW_ijV_j_ttl = kernel_->dW(distance, displacement);
        Vecd dW_ijV_j_e_ij_ttl = dW_ijV_j_ttl * displacement / (distance + TinyReal);

        // correct normal direction, pointing from fluid to shell
        const Vecd n_i = direction_corrector_ * n_[index_i];

        Vecd pos_i_dummy = pos_i + n_i * particle_distance_;
        Vecd displacement_dummy = pos_i_dummy - pos_j;
        Real distance_dummy = displacement_dummy.norm();

        // k1 and k2 are principle curvatures of the shell
        // For 2d, k1 is the curvature, k2 = 0
        const Real k1_i = direction_corrector_ * k1_ave_[index_i];
        const Real k2_i = direction_corrector_ * k2_ave_[index_i];

        int counter = 0;
        while (distance_dummy < kernel_->CutOffRadius())
        {
            counter++;
            const Real factor_1 = 1 + counter * k1_i * particle_distance_;
            const Real factor_2 = 1 + counter * k2_i * particle_distance_;
            if (factor_1 <= 0 || factor_2 <= 0)
                break;
            const Real Vol_i_factor = factor_1 * factor_2;

            Real dW_ijV_j = kernel_->dW(distance_dummy, displacement_dummy) * Vol_i_factor;
            Vecd e_ij = displacement_dummy / distance_dummy;
            dW_ijV_j_ttl += dW_ijV_j;
            dW_ijV_j_e_ij_ttl += dW_ijV_j * e_ij;

            // calculate the position and volume of the next dummy particle
            pos_i_dummy += n_i * particle_distance_;
            displacement_dummy = pos_i_dummy - pos_j;
            distance_dummy = displacement_dummy.norm();
        }

        Vecd e_ij_corrected = dW_ijV_j_e_ij_ttl / dW_ijV_j_ttl;
        Real dW_ijV_j_corrected = dW_ijV_j_ttl * particle_distance_; // from surface area to volume

        // create new neighborhood
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, index_j, W_ij, dW_ijV_j_corrected, e_ij_corrected)
            : initializeNeighbor(neighborhood, distance, index_j, W_ij, dW_ijV_j_corrected, e_ij_corrected);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
ShellNeighborBuilderInnerWithContactKernel::ShellNeighborBuilderInnerWithContactKernel(SPHBody &body, SPHBody &contact_body) : NeighborBuilderInner(body)
{
    // create a reduced kernel with refined smoothing length for shell
    Real smoothing_length = contact_body.getSPHAdaptation().ReferenceSmoothingLength();
    kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(smoothing_length);
    kernel_->reduceOnce();
}
//=================================================================================================//
NeighborBuilderShellSelfContact::
    NeighborBuilderShellSelfContact(SPHBody &body)
    : BaseNeighborBuilderContactShell(body),
      k1_(body.getBaseParticles().registerStateVariableData<Real>("1stPrincipleCurvature")),
      k2_(body.getBaseParticles().registerStateVariableData<Real>("2ndPrincipleCurvature")),
      pos0_(body.getBaseParticles().registerStateVariableDataFrom<Vecd>("InitialPosition", "Position"))
{
    // create a unreduced kernel for shell self contact
    Real smoothing_length = body.getSPHAdaptation().ReferenceSmoothingLength();
    kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(smoothing_length);
}
//=================================================================================================//
void NeighborBuilderShellSelfContact::operator()(Neighborhood &neighborhood,
                                                 const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;

    const Vecd pos_j = list_data_j.second;
    const Vecd displacement = pos_i - pos_j;
    const Real distance = displacement.norm();
    const Real distance0 = (pos0_[index_i] - pos0_[index_j]).norm();

    // only particles within 1*dp distance should have force
    if (distance < kernel_->CutOffRadius() && distance0 > kernel_->CutOffRadius())
    {
        Real W_ijV_j_ttl = kernel_->W(distance, displacement);
        Real dW_ijV_j_ttl = kernel_->dW(distance, displacement);
        Vecd dW_ijV_j_e_ij_ttl = dW_ijV_j_ttl * displacement / (distance + TinyReal);

        // correct normal direction, make sure it points from index i to index j
        const Real direction_corrector = -SGN(displacement.dot(n_[index_j]));
        const Vecd n_j = n_[index_j] * direction_corrector; // normal direction of shell particle in cell linked list with index j

        Vecd pos_j_dummy = pos_j + n_j * particle_distance_;
        Vecd displacement_dummy = pos_i - pos_j_dummy;
        Real distance_dummy = displacement_dummy.norm();

        const Real k1_j = k1_[index_j] * direction_corrector; // mean curvature with corrected sign
        const Real k2_j = k2_[index_j] * direction_corrector; // mean curvature with corrected sign

        int counter = 0;
        // only add particles within 1*dp distance
        while (direction_corrector != 0 && distance_dummy < kernel_->CutOffRadius())
        {
            counter++;
            const Real factor_1 = 1 + counter * k1_j * particle_distance_;
            const Real factor_2 = 1 + counter * k2_j * particle_distance_;
            if (factor_1 <= 0 || factor_2 <= 0)
                break;
            const Real Vol_j_dummy = factor_1 * factor_2;
            Real dW_ijV_j = kernel_->dW(distance_dummy, displacement_dummy) * Vol_j_dummy;
            Vecd e_ij = displacement_dummy / distance_dummy;
            W_ijV_j_ttl += kernel_->W(distance_dummy, displacement_dummy) * Vol_j_dummy;
            dW_ijV_j_ttl += dW_ijV_j;
            dW_ijV_j_e_ij_ttl += dW_ijV_j * e_ij;

            // calculate the position and volume of the next dummy particle
            pos_j_dummy += n_j * particle_distance_;
            displacement_dummy = pos_i - pos_j_dummy;
            distance_dummy = displacement_dummy.norm();
        }

        Vecd e_ij_corrected = dW_ijV_j_e_ij_ttl / dW_ijV_j_ttl;
        Real W_ij_corrected = W_ijV_j_ttl * particle_distance_ / thickness_[index_j];
        // change surface mass and area to volumetric mass and volume
        Real dW_ijV_j_corrected = dW_ijV_j_ttl * particle_distance_;

        // create new neighborhood
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected)
            : initializeNeighbor(neighborhood, distance, index_j, W_ij_corrected, dW_ijV_j_corrected, e_ij_corrected);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
NeighborBuilderSurfaceContactFromShell::NeighborBuilderSurfaceContactFromShell(SPHBody &body, SPHBody &contact_body, bool normal_correction)
    : BaseNeighborBuilderContactFromShell(body, contact_body, normal_correction)
{
    Real source_smoothing_length = body.getSPHAdaptation().ReferenceSmoothingLength();
    Real target_smoothing_length = contact_body.getSPHAdaptation().ReferenceSmoothingLength();
    kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(0.5 * (source_smoothing_length + target_smoothing_length));
}
//=================================================================================================//
NeighborBuilderSurfaceContactFromSolid::NeighborBuilderSurfaceContactFromSolid(SPHBody &body, SPHBody &contact_body)
    : NeighborBuilderSurfaceContact(body, contact_body)
{
    Real dp_1 = body.getSPHBodyResolutionRef();
    Real dp_2 = contact_body.getSPHBodyResolutionRef();
    offset_W_ij_ = kernel_->W(0.5 * (dp_1 + dp_2), ZeroVecd);
}
//=================================================================================================//
void NeighborBuilderSurfaceContactFromSolid::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                            const Vecd &displacement, size_t index_j)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(std::max(kernel_->W(distance, displacement) - offset_W_ij_, Real(0)));
    neighborhood.dW_ij_.push_back(kernel_->dW(distance, displacement));
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(kernel_->e(distance, displacement));
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilderSurfaceContactFromSolid::initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                                                                const Vecd &displacement, size_t index_j)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = std::max(kernel_->W(distance, displacement) - offset_W_ij_, Real(0));
    neighborhood.dW_ij_[current_size] = kernel_->dW(distance, displacement);
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = kernel_->e(distance, displacement);
}
//=================================================================================================//
void NeighborBuilderSurfaceContactFromSolid::operator()(Neighborhood &neighborhood,
                                                        const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    if (distance < kernel_->CutOffRadius())
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j);
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
NeighborBuilderSplitInnerAdaptive::
    NeighborBuilderSplitInnerAdaptive(SPHBody &body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()),
      h_ratio_(body.getBaseParticles().getVariableDataByName<Real>("SmoothingLengthRatio")),
      level_(body.getBaseParticles().getVariableDataByName<int>("SmoothingLengthLevel")) {}
//=================================================================================================//
void NeighborBuilderSplitInnerAdaptive::
operator()(Neighborhood &neighborhood, const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;

    // Only add neighbors when dp_i >= dp_j
    // When dp_i == dp_j, only add neighbors when i < j
    size_t level_i = level_[index_i];
    if (size_t level_j = level_[index_j]; level_i > level_j || (level_i == level_j && index_i >= index_j))
        return;

    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    Real i_h_ratio = h_ratio_[index_i];
    Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[index_j]);
    Real cutoff_radius = kernel_->CutOffRadius(h_ratio_min);
    if (distance < cutoff_radius && index_i != index_j)
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, i_h_ratio, h_ratio_min)
            : initializeNeighbor(neighborhood, distance, displacement, index_j, i_h_ratio, h_ratio_min);
        neighborhood.current_size_++;
    }
};
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
