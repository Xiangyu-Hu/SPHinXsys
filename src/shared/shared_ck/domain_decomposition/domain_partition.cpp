#include "domain_partition.h"

#include <cmath>
#include <iostream>

namespace SPH
{
//=================================================================================================//
DomainPartition::DomainPartition(const BoundingBoxd &system_bounds, const Arrayi &grid_dims,
                                 Real halo_width, int rank, int world_size)
    : system_bounds_(system_bounds), grid_dims_(grid_dims), halo_width_(halo_width),
      rank_(rank), world_size_(world_size)
{
    int total_subdomains = 1;
    for (int i = 0; i != Dimensions; ++i)
    {
        if (grid_dims_[i] < 1)
        {
            std::cout << "\n ERROR: DomainPartition grid dimension " << i
                      << " is " << grid_dims_[i] << ", must be at least 1!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        total_subdomains *= grid_dims_[i];
    }

    if (total_subdomains != world_size_)
    {
        std::cout << "\n ERROR: DomainPartition grid holds " << total_subdomains
                  << " subdomains but there are " << world_size_ << " ranks!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    if (rank_ < 0 || rank_ >= world_size_)
    {
        std::cout << "\n ERROR: DomainPartition rank " << rank_
                  << " is outside [0, " << world_size_ << ")!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    for (int i = 0; i != Dimensions; ++i)
    {
        subdomain_extent_[i] =
            (system_bounds_.upper_[i] - system_bounds_.lower_[i]) / Real(grid_dims_[i]);

        // A subdomain thinner than the halo would need data from ranks beyond its immediate
        // neighbours, which the 26-neighbour topology below cannot express.
        if (subdomain_extent_[i] < halo_width_)
        {
            std::cout << "\n ERROR: DomainPartition subdomain extent " << subdomain_extent_[i]
                      << " along axis " << i << " is thinner than the halo width "
                      << halo_width_ << "!" << std::endl;
            std::cout << "\n Use fewer ranks along this axis, or a coarser resolution." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }

    subdomain_bounds_ = subdomainBoundsOfRank(rank_);
    haloed_bounds_ = expandByHalo(subdomain_bounds_);
    buildNeighbourTopology();
}
//=================================================================================================//
Arrayi DomainPartition::gridCoordsFromRank(int rank) const
{
    Arrayi coords = Arrayi::Zero();
    int remainder = rank;
    for (int i = 0; i != Dimensions; ++i)
    {
        coords[i] = remainder % grid_dims_[i];
        remainder /= grid_dims_[i];
    }
    return coords;
}
//=================================================================================================//
int DomainPartition::rankFromGridCoords(const Arrayi &coords) const
{
    int rank = 0;
    int stride = 1;
    for (int i = 0; i != Dimensions; ++i)
    {
        if (coords[i] < 0 || coords[i] >= grid_dims_[i])
        {
            return -1; // off the grid: no such neighbour
        }
        rank += coords[i] * stride;
        stride *= grid_dims_[i];
    }
    return rank;
}
//=================================================================================================//
BoundingBoxd DomainPartition::subdomainBoundsOfRank(int rank) const
{
    Arrayi coords = gridCoordsFromRank(rank);
    Vecd lower = Vecd::Zero();
    Vecd upper = Vecd::Zero();
    for (int i = 0; i != Dimensions; ++i)
    {
        lower[i] = system_bounds_.lower_[i] + Real(coords[i]) * subdomain_extent_[i];
        upper[i] = lower[i] + subdomain_extent_[i];
    }
    return BoundingBoxd(lower, upper);
}
//=================================================================================================//
BoundingBoxd DomainPartition::expandByHalo(const BoundingBoxd &bounds) const
{
    Vecd halo = Vecd::Constant(halo_width_);
    return BoundingBoxd(bounds.lower_ - halo, bounds.upper_ + halo);
}
//=================================================================================================//
void DomainPartition::buildNeighbourTopology()
{
    neighbour_ranks_.clear();
    neighbour_haloed_bounds_.clear();

    Arrayi this_coords = gridCoordsFromRank(rank_);

    int number_of_offsets = 1;
    for (int i = 0; i != Dimensions; ++i)
    {
        number_of_offsets *= 3; // each axis contributes -1, 0, +1
    }

    for (int k = 0; k != number_of_offsets; ++k)
    {
        Arrayi offset = Arrayi::Zero();
        int remainder = k;
        bool is_self = true;
        for (int i = 0; i != Dimensions; ++i)
        {
            offset[i] = remainder % 3 - 1;
            remainder /= 3;
            if (offset[i] != 0)
            {
                is_self = false;
            }
        }
        if (is_self)
        {
            continue;
        }

        int neighbour_rank = rankFromGridCoords(this_coords + offset);
        if (neighbour_rank < 0)
        {
            continue; // subdomain boundary coincides with the system boundary
        }

        neighbour_ranks_.push_back(neighbour_rank);
        neighbour_haloed_bounds_.push_back(expandByHalo(subdomainBoundsOfRank(neighbour_rank)));
    }
}
//=================================================================================================//
int DomainPartition::ownerRank(const Vecd &position) const
{
    if (!system_bounds_.checkContain(position))
    {
        return -1;
    }

    Arrayi coords = Arrayi::Zero();
    for (int i = 0; i != Dimensions; ++i)
    {
        coords[i] = int(std::floor((position[i] - system_bounds_.lower_[i]) / subdomain_extent_[i]));
        // A position exactly on the upper system bound floors to grid_dims_[i];
        // it belongs to the last subdomain along that axis.
        coords[i] = SMAX(0, SMIN(coords[i], grid_dims_[i] - 1));
    }
    return rankFromGridCoords(coords);
}
//=================================================================================================//
void DomainPartition::haloTargetRanks(const Vecd &position, StdVec<int> &targets) const
{
    for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
    {
        if (neighbour_haloed_bounds_[k].checkContain(position))
        {
            targets.push_back(neighbour_ranks_[k]);
        }
    }
}
//=================================================================================================//
bool DomainPartition::isHaloCandidate(const Vecd &position) const
{
    for (size_t k = 0; k != neighbour_haloed_bounds_.size(); ++k)
    {
        if (neighbour_haloed_bounds_[k].checkContain(position))
        {
            return true;
        }
    }
    return false;
}
//=================================================================================================//
} // namespace SPH
