/**
 * @file    domain_partition.h
 * @brief   Cartesian decomposition of the system domain across MPI ranks.
 * @details This is pure geometry and topology: it knows which region of space each rank
 *          owns, which ranks are its neighbours, and which ranks need a copy of a given
 *          particle as a halo particle. It contains no MPI calls and no particle data,
 *          so it can be unit-tested standalone without launching a distributed job.
 *
 *          The halo width must be at least the largest kernel cutoff radius, since that
 *          is the distance over which particle interactions reach across a subdomain
 *          boundary. It is the same length scale that sizes the cell-linked-list cells.
 * @author  SPHinXsys multi-GPU draft
 */

#ifndef DOMAIN_PARTITION_H
#define DOMAIN_PARTITION_H

#include "data_type.h"

namespace SPH
{
/**
 * @class DomainPartition
 * @brief Splits a bounding box into a regular grid of subdomains, one per rank.
 *
 * Ranks are laid out in row-major order over the Cartesian grid, so rank r has grid
 * coordinates given by unflattening r with strides (1, n_x, n_x * n_y).
 */
class DomainPartition
{
  public:
    /**
     * @param system_bounds  Bounds of the whole simulation, as given to SPHSystem.
     * @param grid_dims      Number of subdomains along each axis; product must equal world_size.
     * @param halo_width     Halo thickness, >= largest kernel cutoff radius.
     * @param rank           This process's rank.
     * @param world_size     Total number of ranks.
     */
    DomainPartition(const BoundingBoxd &system_bounds, const Arrayi &grid_dims,
                    Real halo_width, int rank, int world_size);

    int Rank() const { return rank_; };
    int WorldSize() const { return world_size_; };
    Real HaloWidth() const { return halo_width_; };
    const Arrayi &GridDims() const { return grid_dims_; };

    /** Region this rank owns. Particles outside it must be migrated away. */
    const BoundingBoxd &SubdomainBounds() const { return subdomain_bounds_; };
    /** Owned region grown by the halo width; the region this rank needs data for. */
    const BoundingBoxd &HaloedBounds() const { return haloed_bounds_; };

    /** Subdomain of an arbitrary rank, without its halo. */
    BoundingBoxd subdomainBoundsOfRank(int rank) const;

    /** Rank owning a position, or -1 if the position is outside the system bounds. */
    int ownerRank(const Vecd &position) const;
    bool isOwnedByThisRank(const Vecd &position) const { return ownerRank(position) == rank_; };

    /**
     * Ranks adjacent to this one (up to 8 in 2D, 26 in 3D), excluding this rank and
     * excluding off-grid directions. Fixed for the lifetime of the partition, so the
     * communication topology is established once at setup.
     */
    const StdVec<int> &NeighbourRanks() const { return neighbour_ranks_; };

    /**
     * Neighbour ranks that need @p position as a halo particle, i.e. those whose haloed
     * region contains it. Appends to @p targets without clearing it.
     *
     * Cheap in practice: only particles near a subdomain boundary are ever tested, and
     * the loop is over the (at most 26) precomputed neighbours.
     */
    void haloTargetRanks(const Vecd &position, StdVec<int> &targets) const;

    /** True if the position lies within the halo shell of any neighbour, i.e. it must be sent. */
    bool isHaloCandidate(const Vecd &position) const;

  protected:
    BoundingBoxd system_bounds_;
    Arrayi grid_dims_;
    Real halo_width_;
    int rank_;
    int world_size_;

    Vecd subdomain_extent_; /**< size of one subdomain along each axis */
    BoundingBoxd subdomain_bounds_;
    BoundingBoxd haloed_bounds_;
    StdVec<int> neighbour_ranks_;
    StdVec<BoundingBoxd> neighbour_haloed_bounds_; /**< parallel to neighbour_ranks_ */

    Arrayi gridCoordsFromRank(int rank) const;
    int rankFromGridCoords(const Arrayi &coords) const;
    BoundingBoxd expandByHalo(const BoundingBoxd &bounds) const;
    void buildNeighbourTopology();
};
} // namespace SPH
#endif // DOMAIN_PARTITION_H
