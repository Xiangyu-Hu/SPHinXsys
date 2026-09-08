/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	domain_decomposition.h
 * @brief 	Geometric decomposition of the computational domain over the devices.
 * @details A slab (one dimensional) decomposition is used: the domain is cut by
 *          planes normal to one axis, one slab per subdomain, so that a subdomain has
 *          at most two neighbors. That keeps the exchange pattern, the halo bookkeeping
 *          and the load balancing simple, which is what matters for a first
 *          implementation; a recursive or Cartesian decomposition can replace it
 *          later behind the same SubdomainMap interface.
 * @author	Niki Loppi
 */

#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "base_data_package.h"
#include "subdomain_scope.h"
#include "sphinxsys_containers.h"

namespace SPH
{
/**
 * @class SubdomainMap
 * @brief Device-copyable description of which subdomain a position belongs to.
 * @details Trivially copyable and small, so it is captured by value into kernels
 *          rather than being reached through a USM pointer.
 */
class SubdomainMap
{
  public:
    SubdomainMap() = default;

    int number_of_subdomains_ = 1;
    int split_axis_ = 0;
    Real halo_width_ = Real(0);
    /** Cut planes along split_axis_, ascending. cut_plane_[0] and
     *  cut_plane_[number_of_subdomains_] are the outer bounds of the domain. */
    Real cut_plane_[execution::MaxSubdomains + 1] = {};

    /** Owner of a position. Positions outside the domain are clamped to the end
     *  subdomains, so a particle leaving the domain is never left unowned. */
    int subdomainOf(const Vecd &position) const
    {
        const Real coordinate = position[split_axis_];
        for (int subdomain = 0; subdomain < number_of_subdomains_ - 1; ++subdomain)
        {
            if (coordinate < cut_plane_[subdomain + 1])
            {
                return subdomain;
            }
        }
        return number_of_subdomains_ - 1;
    };

    /** Whether a position lies within the halo band that `subdomain` needs, that is
     *  its own slab widened by the interaction cut-off on both sides. */
    bool inHaloBandOf(const Vecd &position, int subdomain) const
    {
        const Real coordinate = position[split_axis_];
        return coordinate >= cut_plane_[subdomain] - halo_width_ &&
               coordinate < cut_plane_[subdomain + 1] + halo_width_;
    };

    /** Neighbor on the given side, or -1 at the ends of the domain.
     *  side 0 is the lower (left) neighbor, side 1 the upper (right) one. */
    int neighborOf(int subdomain, int side) const
    {
        const int neighbor = side == 0 ? subdomain - 1 : subdomain + 1;
        return neighbor >= 0 && neighbor < number_of_subdomains_ ? neighbor : -1;
    };

    static int oppositeSide(int side) { return 1 - side; };
};

/**
 * @class SlabDecomposition
 * @brief Host side owner of the cut planes, including their dynamic adjustment.
 */
class SlabDecomposition
{
  public:
    /**
     * @param bounds        domain bounds, normally the SPHSystem bounds
     * @param halo_width    interaction cut-off; the halo band has this width
     * @param number_of_subdomains  the number of devices, or of host subdomains
     * @param split_axis    axis to cut along; by default the longest one
     */
    SlabDecomposition(const BoundingBoxd &bounds, Real halo_width,
                      int number_of_subdomains, int split_axis = -1);

    const SubdomainMap &getSubdomainMap() const { return subdomain_map_; };
    int NumberOfSubdomains() const { return subdomain_map_.number_of_subdomains_; };
    int SplitAxis() const { return subdomain_map_.split_axis_; };
    Real HaloWidth() const { return subdomain_map_.halo_width_; };
    Real CutPlane(int index) const { return subdomain_map_.cut_plane_[index]; };

    /**
     * Move the cut planes so that each subdomain holds a comparable share of the work.
     * @param particles_per_subdomain measured load, currently the owned particle count.
     * @param relaxation fraction of the correction applied per call, damping the
     *        oscillation that a fully applied correction would cause.
     * @return true when any cut plane moved by more than a tolerance, in which case
     *         the caller must trigger a migration before the next step.
     *
     * The estimate assumes the particle density is locally uniform along the split
     * axis, which is what makes a one dimensional rebalance cheap: no particle
     * histogram is needed, only the per-subdomain counts.
     */
    bool rebalance(const StdVec<UnsignedInt> &particles_per_subdomain,
                   Real relaxation = Real(0.5));

    /** Minimum slab thickness, enforced so that a subdomain never becomes thinner
     *  than its own halo, which would make a particle appear in the halo of a
     *  non-adjacent subdomain and break the two-neighbor assumption. */
    Real MinimumSlabThickness() const { return Real(2) * subdomain_map_.halo_width_; };

    std::string describe() const;

  protected:
    BoundingBoxd bounds_;
    SubdomainMap subdomain_map_;
};
} // namespace SPH
#endif // DOMAIN_DECOMPOSITION_H
