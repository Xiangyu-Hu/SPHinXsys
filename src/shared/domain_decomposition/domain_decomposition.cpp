#include "domain_decomposition.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>

namespace SPH
{
//=================================================================================================//
SlabDecomposition::SlabDecomposition(const BoundingBoxd &bounds, Real halo_width,
                                     int number_of_subdomains, int split_axis)
    : bounds_(bounds)
{
    if (number_of_subdomains < 1 || number_of_subdomains > execution::MaxSubdomains)
    {
        std::cout << "\n Error: SlabDecomposition supports 1 to " << execution::MaxSubdomains
                  << " subdomains, " << number_of_subdomains << " requested. \n";
        exit(1);
    }

    if (split_axis < 0)
    { // cut along the longest extent, which gives the smallest interface area
        const Vecd extent = bounds.upper_ - bounds.lower_;
        split_axis = 0;
        for (int axis = 1; axis < Dimensions; ++axis)
        {
            if (extent[axis] > extent[split_axis])
            {
                split_axis = axis;
            }
        }
    }

    subdomain_map_.number_of_subdomains_ = number_of_subdomains;
    subdomain_map_.split_axis_ = split_axis;
    subdomain_map_.halo_width_ = halo_width;

    const Real lower = bounds.lower_[split_axis];
    const Real upper = bounds.upper_[split_axis];
    const Real slab_thickness = (upper - lower) / Real(number_of_subdomains);
    for (int index = 0; index <= number_of_subdomains; ++index)
    {
        subdomain_map_.cut_plane_[index] = lower + Real(index) * slab_thickness;
    }
    // The outer planes must not clip the domain, positions exactly on the upper bound
    // included; subdomainOf() clamps outside positions to the end subdomains anyway.
    subdomain_map_.cut_plane_[0] = lower;
    subdomain_map_.cut_plane_[number_of_subdomains] = upper;

    if (number_of_subdomains > 1 && slab_thickness < MinimumSlabThickness())
    {
        std::cout << "\n Warning: slab thickness " << slab_thickness
                  << " is below twice the halo width " << subdomain_map_.halo_width_
                  << "; the halo would reach beyond the adjacent subdomain. Use fewer "
                  << "subdomains or a larger domain. \n";
    }
}
//=================================================================================================//
bool SlabDecomposition::rebalance(const StdVec<UnsignedInt> &particles_per_subdomain,
                                  Real relaxation)
{
    const int number_of_subdomains = subdomain_map_.number_of_subdomains_;
    if (number_of_subdomains < 2)
    {
        return false;
    }

    const UnsignedInt total_particles =
        std::accumulate(particles_per_subdomain.begin(), particles_per_subdomain.end(), UnsignedInt(0));
    if (total_particles == 0)
    {
        return false;
    }
    const Real target_share = Real(1) / Real(number_of_subdomains);

    // Cumulative particle fraction at each cut plane, then move the plane towards the
    // position where the fraction would equal the uniform target. Only interior planes
    // move; the outer ones are the domain bounds.
    StdVec<Real> cumulative_fraction(number_of_subdomains + 1, Real(0));
    for (int subdomain = 0; subdomain < number_of_subdomains; ++subdomain)
    {
        cumulative_fraction[subdomain + 1] =
            cumulative_fraction[subdomain] +
            Real(particles_per_subdomain[subdomain]) / Real(total_particles);
    }

    bool moved = false;
    const Real minimum_thickness = MinimumSlabThickness();
    StdVec<Real> new_planes(subdomain_map_.cut_plane_,
                            subdomain_map_.cut_plane_ + number_of_subdomains + 1);
    for (int plane = 1; plane < number_of_subdomains; ++plane)
    {
        const Real target_fraction = Real(plane) * target_share;
        const Real current_fraction = cumulative_fraction[plane];
        const Real local_share = cumulative_fraction[plane] - cumulative_fraction[plane - 1];
        if (local_share <= Real(0))
        {
            continue; // empty slab, no density estimate available
        }
        const Real slab_thickness =
            subdomain_map_.cut_plane_[plane] - subdomain_map_.cut_plane_[plane - 1];
        // Density along the split axis, approximated as uniform inside the slab.
        const Real correction = (target_fraction - current_fraction) / local_share * slab_thickness;
        new_planes[plane] += relaxation * correction;
    }

    // Enforce monotonicity and the minimum thickness before committing.
    for (int plane = 1; plane < number_of_subdomains; ++plane)
    {
        const Real lower_limit = new_planes[plane - 1] + minimum_thickness;
        const Real upper_limit = subdomain_map_.cut_plane_[number_of_subdomains] -
                                 Real(number_of_subdomains - plane) * minimum_thickness;
        new_planes[plane] = std::min(std::max(new_planes[plane], lower_limit), upper_limit);
        if (std::abs(new_planes[plane] - subdomain_map_.cut_plane_[plane]) >
            Real(1.0e-6) * minimum_thickness)
        {
            moved = true;
        }
        subdomain_map_.cut_plane_[plane] = new_planes[plane];
    }
    return moved;
}
//=================================================================================================//
std::string SlabDecomposition::describe() const
{
    std::ostringstream stream;
    stream << "SlabDecomposition: " << subdomain_map_.number_of_subdomains_
           << " subdomains along axis " << subdomain_map_.split_axis_
           << ", halo width " << subdomain_map_.halo_width_ << "\n";
    for (int subdomain = 0; subdomain < subdomain_map_.number_of_subdomains_; ++subdomain)
    {
        stream << "  [" << subdomain << "] ["
               << subdomain_map_.cut_plane_[subdomain] << ", "
               << subdomain_map_.cut_plane_[subdomain + 1] << ")\n";
    }
    return stream.str();
}
//=================================================================================================//
} // namespace SPH
