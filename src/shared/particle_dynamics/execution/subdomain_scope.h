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
 * @file 	subdomain_scope.h
 * @brief 	Thread-local subdomain selection for domain decomposed execution.
 * @details A domain decomposed run splits the particles into subdomains, each with its
 *          own replica of every particle variable, of every computing kernel and of the
 *          particle counters. Which replica a piece of code sees is not passed as an
 *          argument; it is bound to the calling host thread by a SubdomainScope and read
 *          back through currentSubdomainID(). That indirection is what lets the whole
 *          of the CK dynamics stay unchanged between a single-subdomain run, a multi-GPU
 *          run and a multi-subdomain CPU run.
 *
 *          A subdomain is a unit of decomposition, not a unit of hardware. Under
 *          ParallelMultiDevicePolicy one subdomain maps to one GPU; under
 *          ParallelMultiHostPolicy it maps to a host thread, or to nothing at all when
 *          the runner is in sequential mode. Keeping one name for the concept is
 *          deliberate: a bug reproduced on the CPU path is then literally the same bug.
 *
 *          This header carries no backend dependency and is included from the shared
 *          part of the library. In a non-decomposed run currentSubdomainID() is always
 *          0 and every construct here degenerates to the previous behavior.
 * @author	Niki Loppi
 */

#ifndef SUBDOMAIN_SCOPE_H
#define SUBDOMAIN_SCOPE_H

namespace SPH
{
namespace execution
{
/** Compile-time upper bound on the number of subdomains in one process. It sizes the
 *  small per-subdomain pointer arrays held by each variable, so it is kept small. */
constexpr int MaxSubdomains = 8;

/** Subdomain driven by the calling host thread. */
inline int &currentSubdomainIDRef()
{
    static thread_local int current_subdomain_id = 0;
    return current_subdomain_id;
}

inline int currentSubdomainID() { return currentSubdomainIDRef(); }

/** Number of subdomains in the run. Set once during initialization, by
 *  DeviceEnvironment on the SYCL path or by SubdomainRunner on the host path. */
inline int &numberOfSubdomainsRef()
{
    static int number_of_subdomains = 1;
    return number_of_subdomains;
}

inline int numberOfSubdomains() { return numberOfSubdomainsRef(); }

inline bool isMultiSubdomain() { return numberOfSubdomains() > 1; }

/** Set while the calling thread executes inside a fan-out. Dynamics compose: an
 *  interaction algorithm's exec() fans out and then calls exec() of its pre- and
 *  post-processes, which would fan out again. A nested fan-out must degenerate to a
 *  plain call on the already bound subdomain, both to avoid re-entering the worker
 *  pool and because the work is already distributed. */
inline bool &insideFanOutRef()
{
    static thread_local bool inside_fan_out = false;
    return inside_fan_out;
}

inline bool insideFanOut() { return insideFanOutRef(); }

/**
 * @class SubdomainScope
 * @brief RAII binding of the calling host thread to one subdomain.
 * @details Nested scopes are allowed; the previous binding is restored on exit, which
 *          lets code temporarily reach into another subdomain's replica, as the halo
 *          exchange does when it reads a neighbor's send buffer.
 */
class SubdomainScope
{
  public:
    explicit SubdomainScope(int subdomain_id)
        : previous_subdomain_id_(currentSubdomainIDRef())
    {
        currentSubdomainIDRef() = subdomain_id;
    };
    ~SubdomainScope() { currentSubdomainIDRef() = previous_subdomain_id_; };

    SubdomainScope(const SubdomainScope &) = delete;
    SubdomainScope &operator=(const SubdomainScope &) = delete;

  private:
    int previous_subdomain_id_;
};
} // namespace execution
} // namespace SPH
#endif // SUBDOMAIN_SCOPE_H
