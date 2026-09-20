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
 * @file 	subdomain_runner.h
 * @brief 	Execution of a fan-out body over the subdomains of a host side run.
 * @details Two modes, and the choice between them is the whole point of the CPU
 *          decomposition path:
 *
 *          - Sequential: the subdomains are visited one after another on the calling
 *            thread. Fully deterministic and single threaded, so a debugger steps
 *            through the decomposition logic, a data race cannot mask a logic error,
 *            and a failure is bit reproducible. This is the default.
 *          - Threaded: one persistent host thread per subdomain, with the same barrier
 *            structure as the multi-GPU path. This is what reproduces the ordering
 *            bugs of the real thing, and is meant to be run under ThreadSanitizer.
 *
 *          Running both and comparing is the intended workflow: a difference between
 *          sequential and threaded is a synchronization bug, a difference between one
 *          subdomain and several is a decomposition bug. Separating those two questions
 *          is far cheaper here than on a multi-GPU node.
 * @author	Niki Loppi
 */

#ifndef SUBDOMAIN_RUNNER_H
#define SUBDOMAIN_RUNNER_H

#include "subdomain_scope.h"

#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace SPH
{
namespace execution
{
/**
 * @class SubdomainWorkerPool
 * @brief One persistent host thread per subdomain, all executing the same callable.
 * @details Persistent because creating threads per call would dominate the runtime:
 *          a single acoustic step issues dozens of fan-outs. Shared by the host path
 *          and, on the SYCL path, by DeviceEnvironment.
 */
class SubdomainWorkerPool
{
  public:
    explicit SubdomainWorkerPool(int number_of_workers);
    ~SubdomainWorkerPool();

    /** Run body(subdomain_id) on every worker and block until all have returned.
     *  The first exception escaping any worker is rethrown on the caller. */
    void run(const std::function<void(int)> &body);

  private:
    void workerLoop(int subdomain_id);

    int number_of_workers_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable start_condition_;
    std::condition_variable done_condition_;
    const std::function<void(int)> *body_ = nullptr;
    std::size_t generation_ = 0;
    int pending_ = 0;
    bool shutting_down_ = false;
    std::exception_ptr first_exception_;
};

/**
 * @class SubdomainRunner
 * @brief Host side driver of a decomposed run.
 */
class SubdomainRunner
{
  public:
    enum class Mode
    {
        Sequential, /**< deterministic, single threaded; the default */
        Threaded    /**< one thread per subdomain, mirroring the multi-GPU path */
    };

    SubdomainRunner(const SubdomainRunner &) = delete;
    void operator=(const SubdomainRunner &) = delete;

    static SubdomainRunner &getInstance()
    {
        static SubdomainRunner instance;
        return instance;
    }

    /** Must be called before any particle data is generated, since it fixes how many
     *  replicas each variable allocates. */
    void initialize(int number_of_subdomains, Mode mode = Mode::Sequential);

    int NumberOfSubdomains() const { return number_of_subdomains_; };
    Mode getMode() const { return mode_; };
    std::string describe() const;

    /** Run body(subdomain_id) once per subdomain, with the subdomain bound. */
    template <class Body>
    void forEachSubdomain(const Body &body)
    {
        if (number_of_subdomains_ == 1 || insideFanOut())
        { // single subdomain, or a nested fan-out already bound to one
            body(currentSubdomainID());
            return;
        }

        if (mode_ == Mode::Sequential)
        {
            for (int subdomain_id = 0; subdomain_id < number_of_subdomains_; ++subdomain_id)
            {
                SubdomainScope scope(subdomain_id);
                insideFanOutRef() = true;
                try
                {
                    body(subdomain_id);
                }
                catch (...)
                {
                    insideFanOutRef() = false;
                    throw;
                }
                insideFanOutRef() = false;
            }
            return;
        }

        std::function<void(int)> wrapped = [&](int subdomain_id)
        { body(subdomain_id); };
        worker_pool_->run(wrapped); // the workers are permanently marked as in-fan-out
    };

  private:
    SubdomainRunner() = default;
    ~SubdomainRunner() = default;

    int number_of_subdomains_ = 1;
    Mode mode_ = Mode::Sequential;
    std::unique_ptr<SubdomainWorkerPool> worker_pool_;
};

inline SubdomainRunner &subdomain_runner = SubdomainRunner::getInstance();
} // namespace execution
} // namespace SPH
#endif // SUBDOMAIN_RUNNER_H
