#include "subdomain_runner.h"

#include <iostream>
#include <sstream>

namespace SPH
{
namespace execution
{
//=================================================================================================//
SubdomainWorkerPool::SubdomainWorkerPool(int number_of_workers)
    : number_of_workers_(number_of_workers)
{
    workers_.reserve(number_of_workers_);
    for (int subdomain_id = 0; subdomain_id < number_of_workers_; ++subdomain_id)
    {
        workers_.emplace_back([this, subdomain_id]
                              { workerLoop(subdomain_id); });
    }
}
//=================================================================================================//
SubdomainWorkerPool::~SubdomainWorkerPool()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        shutting_down_ = true;
        ++generation_;
    }
    start_condition_.notify_all();
    for (auto &worker : workers_)
    {
        if (worker.joinable())
        {
            worker.join();
        }
    }
}
//=================================================================================================//
void SubdomainWorkerPool::workerLoop(int subdomain_id)
{
    SubdomainScope scope(subdomain_id); // binding is permanent for this worker
    insideFanOutRef() = true;           // so that nested exec() calls do not fan out again
    std::size_t last_generation = 0;
    while (true)
    {
        const std::function<void(int)> *body = nullptr;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            start_condition_.wait(lock, [&]
                                  { return generation_ != last_generation; });
            last_generation = generation_;
            if (shutting_down_)
            {
                return;
            }
            body = body_;
        }

        std::exception_ptr caught;
        try
        {
            (*body)(subdomain_id);
        }
        catch (...)
        {
            caught = std::current_exception();
        }

        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (caught && !first_exception_)
            {
                first_exception_ = caught;
            }
            --pending_;
        }
        done_condition_.notify_one();
    }
}
//=================================================================================================//
void SubdomainWorkerPool::run(const std::function<void(int)> &body)
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        body_ = &body;
        pending_ = number_of_workers_;
        first_exception_ = nullptr;
        ++generation_;
    }
    start_condition_.notify_all();

    std::unique_lock<std::mutex> lock(mutex_);
    done_condition_.wait(lock, [&]
                         { return pending_ == 0; });
    body_ = nullptr;
    if (first_exception_)
    {
        std::exception_ptr to_rethrow = first_exception_;
        first_exception_ = nullptr;
        lock.unlock();
        std::rethrow_exception(to_rethrow);
    }
}
//=================================================================================================//
void SubdomainRunner::initialize(int number_of_subdomains, Mode mode)
{
    if (number_of_subdomains < 1 || number_of_subdomains > MaxSubdomains)
    {
        std::cout << "\n Error: SubdomainRunner supports 1 to " << MaxSubdomains
                  << " subdomains, " << number_of_subdomains << " requested. \n";
        exit(1);
    }
    if (numberOfSubdomains() != 1 && numberOfSubdomains() != number_of_subdomains)
    {
        std::cout << "\n Warning: the number of subdomains is being changed after it was "
                  << "already used; variables allocated so far keep their old number of "
                  << "replicas. Call initialize() before generating particles. \n";
    }

    number_of_subdomains_ = number_of_subdomains;
    mode_ = mode;
    numberOfSubdomainsRef() = number_of_subdomains;

    worker_pool_.reset();
    if (mode_ == Mode::Threaded && number_of_subdomains_ > 1)
    {
        worker_pool_ = std::make_unique<SubdomainWorkerPool>(number_of_subdomains_);
    }

    std::cout << describe() << std::endl;
}
//=================================================================================================//
std::string SubdomainRunner::describe() const
{
    std::ostringstream stream;
    stream << "SPHinXsys subdomain runner: " << number_of_subdomains_ << " subdomain(s), "
           << (mode_ == Mode::Sequential ? "sequential" : "threaded") << " mode";
    return stream.str();
}
//=================================================================================================//
} // namespace execution
} // namespace SPH
