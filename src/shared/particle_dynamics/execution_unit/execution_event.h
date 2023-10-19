#ifndef SPHINXSYS_EXECUTIONEVENT_HPP
#define SPHINXSYS_EXECUTIONEVENT_HPP

#include <sycl/sycl.hpp>

namespace SPH::execution
{
class ExecutionEvent
{
  public:
    ExecutionEvent() = default;
    ExecutionEvent(const sycl::event &event);
    ExecutionEvent(const std::vector<sycl::event> &event);

    ExecutionEvent(const ExecutionEvent &executionEvent) = default;
    ExecutionEvent(ExecutionEvent &&executionEvent) = default;

    ExecutionEvent &operator=(const ExecutionEvent &event) = default;
    ExecutionEvent &operator=(ExecutionEvent &&event) = default;

    ExecutionEvent &operator=(sycl::event event);

    const std::vector<sycl::event> &getEventList() const;

    ExecutionEvent &add(sycl::event event);
    ExecutionEvent &add(const ExecutionEvent &event);

    void wait();
    ExecutionEvent &then(std::function<void()> &&func);

  private:
    std::vector<sycl::event> event_list_;
};
} // namespace SPH::execution

#endif // SPHINXSYS_EXECUTIONEVENT_HPP
