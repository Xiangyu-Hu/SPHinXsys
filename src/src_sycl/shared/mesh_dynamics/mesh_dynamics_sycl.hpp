namespace SPH
{
template <typename FunctionOnData>
void BaseMeshDynamics::package_parallel_for(const execution::ParallelDevicePolicy &par_device, const FunctionOnData &function)
{
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue.submit([&](sycl::handler &cgh)
                      { cgh.parallel_for(execution_instance.getUniformNdRange(num_grid_pkgs_ - 2), [=](sycl::nd_item<1> index)
                                        {
                                if(index.get_global_id(0) + 2< num_grid_pkgs_)
                                    function(index.get_global_id(0) + 2); }); })
        .wait_and_throw();
}
}