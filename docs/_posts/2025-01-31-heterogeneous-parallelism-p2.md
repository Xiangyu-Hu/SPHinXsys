---
layout: post
title:  "Towards heterogeneous parallelism for SPHinXsys (Part 2)"
date:   2025-01-31
categories: high performance computing 
---
Xiangyu Hu and Alberto Guarnieri

## Implementing SYCL kernels in SPHinXsys

The long term aim for SPHinXsys development is
employing unified codebase for the SPH methods
so the multi-physics simulation can be carried on
heterogeneous computing systems in a collaborative way.
Therefore each SPH algorithm can be executed on host or device
in sequential or parallel.

Before implementing SYCL kernels,
SPHinXsys already has defined two execution policies,
namely `sequenced_policy` and `parallel_policy`
for CPU computing.
We deliberately implement the sequential execution on CPU
because it is essential for easy debug propose.
As SPHinXsys splits the execution,
namely `particle_for` and `particle_reduce`,
and SPH methods, i.e. the classes inherited from the base
`local_dynamics`
with the same source code for both executions.

Our basic concept on implementing SYCL kernels is to extend
the same programming pattern for CPU computing
by adding an extra execution policy,
namely `parallel_device_policy`  
to identify that the same source code of SPH methods
will be executed in parallel on the device.

Though the concept is straightforward,
its implementation is not trivial since it relies
on the an important assumption that the codes written in
computing kernels should work for all execution policies.
However, the source code of SPH methods previously used for CPU computing
generally are not executable on device
due to the widely used references, virtual functions
and standard library dynamic memory allocations.
Therefore, an extended program pattern
based on SYCL unified share memory (USM) is developed.
The extended pattern has three layer,
i.e. outer, kernel-shell and kernel,
dependent on their distances to computing kernel
which is the only one executed by the execution algorithms.

The outer layer is the same as before and defines the physical problem.
Within this layer, the global variables are defined,
the preprocess is done by generating particles respect to geometric information.
Third-parties environment and methods,
e.g. Simbody environment and methods, are referred directly.
The kernel-shell layer is the SPH method definition layer,
which can be obtained from the CPU computing code with minimum modification.
This layer is used to define individual SPH method.
At this layer, the interface data structures,
namely `discrete_` and `singular_variables`
are used to transfer the outer-layer data to those used in computing kernel.
The kernel layer defines all the computing via kernels,
and is obtained by splitting out the computing function from
original method classes.
The objects of this layer are only instantiated and dispatched
at the beginning of computing with
the help of the kernel-shell layer interface and execution implementation.
Note that the data transformed from kernel-shell layer has the form of raw pointer,
and can be accessed on host and device just as arrays.

## Cell linked-list, direct search and sorting

Since the limitation of device hardware,
SYCL specification does not allow the usage of
`concurrent_vector`,
which is has been used in SPHinXys for constructing the cell linked-list for
particle neighbor searching algorithms.
To achieve high performance on constructing cell linked-list
on device concurrently,
atomic operations are used to avoid the possible thread conflicts.

Again, due to the memory limit of GPU,
the full particle configuration,
including neighbor list, kernel values
used in SPhinXsys can not be used anymore.
Therefore, direct search is used,
i.e. the kernel values will be computed for all particles
within the nearest cells in the computing kernels directly.

Similarly,
due to poor performance of GPU on recursive algorithms,
the particle sorting executed by device uses
`radix_sort` other than `quick_sort`.
Note that, this is the only computing kernel that
used different codes for device execution.

In the next post, I will present the performance evaluation
for SYCL kernels implemented in SPHinXsys.

<script src="https://giscus.app/client.js"
        data-repo="Xiangyu-Hu/SPHinXsys"
        data-repo-id="MDEwOlJlcG9zaXRvcnkxODkwNzAxNDA="
        data-category="Announcements"
        data-category-id="DIC_kwDOC0T7PM4CPNAR"
        data-mapping="pathname"
        data-strict="0"
        data-reactions-enabled="1"
        data-emit-metadata="0"
        data-input-position="bottom"
        data-theme="light"
        data-lang="en"
        crossorigin="anonymous"
        async>
</script>
