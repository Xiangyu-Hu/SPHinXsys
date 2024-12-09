---
layout: post
title:  "Multi-layer data structure in SPHinXsys"
date:   2024-11-27
categories: sphinxsys update
---

Designed for heterogeneous computing architecture,
three layers, i.e. outer, kernel-shell and kernel,
dependent on their distances to computing kernel, are used in SPHinXsys.
Accordingly, the characteristics of data structure are different at each layer.

The outer layer is the problem-definition layer, at which the physical problem for the simulation or computing is defined.
The global variables are defined, the preprocess is done by generating particles respect to geometric information.
Third-parties environment and methods, e.g. Simbody environment and methods, are defined.
The data structure from standard C++ library and other third-parties libraries can be directly used for this propose.

The kernel-shell layer is the SPH method-definition layer,
which is used to define individual SPH method,
including sub-methods, e.g. those for particle adaptation and material properties.
At this layer, the interface data structures, namely, discrete and singular variables are
used to transfer the outer-layer data to those used in computing kernel.
These interface are also containers of the transformed data,
which are trivially-copyable and can be used on a computing device, typically GPU.

The kernel layer defines all the computing via computing kernels running on the device,
and only instantiates at the beginning of computing
with the help of the kernel-shell layer interface and execution implementation.
The data transformed from kernel-shell layer has the form of pointer.
For that from discrete variable, it can be accessed just as arrays.

Still, one can use other constant trivially-copyable data structures in any of the 3 layers.
Beside simple types, one typical more complex example is the RiemannSolver,
which is a trivially-copyable constant and can be use in computing kernels directly.

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
