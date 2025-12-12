---
layout: post
title:  "Contiguous Storage of Grid Data for Heterogeneous Computing (Part 1)"
date:   2025-12-11
categories: high performance computing 
---
Xiangyu Hu and Fan Gu

<p align="center"><img src="{{site.baseurl}}/assets/img/core-packages.png" alt="Heading Figure" style="max-width:100%; height:auto;">
<center>Fig. 1. Particle generation from an extruded prism.
  The gray surface indicates the clipped zero level set, the framed boxes indicate the core data package blocks used to describe the surface and the small spheres indicate the SPH particles generated and physically relaxed for later numerical simulation.</center> </p>

## Introduction

Many numerical methods in scientific computing, computer graphics,
and computational physics rely on structured Cartesian grids to discretize
scalar and vector fields. These grids are widely used in 
level set methods
[(Osher at al. 2001)](https://doi.org/10.1006/jcph.2000.6636),
computational fluid dynamics (CFD) solvers
[(Enright et al. 2002, ](https://doi.org/10.1006/jcph.2002.7166) [Han et al. 2014)]( https://doi.org/10.1016/j.jcp.2013.12.061),
and volumetric data processing pipelines [(Museth 2013)](https://doi.org/10.1145/2487228.2487235).
While Cartesian grids offer simplicity and regularity for finite difference and finite volume methods, large-scale simulations often involve highly sparse domains where only a small fraction of the grid is actively used at any time.

Naively storing and updating such sparse data results in excessive memory consumption 
and suboptimal memory access patterns, 
particularly on modern hardware architectures 
with hierarchical memory system and parallel processing constraints. 
To address this, data structures such as OpenVDB 
[(Museth 2013)](https://doi.org/10.1145/2487228.2487235). 
and SPGgrid [(Setaluri et al. 2014)](https://doi.org/10.1145/2661229.2661269). have been developed to reduce memory overhead 
and improve computational performance. 
These frameworks leverage hierarchical spatial partitioning or page-based layouts to exploit sparsity. 
However, they are primarily optimized for CPU-based execution, 
and often suffer from performance degradation when ported to Graphic Processing Unit (GPU) architectures. 
Key challenges include high random-access latency and limited concurrency during data update.

GPU has become prevalent in accelerating a wide range of numerical simulations 
due to their high parallel processing capability. 
Despite the high computing power provided, 
leveraging GPU typically requires the adoption of specialized programming models, 
such as CUDA(NVIDIA), HIP(AMD), SYCL, and others, 
which presents a steep learning curve for end-users.

In our prior implementations in the open-source SPH (Smoothed Particle Hydrodynamics) 
multi-physics library SPHinXsys
[(Han et al. 2014,]( https://doi.org/10.1016/j.jcp.2013.12.061)
[(Zhnag et al. 2021)]( https://doi.org/10.1016/j.cpc.2021.108066), 
memory allocation for activated grid regions was managed via dynamic memory pools. 
While this approach yields efficient memory usage, 
the resulting data structures are not inherently compatible with GPU execution. 
In this work, we present a redesigned storage architecture aimed at improving computational efficiency and enabling GPU compatibility. 
The new design is optimized for unified usage across both host and device platforms, 
abstracting underlying implementation details from the user.

This work introduces an optimized data structure tailored for GPU execution, 
along with a newly developed computing kernels. 
The key contributions are:

- Abstraction of GPU-specific and SYCL-specific implementation details, thereby simplifying usage for end-users.

- Unified codebase enabling seamless execution on both CPU and GPU platforms through a single code logic.

- Heterogeneous computing is achieved for both sparse-grid and grid-particle coupling numerical simulations.

The proposed method incorporates computing kernels, data structures and execution strategies. 
The computing kernel orchestrates the computational process, 
while data structure manage storage and synchronization of data across CPU and GPU. 
Upon execution, the computing kernel automatically retrieves the appropriate host 
or device variable based on the selected execution strategy, 
thus ensuring consistency and portability.

## Design Objectives and Main Components

The present objective of using sparse-grid storage in SPHinXsys is two folded.
One is to represent the body surfaces with a level set field, 
similar to the application of visualization,
which is used to generate SPH particles for numerical simulation. 
The other is level-set based operations, 
which are both computation and memory intensive, 
including small-feature cleaning, 
water-tight-surface ensuring and 
convolution integrals with SPH smoothing kernel
[(Yu et al. 2023)](https://doi.org/10.1016/j.cpc.2023.108744).
Note that these operations are essential for dynamically reorganizing 
the SPH particles to be body-fitted and well distributed 
for stable and accurate numerical simulations.

Therefore, the present method is composed of three primary components as 
sparse-grid (narrow-band) storage, access and computation, respectively. 
While the first component is managed by the `MeshWithGridDataPackage` class,
the second provides the functions for computing stencils 
of interpolations and differentials based on the access function `NeighbourIndexShift`.
Two patterns of computation, 
i.e. `ALLMeshDynamics` and `MeshPackageDynamics`,
are provided to transverse the entire Cartesian grid 
and activated cells, respectively, 
with execution policies (e.g., sequential, 
or parallel executions on host and device) to be applied.

In the next post, 
I will give the details on the implementation of these three components.

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
