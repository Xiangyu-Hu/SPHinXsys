---
layout: post
title:  "Contiguous Storage of Grid Data for Heterogeneous Computing (Part 4)"
date:   2026-01-13
categories: high performance computing 
---
Xiangyu Hu and Fan Gu

<p align="center"><img src="{{site.baseurl}}/assets/img/multi-resolution-particles.png" 
alt="Heading Figure" style="max-width:100%; height:auto;">
<center> Fig. 4. Multi-resolution particles generated using multi-resolution level-set field. </center> </p>

## Performance Evaluation

Compared to the previous implementation on SPHinXsys, memory usage has been significantly reduced, in the present work, by eliminating the need to store a data address for each individual data entry within an activated cell. As a result, memory usage is reduced by 6×6×8=288 bytes per cell per mesh variable in each layer. Given the limited memory resources on GPU devices, this optimization is particularly valuable.

| Test Setting                 | OpenVDB | SPGrid | SPHinXsys |
|------------------------------|---------|--------|-----------|
|Sequential, 1 thread          |79.563   |77.2598 | 22.948    |
|Sequential, 4 threads         |32.322   |29.8752 |7.429      |
|Stencil, 1 thread             |1013.162 |229.572 |59.972     |
|Stencil, 4 threads            |303.902  |68.8437 |21.378     |

<p> <center>Tab. 1. Timing for operations using SPGrid, OpenVDB, and SPHinXsys.</center> </p>

To evaluate the computational performance with previous popular sparse-grid methods, namely OpenVDB and SPGrid, a shelled sphere is used as the shared test geometry across all methods. The shell is centered at (0.5, 0.5, 0.5) and has an inner radius of 0.3 and outer one of 0.31, with a resolution of 1/1024. A sequential access and a stencil operation are carried out on all the activated data.

Note that, for sequential access, the `LeafManager` is employed for OpenVDB. It is initialized with the target grid, after which each leaf node is iterated. Within each leaf, the offset filter provided by OpenVDB is used to iterate through active data points. The sequential benchmark involves making minor changes to each value stored in an activated cell.

For stencil operations, the seven-point Laplacian operator provided by OpenVDB is used. A comparable stencil computation is implemented by using `NeighbourIndexShift` for SPHinXsys. In both cases, all active data are visited, and a Laplacian operation is applied at each active data too. Note that, due to the lack of reference on GPU execution, here, we only provide the CPU performance of the present method. Also note that, for this test OpenVDB requires only one-third of the data storage compared to SPHinXsys. This difference is expected, as OpenVDB initializes level sets on a data-wise basis, allowing each data to be individually categorized. On the other hand, when data access and stecil operations are tested, relative more data will be accessed or modified in the present method.

The evaluation results for both single and multi-thread execution are shown in Tab. 1, which shows, in both scenarios, the presented method outperforms OpenVDB and SPGrid. This improvement is attributed to the overhead of tree traversal required during neighbor lookup in their methods. Although OpenVDB implements a software cache to mitigate this issue, the traversal process remains a significant bottleneck, especially when frequent neighbor access is required for operations like the Laplacian. In addition, the prefetch enabled by sequential access of contiguous memory storage also significantly influence the results.

## Numerical Examples

Here, we demonstrate two typical applications, both as pre-processing tool for SPH simulation, of the present method.

In the first case, SPH particles are generated from a gear model in STL file for a gearbox simulation, as shown in Fig. 5. While the STL model visually is perfectly water-tight, there are many triangle faces are not strictly connected and a straightforward application of the signed-distance function from a triangle mesh will leads to many particles generated (leaked) away from the surface due to the contain condition is wrongly predicted.

<p align="center"><img src="https://arxiv.org/html/2512.11473v1/gear.jpeg" 
alt="Heading Figure" sstyle="max-width:100%; height:auto;">
<center> Fig. 5. STL triangle mesh of a gear. The leaking regions are indicated with yellow lines. </center> </p>

With the present method, after the initial level set is generated, the contain consistency correction algorithm is applied and run on GPU with high computational efficiency. The corrected level set field and final generated SPH particles is shown in Fig. 6.

<p align="center"><img src="https://arxiv.org/html/2512.11473v1/x1.png" 
alt="Heading Figure" sstyle="max-width:100%; height:auto;">
<center> Fig. 6. Corrected level-set field, zero-level contour and the finally generated SPH particles. The zero-level contour is colored with green. </center> </p>

Note that gear surface is represented with a narrow band of (core) data packages, which is able to achieve memory efficiency and is able to handle those level-set based algorithms evolving data locations considerably far from the surface.

Second example is to handle a industrial design of heat exchanger, as shown in Fig. 7, where very small geometry features are presented near the connections between pipes and plates. Note that, while these features can be captured accurately using high-resolution level-set field, they are not necessarily for a optimization-oriented numerical simulations focusing on the heat-transfer performance only.

<p align="center"><img src="https://arxiv.org/html/2512.11473v1/x2.png" 
alt="Heading Figure" sstyle="max-width:100%; height:auto;">
<center> Fig. 7. Geometric model of a heat exchanger. (left) the entire model. (right) a cross-section view indicting very small features.  </center> </p>

As shown in Fig. 8, these small features are removed by applying the cleaning algorithm. With the present method, this operation can be carried out on GPU, and decreases computation time greatly (with about 5 times speed up). As shown in the finally generated SPH particles, a regular particle distribution without single layer or filament of particles generated and hence ensure stable SPH simulations.

<p align="center"><img src="https://arxiv.org/html/2512.11473v1/x3.png" 
alt="Heading Figure" sstyle="max-width:100%; height:auto;">
<center> Fig. 8.SPH particles for a heat exchanger model. (left) the entire model. (right) a detail-view near the region where small features are cleaned.</center> </p>

## Limitations and Future Work

It is noteworthy that, compared previous sparse-grid method, The present approach marks an entire data package as activated if its center overlaps with the body surface, regardless of the actual volume intersected. Such tagging process can be further optimized to reduce memory redundancy while maintaining computational accuracy. Another future work would be developing GPU-enabled sign-distance-from-triangle-mesh algorithm, so that the full sparse-grid method can be run on GPU for maximum computational efficiency.

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
