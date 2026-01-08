---
layout: post
title:  "Contiguous Storage of Grid Data for Heterogeneous Computing (Part 3)"
date:   2026-01-08
categories: high performance computing 
---
Xiangyu Hu and Fan Gu

<p align="center"><img src="{{site.baseurl}}/assets/img/multi-resolution-level-set.png" alt="Heading Figure" sstyle="max-width:100%; height:auto;">
<center> Fig. 3. Multi-resolution level-set field shown with successively doubled resolutions from top to boottom and left to right. </center> </p>

## Mesh Local Dynamics

The class of a specific mesh local dynamics is decomposed into two components: variable management and a computing kernel. Variables involved in an operation must either be stored as copyable small objects or managed by the variables mentioned above. While the small objects are copied directly to the device alongside computing kernel instance, larger data structures are managed via the functions of `DiscreteVariable`, as will be discussed in the next section, which ensures correct allocation and synchronization across host and device.

Each local dynamics has its own dedicated computing kernel definition as a nested class, ensuring modularity and facilitating kernel reuse. Therefore, every computational operation is thus encapsulated in a class form, enabling the implementation class to construct an executable instance.

During initialization, the mesh local dynamics object is instantiated and configured with all required variables. However, the computing kernel is only initialized at the first request, in which the host-device transfer of data takes place according to the execution policy. Computing is then carried out transparently on the target platform as determined by the policy.

### Implementation Class and Data Transfer

An implementation class is responsible for the computing kernel transfer and provide the mesh dynamics with a valid computing kernel based on execution policy. As shown in the previous code snippets, the computing kernel will be created in the implementation object by a getter function at the first call. According to the execution policy chosen, the implementation object is responsible to copy the computing kernel inside local dynamics instance to the device to execute. As computing kernels take only copyable objects or raw variable pointers, the kernel can be copied to the device directly with no extra attention.

To enable computation on GPU, all relevant data must be resident on the device in advance. In SYCL, two primary memory access mechanisms are available. The first utilizes `sycl::buffer` and accessors. SYCL automatically manages data transfers between host and device, ensuring that data is available when accessed. This implicit data management simplifies development but can reduce control and transparency.

The present work adopts mainly the second unified share memory (USM) mechanism, due to its explicit control and pointer access can be used in the same way on host and device with plain `cpp`. The pointer-based nature of USM simplifies the function interface, as only raw pointers need to be handled. This eliminates the need for buffer management in user-defined functions, allowing them to work directly without aware of the SYCL types and functions.

For data arrays that requires updates, such as mesh, background-mesh and meta variables, they are managed by the `DiscreteVariable` class. Each `DiscreteVariable` includes a `DeviceOnlyVariable` instance that manages the device-side allocation, transfer, synchronization, and de-allocation, ensuring correct memory management and preventing leaks. Note that, device-side memory is only allocated and data transferred when the device pointer is first requested. Again the latter is generated only when a computing kernel is first requested by the loop function mesh_for package_for, as shown in the previous code snippets, with device execution policies.

For singular global parameters, such as total number of data packages, memory is allocated using `malloc_shared` from the USM mechanism. This shared memory is accessible from both host and device, with synchronization handled automatically by SYCL.

## Multi-resolution Level-set

The level set is constructed from an input geometry (STL file) using a multi-resolution approach, as shown in Fig. 3, in order to ensure water-tight surface when the input geometries are leaking due to the topological inconsistency. Given a target resolution, several layers with successively doubled resolutions are established. Each layer corresponds to a distinct spatial resolution as an instance of `MeshWithGridDataPackages`, which encapsulates both level-set data and metadata for that layer. The coarsest mesh is initialized first, and subsequent layers that each with doubled resolution are initialized based on the preceding coarser layer. The geometric shape and boundary conditions remain consistent across all layers; but the finer layer contains 4 (2D) or 8 (3D) times as many logical cells as its coarser counterpart. This successive initialization continues until the desired resolution is achieved. The number of required layers is computed in advance.

### Initialization Process

The present initialization process is mainly carried out on CPUs due to the present signed-distance function from a triangle mesh 
[(Bærentzen at al. 2005)](https://doi.org/10.1109/VG.2005.194111), lacking of GPU implementation. Specifically, the initialization process is composed of the following steps.

* Initial cell tagging for the core data packages. For the coarsest layer, the distance from each cell of the background mesh to the surface is evaluated and whose distances are less than the cell size is activated. For successive layers, only the cells covered by the coarse core cells will be evaluated. Using Intel’s TBB (Threading building Blocks) concurrent vector, the cell index and packages type are inserted.

* Tagging the cell for inner data packages. This is done by checking if any nearest neighbor of a mesh cell activated by a core data package. Again, the cell index and packages type are inserted to the concurrent vector.

* Sorting for locality. During tagging, activated cells are identified in parallel. This results that the indexes of activated cells are in arbitrary order, which may degrade cache locality. To address this, a sorting step is introduced to reorder cell indexes according a predefined sequence, typically based on spatial location. After sorting, the cell indexes and data-package types are then copied to newly allocated meta variables. So that they can be used later in GPU computing.

* After that the level-set value is evaluated for all data packages using the sign-distance function from triangle mesh [(Bærentzen at al. 2005)](https://doi.org/10.1109/VG.2005.194111).

* Finally, the neighborhood of all activated cells are defined. Note that, in order to achieve full consistency for data queries, the far-field packages’ neighbors are set to themselves.

Note that, after the level-set value is evaluated for all data packages, the computation can be carried out either on CPU or GPU according ro the chosen execution policies. Also note that, the level-set value evaluation can be carried before the sorting, which then can be run on GPU but need to reorder the level-set data packages evaluated in arbitrary order.

### Consistency Correction and Small Feature Cleaning

Due to the possible topological inconsistency of STL file, directly application of signed-distance function [(Bærentzen at al. 2005)](https://doi.org/10.1109/VG.2005.194111) to identify the contain condition (or the sign) of level set may go wrong. Such issue may finally leads to the generation of SPH particles wrong. To avoid this issue, only the sign of level set for those data points very close to the surface is directly used, those at other locations are obtained by a two-step diffusion process from the near interface to the entire domain. The first coarse step is on the mesh cells and the second refined one is on the data packages. Note that only on the coarsest layer all cells are evaluated. For refined layers, the operation is limited to inner cells or data packages.

Very often the geometry which is originally generated for manufacturing includes many small features which are not necessary for computational fluid or solid dynamics (CFD or CSD) simulations and may lead to numerical instabilities if not cleaned. In the present work, we reimplemented the level-set cleaning algorithms (only on the finest layer) in [(Yu et al. 2023)](https://doi.org/10.1016/j.cpc.2023.108744) so that it can be run on GPU.

Note that, these algorithms heavily replies on the indirect access to neighbor data packages. As shown in the following,

```cpp
template <int PKG_SIZE, typename RegularizeFunction>
Vec2d regularizedCentralDifference(
PackageData<Real, PKG_SIZE> *input, const CellNeighborhood2d &neighborhood,
const Array2i &data_index, const RegularizeFunction &regularize_function)
{
    DataPackagePair center =
    NeighbourIndexShift<PKG_SIZE>(data_index, neighborhood);
    DataPackagePair x1 = NeighbourIndexShift<PKG_SIZE>(
        data_index + Array2i(1, 0), neighborhood);
    DataPackagePair x2 = NeighbourIndexShift<PKG_SIZE>(
        data_index + Array2i(-1, 0), neighborhood);
    DataPackagePair y1 = NeighbourIndexShift<PKG_SIZE>(
        data_index + Array2i(0, 1), neighborhood);
    DataPackagePair y2 = NeighbourIndexShift<PKG_SIZE>(
        data_index + Array2i(0, -1), neighborhood);
    Real dphidx_p = input[x1.first](x1.second) - input[center.first](center.second);
    Real dphidx_m = input[center.first](center.second) - input[x2.first](x2.second);
    Real dphidx = regularize_function(dphidx_p, dphidx_m);
    Real dphidy_p = input[y1.first](y1.second) - input[center.first](center.second);
    Real dphidy_m = input[center.first](center.second) - input[y2.first](y2.second);
    Real dphidy = regularize_function(dphidy_p, dphidy_m);
    return Vec2d(dphidx, dphidy);
} 
```

a typical regularized central difference scheme for gradient evaluation can be implemented with the help of `NeighbourIndexShift`. Also note that efficient high-order finite difference can be extended easily due to the relative large size of the data package.

### Integral with SPH Smoothing Kernel and Grid-Particle Coupling

Another algorithm in [(Yu et al. 2023)](https://doi.org/10.1016/j.cpc.2023.108744) is computing the integral with SPH smoothing kernel on the sparse-grid. This is essential for the physical relaxation of the SPH particles  [(Zhu et al. 2021)](https://link.springer.com/article/10.1007/s42241-021-0031-y), which is also reimplemented for GPU execution. The SPH particle relaxation is a typical grid-particle coupling algorithm in which the integral field is interpolated by bi- or tri-linear interpolation to the particle’s position and used in the form of surface force to drive the particle. Note that, the interpolation operation is different from the previous local mesh dynamics as it, other than looping on the cells or data packages, but evolves position-based random memory access of a data package and may be its neighbors, which again heavily relies on the usage of `NeighbourIndexShift`.

In the next post, 
I will continue to give the details on performance evaluation and numerical examples.

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
