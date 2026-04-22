# Context
In engineering applications, we often encounter problems with only a small portion of the model needs refinement. For example, the folding part of a sheath requires a high resolution due to the complicated geometry, whereas the other parts can be modeled with lower resolutions. If we use a uniform refined particle spacing over the whole domain, the computational cost will be unaffordable. Hence, we introduce a multi-resolution algorithm for volumetric solid. 

## SPH with adaptive resolution
The discretization scheme of SPH relies on the gradient of kernel function $\frac{\partial }{\partial r}W(r_{ij}, h_{ij})$. For a pair of particles with different smoothing length $h_{i}$ and $h_j$, the term $h_{ij}$ can be chosen as asymmetric, e.g., $h_{ij}=h_i$ (the gather formulation), or $h_{ij}=h_j$ (the scatter formulation), or symmetric, e.g., $h_{ij}=(h_i+h_j)/2$ (the average support). To be consistent with the fluid multi-resolution implemented in SPHinXsys, here we select a symmetric form $h_{ij} = \max(h_i, h_j)$.

Since the total Lagrangian form is used for solid solvers, the inner neighbor configurations only need to be built once in the initial frame. The multi-resolution particle distribution stays unchanged during the simulation, so no splitting-merging is considered.

With the properly defined adaptive kernel value, most algorithms used for solid can be directly applied to the multi-resolution object without changes.

## Multi-level cell-linked lists
The refinement resolution is decided by the refinement level $l$, with $dp^{l} = 2^{-l} dp^0$, where $dp^0$ is the largest particle spacing. For a highest refinement level, $l+1$ mesh levels are created, with grid size equal to the cut-off radius $\Delta^l = 2.3dp_{l}$. A particle with cut-off radius $R$ will be added to mesh level $l$ if it satisfies:

$$R - \Delta^l < eps$$

where $eps$ is a small value. The mesh level of a solid with refinement level $l=1$ is illustrated in the graph below. The refined particles with $dp=0.5dp^0$ fall in mesh level 1, while the coarse particles with $dp=dp^0$ and transition particles with $0.5dp^0 < dp < dp^0$ fall in mesh level 0.

![img1](./multi_cell_level.svg)

## Damping with adpative resolution
The most tricky part of solid adpative resolution is how to implement the operating splitting damping algorithm. The damping method used for uniform resolution is introduced in [Zhu2022](https://doi.org/10.1016/j.jcp.2022.111105). Since the damping is pairwise, the variable value of a particle i and its inner neighbor particle j needs to be modified together (see `Damping<Inner<Pairwise>, DataType, DampingRateType>::interaction()`). However, as the interaction function is running over all particle ids in parallel, other particles might try to modify the variable of j at the same time, leading to memory conflict. 

To solve this problem, we need to make sure that those particles running in parallel don't share any neighbors. This is achieved by dividing the cells into 9 (for 2d) or 27 (for 3d). As there is a distance of $2\Delta$ between splitting cells with the same index, particles in those cells won't have common neighbors, hence all the cells with same index i can run in parallel. 

![img2](./split_cell.svg)

Pseudo code:

```
for (int i in {0,1,2,3,4,5,6,7,8})
{
    auto data = cell[i].get_data_in_cell();
    parallel_exec(data);
}
```

However, for the cells of different mesh levels, it is difficult for us to know their position. A particle in split cell 1 of mesh level 0 and a particle in split cell 1 of mesh level 1 may share neighbors. As can be seen from the graph below, split mesh 1 of two different levels are adjacent to each other.

Moreoever, since the cut-off radius $R_{ij}$ is determined by $\max(R_i)$, two particles in the same splitting cell of level $l$ can share neighbors in mesh level $l-1$, e.g., particle i and j falling in split cell 1 share the same neighbor k.

![img2](./multi_level_split_cell.svg)

To fix this problem, we propose the following method:

*  A new inner relation is used for damping, where particle j is only added to the inner neighbor of i if it satisfies: (1) i and j are in the same mesh level and j > i, or (2) mesh level of j is higher than mesh level of i

 By this rule, a particle can only see neighbors with the same level spacing or smaller spacing, i.e., particle i and j are in the neighborhood of k, but k is not in the neighborhood of i and j, thus, i and j no longer share neighbors. To be consistent with the cross-level neighbors, particles of the same level should also only see each other once, hence, only particles with a larger index are included in the neighborhood.

*  Loop over split cells of different mesh levels in sequence to avoid adjacent cells between different levels. 

Pseudo code:
```
for (int l in {0,1,..,mesh_level_max})
{
    auto cell = get_cell(l);
    for (int i in {0,1,2,3,4,5,6,7,8})
    {
        auto data = cell[i].get_data_in_cell();
        parallel_exec(data);
    }
}
```

Note that as the new inner relation only has half of the neighbors, with the same physical viscosity, the speed to reach quasi-steady state will be slower than the old method. To speed up this process, I suggest users to increase the physical viscosity by twice.

## Detailed implementation

Those new classes are implemented for multi-resolution solid:

* `Integration1stHalfPK2RightCauchy`: the pair numerical damping used in `Integration1stHalf` is only applicable to uniform resolution, so a non-pairwise right Cauchy numerical damping used in earlier versions of SPHinXsys is added back
* `AdaptiveSplittingInnerRelation`: the inner relation used for damping, where particle pairs can only see each other once
* `void MultilevelCellLinkedList::particle_for_split(const execution::ParallelPolicy &, const LocalDynamicsFunction &)`: split cell loop algorithm for multi-level cell-linked lists

## Limitations
As the split cells of different mesh levels can only run in sequence, we will have $9(l+1)$ or $27(l+1)$ split cells with a refinement level of $l$, which means that the higher the refinement level is, the worse the parallelization of the algorithm becomes. As a result, despite the decrease in particle number, the multi-resolution spends a longer time in damping than the fully-refined uniform resolution case. For a large scale simulation, damping only takes up a small amount of computational time, so the low efficiency of damping might be acceptable, but this remains to be a critical problem.

For future studies, we may explore using atomic to prevent memory conflict at the interface cells of different mesh levels, which can enable us to run the split cells of different mesh levels in parallel.

## Examples

A test on 3d cantilever beam can be found in `3d_examples/cantilever_beam`. Simulation ran on 12th Gen Intel(R) Core(TM) i7-12700H with 20 cores. The reference deflection is taken as $-0.02 - 5.86 \times 10^{-3}$, which is the result with the conventional TLSPH method with h/dp=24 from [Zhu2022](https://doi.org/10.1016/j.jcp.2022.111105).

| Refinement level | Particle number [-]| dp_max [m] | dp_min [m] | Computational time [s]| Damping time [s]| Deflection [m] | Error [%] |
|-------|------------------|----------------------------|--------------------------------------|-----------------------|-----------------|-----------------|-----------------|
| 0 | 114688 | 0.00125 | 0.00125 | 564.12 | 135.11 | -0.02605 | 0.73 |
| 1 | 53203  | 0.0025  | 0.00125 | 295.93 | 78.96  | -0.02682 | 3.71 |
| 2 | 46175  | 0.005   | 0.00125 | 246.19 | 67.88  | -0.02867 | 10.86 |

In this example, the damping time of the multi-resolution solid is shorter than the fully-refined uniform one, but it could also take a longer time depending on the particle number and sizes.
