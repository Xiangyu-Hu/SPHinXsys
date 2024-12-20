---
layout: post
title:  "Data type handling in SPHinXsys for heterogeneous architecture"
date:   2024-12-20
categories: sphinxsys update
---

The core design idea of SPHinXsys for heterogeneous architecture
is encapsulating the computing-intensive operations into computing kernels
which can be executed on different devices including CPU and GPU.

Although the idea is straightforward, it relies on the an important assumption
that the codes written in computing kernels should work on all devices.
This is certainly the case when only simple data types are used.
When more complex data types are involved,
as it is the case for SPHinXsys due to third-party libraries,
one need to be careful since these types most probably
are not as robust as the simple data types.

Therefore, in SPHinXsys computing kernels we propose the following guidelines:

1. Use Eigen3 math vectors and matrixes as much as possible.
2. Keep the right-hand-side data type homogenous in each expression.
3. Use explicit type conversion, and only converting one data type at a time.

Here, we provide an example to illustrate the above guidelines.
The following is a computing kernel function:

<https://github.com/Xiangyu-Hu/SPHinXsys/blob/66b9bfefd6e23088d47e1984f8435e8d09ba300e/src/shared/shared_ck/particle_dynamics/solid_dynamics/solid_constraint.hpp#L114-L125>

One can see that there quite some data types are involved.
We first look at the left hand side.
On line 119, the type of the declared variable `Vecd` is Eigen3
math vector with 2 or 3 components.
From line 120 to line 122, `SimTKVec3`,
the 3D math vector used in Simbody library.
On line 123, the return value type `SimTK::SpatialVec` is another
Simbody data type as the combination of torque and force acting on a moblized body.

There are more data types and conversions involved on the right hand side.
On line 119, the data type is device-only unified shared memory (USM) pointers of
`DisreteVariable<Vecd>`.
On lines 120 and 121 there are 4 explicit type conversions:
2 from `Vecd` and that pointed by the USM pointer,
and 2 from `Vec3d` to `SimTKVec3`.
On line 121, beside two conversions,
it involves another value pointed by a device-shared USM pointer for
`SingularVariable<SimTKVec3>`.
On line 122, there is an operation (cross production) for `SimTKVec3`.
On line 123, the constructor of type `SimTK::SpatialVec` is called.

We found that this kernel works well for CPU and SYCL targeting Nvidia GPU
but crash when using SYCL but targeting intel CPU as device,
suggesting there must be something not OK.
Therefore, we have tried to simplifying the operations
according the guidelines. Finally, we have obtained the following code:

<https://github.com/Xiangyu-Hu/SPHinXsys/blob/4e91c4fdc80b79f4633ba9c5aa8e5f8a1b79c53b/src/shared/shared_ck/particle_dynamics/solid_dynamics/solid_constraint.hpp#L114-L125>

Now, the USM data only use Eigen3 math vector `Vecd`,
`SimTKVec3` is not used anymore, the right-hand-side is homogenous
and conversions are explicit and not mixed with other operations.
After these changes, the code works now for all.

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
