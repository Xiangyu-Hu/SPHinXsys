---
layout: post
title:  "Towards heterogeneous parallelism for SPHinXsys (Part 1)"
date:   2025-01-30
categories: high performance computing 
---
Xiangyu Hu and Alberto Guarnieri

## Introduction

This blog is based on the presentation of the same title
to be presented on [19<sup>th</sup> SPHERIC World Conference 2025](https://spheric2025.upc.edu/).

As smoothed particle hydrodynamics (SPH) is a typical particle-based method,
many SPH libraries suffer from high computational costs
due to intensive particle interactions.
Compared to the multi-cores CPU system,
the many-cores device (typically GPU) system
provides a much higher performance/cost ratio
for the intermediate (multi-million particles) scale simulations,
which are frequently encountered in industrial applications.
On the other hand, since the versatile functionality of CPU system,
it is expected that SPH libraries are using heterogeneous parallelism
so that they are able to take advantage of the both.

Currently, some SPH libraries already offer GPU support,
either through vendor-specific implementations,
such as [DualSPHysics](https://dual.sphysics.org/),
or heterogeneous programming using the OpenCL standard,
such as [AQUAgpusph](http://canal.etsin.upm.es/aquagpusph/)
and [pySPH](https://pysph.readthedocs.io/).
An important issue of these implementations is
that they require computing kernels to be defined,
as a low-level specification,
with specific directives separated from the rest of the code,
rendering it less easily programmable and maintainable.

## Open-source multi-physics library SPHinXsys

[SPHinXsys](https://www.sphinxsys.org/) is an open-source C++ SPH multi-physics simulations library.
It addresses the complexities of fluid dynamics, structural mechanics,
fluid-structure interactions, thermal analysis,
chemical reactions and AI-aware optimizations.
SPHinXsys has offered CPU parallelism
using [Intel’s Threading Building Blocks (TBB) library](https://github.com/uxlfoundation/oneTBB)
and several features suitable for open-source based development.

First, since the parallel execution is encapsulated into
the low-level classes decoupled from the SPH method,
the developers, as typical mechanical engineers,
only need to work the numerical discretization
without concerning parallelization.
Second, variable testing approaches, including unit test, google test
and regression test, have been incorporated
so that refraction, adding new features and other modifications
can be implemented without the worries about the already established functionalities.
Third, automated cross-platform
[continuous integration/development (CI/CD) tests](https://github.com/Xiangyu-Hu/SPHinXsys/blob/master/.github/workflows/ci.yml)
have being carried out on the open-source development
platform frequently for any modification of
the main branch of SPHinXsys repository.
In addition, SPHinXsys makes the uses of several third particle libraries,
such as [Simbody library](https://github.com/simbody/simbody)
for multi-body dynamics,
[Pybind11 library](https://github.com/pybind/pybind11)
for generating  python interface used for machine learning and optimization applications.

## SYCL Programming Standard

Different from the other GPU-able SPH libraries,
our work is based on the [SYCL standard](https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html)
for parallel computing.
SYCL specification defines a new single-source
ISO C++ compliant standard with compute acceleration.
It aims to be an high level abstraction that
can be programmed as classical CPU code,
without accelerator-specific directive
and application interface (API) calls.
In this work, the Intel SYCL implementation,
i.e. [Data Parallel C++ (DPC++)](https://github.com/oneapi-src/DPCPP_Reference)
is chosen, because it is a part of Intel’s oneAPI suit
which also includes the TBB library already used in SPHinXsys.

A SYCL computing kernel
(sharing the same term with the smoothing kernel used in SPH)
is a scoped block of code executed on a SYCL device
under the global context `queue` and
the local one `command_group`
initiated by the host.
Note that, as the SYCL host and device both can be CPU,
CI/CD jobs can test SYCL kernels on computers without GPU,
such the standard runner provided by
[Github platform](https://docs.github.com/en/actions/using-github-hosted-runners/using-github-hosted-runners).

In the next post, I will give the details on the implementation of SYCL kernels in SPHinXsys.

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
