---
layout: post
title:  Conflict with vcpkg for SYCL from intel OneAPI"
date:   2025-06-17
categories: sphinxsys update
---

The SYCL compiler provided by `intel-oneapi-base-toolkit` from intel is updated routinely.
Many library files in the toolkit will also be updated.

As a consequence, the SYCL library may be related to some third-party libraries in vcpkg, such as TBB.
When installing the packages from vcpkg and using it as the tool chain for SPHinXsys building,
the old libraries used in the vcpkg building and installation may be considered as essential ones.
In my case, the `mkl` library in 2024 version leads to conflict with the updated version of SYL library
used by SYCL compiler, i.e. `mkl` library in 2025.1 version.  

To solve this conflict, I need to reinstall the vcpkg library with updated OneAPI,
so the new version of `mkl` library is used both for building the package and used in SYCL compiler.

Another solution may be installing the third-party package with static libraries,
so these libraries used for the third-party package are not needed anymore.

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
