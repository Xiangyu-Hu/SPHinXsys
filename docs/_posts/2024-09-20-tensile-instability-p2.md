---
layout: post
title:  "Misinterpretation of tensile instability in SPH solid dynamics (Part 2)"
date:   2024-09-20
categories: sph method
---
Xiangyu Hu

## Rational

Since the method for preventing tensile instability is still not able to eliminate the zigzag artifact (characteristic for hourglass modes),
it is quite rational to propose that both instabilities have contributions.
However, this does not exclude the possibility that only hourglass modes contribute
the instability observed in numerical simulations.

If only the hourglass modes have contribution, one may require that
all the numerical simulations which are identified as tensile instability
become numerically stable if only the method preventing hourglass is applied.

<p align="center"><img src="{{site.baseurl}}/assets/img/zigzag.jpg"alt="Zigzag" height="350"/>
<center>Fig. 2. Illustration of zero energy modes by considering a simple case.
The particles are uniformly distributed along the x-axis,
and change from time step n to time step n + 1
with the velocity perpendicular to the x-axis.</center> </p>

### Origin of hourglass modes

As shown in Fig. 2, we use a simple case to demonstrate the hourglass
(or zero energy) modes.
The particles are uniformly distributed along the $x$-axis,
and the particles change from time step $n$ to time step with the velocity.
Under these circumstances, the original SPH calculations yield
a zero velocity gradient for particle $i$ based on stand SPH discretization,
consequently leading to a negligible stress increment.
In other words, this deformation mode is not resisted and,
naturally, cannot recover even under elastic deformation conditions,
which leads to an inaccurate estimation of the shear stress and shear acceleration.

### Remedy for hourglass modes

An interest phenomenon one can observe is that, as shown in Fig. 1c,
although the artificial stress can eliminate non-physical fractures in some cases,
the artifacts with zigzag patterns appear in the stress profile.
This issue is not unique for the specific cases,
and can be found in other ULSPH methods aimed at mitigating tensile instability.
The zigzag pattern is a typical numerical instability phenomenon well known as hourglass modes,
another fundamental numerical issue related to collocation methods.
Hourglass modes were also found in total Lagrangian SPH (TLSPH) elastic dynamics,
when the deformation is very large, even though TLSPH does not suffer from the tensile instability,
i.e. without particle clustering and non-physical fractures.

### Origin of the numerical instability

Such observation leads to people puzzling which is the original mechanism of the instability.
Is the tensile instability or actually the hourglass modes?.
In the next part, I will introduce the approach in the above journal paper
which is able to resolve this puzzle.

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
