---
layout: post
title:  "Misinterpretation of tensile instability in SPH solid dynamics (Part 3)"
date:   2024-09-25
categories: sph method
---
<head> {% include mathjax.html %} </head>

Xiangyu Hu

## Numerical examples

In this part, we demonstrate that the typical numerical simulations
which identified in the literature suffering from tensile instability,
are numerically stable, i.e. without particle clustering, artificial fracture
or zigzag stress pattern, if the non-hourglass formulation
(denoted as SPH-ENOG) introduced in the last part is employed.
In the following,
the numerical results obtained by the original SPH (denoted as SPH-OG)
method and artificial stress (denoted as SPH-OAS) method are also given as reference.

### 2D oscillating plate

<p align="center"><img src="{{site.baseurl}}/assets/img/oscillation-beam-setup.jpg" alt="Oscillation-plate-setup" height="150"/>
<center>Fig. 3. 2D oscillating plate: case setup.</center> </p>

As shown in Fig. 3, a oscillation 2D plate with one edge fixed is considered.
The length and thickness of the plate are $L$ and $H$ respectively, and the left part is fixed to produce a cantilever plate.

<p align="center"><img src="{{site.baseurl}}/assets/img/oscillation-beam.jpg" alt="Oscillation-plate" height="350"/>
<center>Fig. 4. 2D oscillating plate:
Evolution of particle configuration with time (t=0.05, 0.37 and 0.67)
for (a) SPH-OG, (b) SPH-OAS, and (c) SPH-ENOG.</center> </p>

It is observed form Fig. 4 that SPH-OG leads to severe numerical instabilities
(non-physical fractures and zigzag patterns in Fig. 1).
Numerical fractures occur just at the beginning of the simulation;
the zigzag particle distribution and the non-uniform profile of von Mises stress
indicate the hourglass modes. It is also shown that SPH-OAS, in which the non-physical fractures can be suppressed. However, the zigzag patterns still occur and became visually evident. Clearly, for the results obtained with non-hourglass formulation,
neither non-physical fractures and zigzag patterns appear.
Furthermore, the particle distribution is still uniform, and the stress profile is smooth.

<p align="center"><img src="{{site.baseurl}}/assets/img/oscillation-beam-decay.jpg" alt="Oscillation-plate-decay" height="350"/>
<center>Fig. 5. 2D oscillating plate:
The decay of oscillation magnitude obtained from SPH-OAS and SPH-ENOG.</center> </p>

As shown in Fig. 5, the simulation has been carried out for over 30 periods.
With the present SPH-ENOG, the deflection only decreases marginally,
but the SPH-OAS exhibits rapid energy decay,
thus it cannot be used for long-duration computations.

### 2D colliding rubber rings

<p align="center"><img src="{{site.baseurl}}/assets/img/two-ring--setup.jpg" alt="Colliding-ring-setup" height="350"/>
<center>Fig. 6. 2D oscillating plate: case setup.</center> </p>

The collision of two rubber rings is simulated. As shown in Fig. 6,
two rings with inner radius 0.03 and outer radius 0.04 are moving towards each other
with a initial velocity. Note that, when two rings collide with each other,
a significant tensile force will be generated.

<p align="center"><img src="{{site.baseurl}}/assets/img/two-ring-easy.jpg" alt="Colliding-ring-easy" height="350"/>
<center>Fig. 7. Evolution of particle configuration with time
(t = 0.002, 0.005, 0.008 and 0.012) for 2D colliding rubber rings.
The results are obtained by different SPH methods,
i.e., SPH-OG (left column), SPH-OAS (middle column),
and SPH-ENOG (right column).</center> </p>

Fig. 7 shows the evolution of particle configuration for the SPH-OG, SPH-OAS and the present SPH-ENOG when the initial velocity magnitude is moderate. Clearly, the SPH-OG suffers from serious non-physical fractures and zigzag patterns at the beginning of the computation,
and the calculation process can barely continue. For the SPH-OAS, the non-physical fractures can be suppressed, and the particle distribution is uniform at the initial stage. However, with the passage of time, the zigzag distribution of particle configuration and von Mises stress gradually becomes apparent. While for the SPH-ENOG, the particle and stress distribution are uniform during the whole calculation process, and the numerical instabilities can be completely removed.

<p align="center"><img src="{{site.baseurl}}/assets/img/two-ring-hard.jpg" alt="Colliding-ring-hard" height="350"/>
<center>Fig. 8. Evolution of particle configuration when the collision velocity is high.</center> </p>

Furthermore, the initial velocity is increased. It can be seen from Fig. 8,
not only zigzag patterns, but also numerical fractures appear for the SPH-OAS.
Fortunately, the present SPH-ENOG performs well even at such large initial velocity,
and all numerical instabilities can be perfectly eliminated,
which suggests the stability and robustness of the present SPH-ENOG.

### No tensile instability in SPH elastic dynamics?  

In these two numerical example, even tensile stress is obvious and dominant the dynamics,
the simulation is very stable using hourglass formulation.
Such observations at least is able to clarify that the numerical instability
in the standard SPH elastic solid dynamics is not introduced by tensile instability.
One more brave assumption is that there may be no tensile ability at all.

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
