---
layout: post
title:  "Misinterpretation of tensile instability in SPH solid dynamics"
date:   2024-09-18
categories: sphinxsys update
---
Xiangyu Hu

## Misinterpretation of tensile instability in SPH solid dynamics (Part 1)

This blog is from some parts from the [Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2024.113072) paper "Essentially non-hourglass SPH elastic dynamics", Volume 510, 1 August 2024, 113072.

More than two decades ago, it was found the numerical simulations of typical elastic dynamics problems using updated Lagrangian smoothed particle hydrodynamics (ULSPH) may lead to numerical instability.

![Tensor-instability]({{ baseurl }}/assets/img/tensor-instability.jpg)
<center>Fig. 1. Illustration for (a) no numerical instabilities, (b) non-physical fractures and (c) zigzag pattern in ULSPH simulations of 2D oscillating plates. The particles are colored with von Mises stress.</center>

<br clear="left"/>

As shown in Fig. 1b for a typical example, the numerical instability involving particle clustering and non-physical fractures presented these simulations is often considered as tensile instability, a fundamental numerical issue caused by the tensile stress.

Since such numerical instability is identified, different approaches have been proposed to address this problem. However, these methods also have their own issues, such as failure to maintain conservation properties, low computational efficiency, or limited applicability, being only suitable for specific cases rather than universally applicable. One popular method is the artificial stress method, based on the concept of counteracting the tension to prevent the instability.

An interest phenomenon one can observe is that, as shown in Fig. 1c, although the artificial stress can eliminate non-physical fractures in some cases, the artifacts with zigzag patterns appear in the stress profile. This issue is not unique for the specific cases, and can be found in other ULSPH methods aimed at mitigating tensile instability. The zigzag pattern is a typical numerical instability phenomenon well known as hourglass modes, another fundamental numerical issue related to collocation methods. Hourglass modes were also found in TLSPH elastic dynamics, when the deformation is very large, even though TLSPH does not suffer from the tensile instability.

Such observation leads to people puzzling which is the original mechanism of the instability, i.e. is the tensile instability or actually the hourglass modes. In the next part, I will introduce the approach in the above journal paper which is able to resolve this puzzle.
