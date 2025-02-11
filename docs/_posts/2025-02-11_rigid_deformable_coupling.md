---
layout: post
title:  "Rigid-deformable solid coupling validation"
date:   2025-02-11
categories: test examples
---
Xiangyu Hu and Weiyi Kong

## Introduction
Under some conditions, the stiffer part of a composite-material body can be considered rigid. It could save a lot of computational time if the problem can be solved as a rigid body coupled with a soft elastic body.

The coupling is implemented in SPHinXsys in a simple way: modeling the rigid and the elastic parts as a single body and apply rigid motion constraints to the rigid part using Simbody. The inner forces computed by 1st part integration will be used as the external forces acting on the rigid part.

## Validation case description

The coupling method is verifed by the twisting of an elastic bar coupled with a rigid cube at the end, as displayed in Fig. 1. The length of the elastic and rigid bar is 4m and 1m, respectively, while the height and width are both 1m. The left end of the bar is fixed. A total force $\mathbf{F} = -p(t) \mathbf{e_y}$ and a torque $\mathbf{\tau} = -p(t)h\mathbf{e_x}$ is applied to the rigid body, where $\mathbf{e_x}$ and $\mathbf{e_y}$ denote the unit vector in x and y direction respectively. The force $p(t)$ takes the form:

$$
\begin{equation}
  p(t)=\begin{cases}
    0.05t, & \text{if $t<1$}.\\
    1, & \text{if $t \ge 1$}.
  \end{cases}
\end{equation}
$$

<p align="center"><img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/geometry.svg" alt="rigid-deformable_coupling" height="100"/>
<center>Fig. 1. Geometry of 3-dimensional elastic bar coupled with a rigid cube.</center> </p>

The coordinates of the points are listed in the table below:

| Coordinate [m] | A | B | C | D | E  | F  | G   | H  | I  | J  | K  | L |
|----------------|---|---|---|---|----|----|-----|----|----|----|----|---|
|x               |-2 |2  | -2|-2 | 2  | 2  | -2  |-2  |3   |3   |3   |3  |
|y               |0.5|0.5|0.5|0.5|-0.5|-0.5|-0.5 |-0.5|-0.5|-0.5|0.5 |0.5|
|z               |1  |1  |0  |0  |1   |0   |0    |1   |1   |0   |1   |0  |
<center>Tab. 1. Coordinates of the coupled bar.</center> </p>

## Results

To validate the coupling method, the data is compared with the FEM simulation results obtained by by [FEBio](https://help.febio.org/FEBioTheory/FEBio_tm_3-4-Section-7.10.html). A static solver is utilized in the simulation, with a time step size of 0.05 and total time step of 20. As the static solver is not available in SPHinXsys, an artificial damping is utilized to accelerate the process to reach the steady state. 

Displacements of the four corner points on the coupling interface are used for comparison. Observer 0, 1, 2 and 3 correspond to point C, B, F and E, respectively.

### Convergence analysis
The convergence of the simulation results with mesh/resolution refinement is displayed in the figures below.

<p align="center">
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_0_dimension_x_disp-res.svg" alt="obs_0_x" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_0_dimension_y_disp-res.svg" alt="obs_0_y" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_0_dimension_z_disp-res.svg" alt="obs_0_z" height="200"/>
</p>
<center>Fig. 2. Displacement at observer 0 against total node/particle number.</center> </p>

<p align="center">
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_1_dimension_x_disp-res.svg" alt="obs_1_x" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_1_dimension_y_disp-res.svg" alt="obs_1_y" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_1_dimension_z_disp-res.svg" alt="obs_1_z" height="200"/>
</p>
<center>Fig. 3. Displacement at observer 1 against total node/particle number.</center> </p>

<p align="center">
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_2_dimension_x_disp-res.svg" alt="obs_2_x" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_2_dimension_y_disp-res.svg" alt="obs_2_y" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_2_dimension_z_disp-res.svg" alt="obs_2_z" height="200"/>
</p>
<center>Fig. 4. Displacement at observer 2 against total node/particle number.</center> </p>

<p align="center">
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_3_dimension_x_disp-res.svg" alt="obs_3_x" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_3_dimension_y_disp-res.svg" alt="obs_3_y" height="200"/>
  <img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/obs_3_dimension_z_disp-res.svg" alt="obs_3_z" height="200"/>
</p>
<center>Fig. 5. Displacement at observer 3 against total node/particle number.</center> </p>

With the same element/particle size, the total number of particles in SPHinXsys is higher than the total nodal number of FEBio, as the rigid cube is modelled by eight nodes in FEBio, while in SPHinXsys it is fully discretized by particles.

The convergence trend of FEBio and SPHinXsys is similar, though one reaches the converges from undershoot while the other converges from the overshoot. The results with finest mesh/resolution, i.e., an element/particle size of 0.025m, are used for comparison.

### Comparison

The displacements with a resolution of 0.025m at the four observation points are listed in Table 2.

| observer | SPHinXsys | FEBio |Error [%] |
|-----------|-----|-------|------------|
| 0 | [-0.52749511, -2.61691304, 0.10558404] | [-0.504557,-2.60476,  0.0948022] | 1.06 |
| 1 | [-0.27093921, -2.19197496, -0.03331926] | [-0.244248 , -2.16765  , -0.0442815] | 1.73 |
| 2 | [-1.22377339, -2.13595798,  0.56736396] | [-1.21962 , -2.11663 ,  0.570901]  | 0.80 | 
| 3 | [-0.96827137, -1.7113192 ,  0.43213959] | [-0.959311, -1.67953 ,  0.431817]  | 1.67 |
<center>Tab. 2. Displacements obtained by SPHinXsys and FEBio at observation points.</center> </p>

The error is defined as $|\mathbf{d_{SPH}}-\mathbf{d_{FE}}| / |\mathbf{d_{FE}}|$. The maximum error is within 2%, which proves the consistency between SPH and FEM.

### Performance analysis
By using a physical viscosity $\eta=2\eta_0$, where $\eta_0=\frac{\alpha}{4}\sqrt{\rho E}h$, $\alpha=0.4$, the SPH simulation reaches the steady state at about 4.6s. The simulation is terminated if the difference of displacement at observer 0 between two different output steps is smaller than a threshold $\epsilon$ in five consecutive output steps, i.e., $|\mathbf{d_0^n}-\mathbf{d_0^{n-1}}|<\epsilon, \quad n=N,N-1,\cdots, N-5$. The threshold $\epsilon$ is set to $0.005m$. The loop time of SPH is compared with the total elapse time of FEbio. As shown in Table 3 and Figure 6, the computational time of SPHinXsys is longer than that of FEbio.

|  Total node number of SPH [-] | Total particle number of FEbio |  end time of SPH [s] | loop time of SPHinXsys [s]| total elapse time of FEbio [s] |
|------------------------|-------------------|------------------------------|--------------------------------|--------------------------------|
| 2624 | 2681 | 4.60735|17.6356|2.81975|
| 8784 | 8289 | 4.60628|38.4047|8.95845|
| 20736 | 18793 | 4.60363|105.181|25.0419|
| 40400 | 35729 | 4.60286|255.339|64.9546|
| 69696 | 60633 | 4.60267|539.735|130.239|
<center>Tab. 3. Computational time of SPHinXsys and FEBio.</center> </p>

<p align="center"><img src="{{site.baseurl}}/assets/img/rigid_deformable_coupling/performance_analysis.svg" alt="rigid-deformable_coupling" height="300"/>
<center>Fig. 6. Computational time of SPHinXsys and FEBio against total node/particle number.</center> </p>

## Conclusions

We presents the realization of rigid-deformable coupling of solid for SPHinXsys with a simple method. The algorithm is validated by comparison with FEM results from FEBio. Although the computational time is longer than FEM due to lack of static solver, the error of displacements is less than 2%, which proves the feasiblity of this method.

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