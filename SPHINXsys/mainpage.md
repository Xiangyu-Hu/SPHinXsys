SPHinXsys (pronunciation: s'finksis)
is an acronym from <b>S</b>moothed <b>P</b>article
<b>H</b>ydrodynamics for <b>in</b>dustrial comple<b>X</b> <b>sys</b>tems.
It provides C++ APIs for physical accurate simulation and aims to model coupled
industrial dynamic systems including fluid, solid, multi-body dynamics and
beyond with SPH (smoothed particle hydrodynamics), 
a meshless computational method using particle discretization.

Included physics
-----------------
Fluid dynamics, solid dynamics, fluid-structure interactions (FSI), 
and their coupling to multi-body dynamics (with SIMBody library https://simtk.org) 

SPH method and algorithms
-----------------
SPH is a fully Lagrangian particle method, 
in which the continuum media is discretized into Lagrangian particles
and the mechanics is approximated as the interaction between them
with the help of a kernel, usually a Gaussian-like function. 
SPH is a mesh free method, which does not require a mesh to define 
the neighboring configuration of particles, 
but construct of update it according to the distance between particles.
A remarkable feature of this method is that its computational algorithm 
involves a large number of common abstractions 
which link to many physical systems inherently. 
Due to such unique feature, 
SPH have been used here for unified modeling of both fluid and solid mechanics. 

The SPH algorithms are based on the published work of the  authors.
The algorithms for the discretization of the fluid dynamics equations 
are based on a weakly compressible fluid formulation, 
which is suitable for the problems with incompressible flows, 
and compressible flows with low Mach number (less than 0.3).
The solid dynamics equations are discretized by a total Lagrangian formulation,
which is suitable to study the problems involving linear and non-linear elastic materials.
The FSI coupling algorithm is  implemented in a kinematic-force fashion, 
in which the solid structure surface describes the phase-interface and, 
at the same time, experiences the surface forces imposed 
by the fluid pressure and friction.
 
Material models
-----------------
Newtonian fluids with isothermal linear equation of state. Non-newtonian fluids with Oldroyd-B model.
Linear elastic solid, non-linear elastic solid with Neo-Hookian model and anisotropic muscle model.

Authors
-----------------
Xiangyu Hu, Luhui Han, Chi Zhang, Shuoguo Zhang, Massoud Rezavand

Project Principle Investigator
-----------------
Xiangyu Hu

Acknowledgements
-----------------
German Research Fundation (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1 and HU1527/12-1.

Please cite
-----------------
Luhui Han and Xiangyu Hu, "SPH modeling of fluid-structure interaction", Journal of Hydrodynamics, 2018: 30(1):62-69.