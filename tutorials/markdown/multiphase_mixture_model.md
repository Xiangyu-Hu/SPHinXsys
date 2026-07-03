## Governing equations in the Lagrangian frame moving with average velocity

In this work, we consider an incompressible multiphase flow using weakly compressible fluid model in the average velocity frame. The governing equations consist of the mass conservation and momentum conservation equations for each phase, as well as the equation of state for the mixture.

### Primary and derived variables
We consider two sets of variables, one is for the mixture, the other is for each phase.
For the mixture, we use compression ratio $\beta$ as a primary variable, which is defined as the ratio of the initial volume element to the current volume element, that is $\beta = \frac{V_0}{V}$, where $V_0$ is the initial volume and $V$ is the current volume.
For each phase $k$, we use the phase volume fraction $\phi_k$, 
with $\sum_k \phi_k = 1$, and the phase velocity $\mathbf{v}_k$ as the primary variables.

The mixture density $\rho$ can be computed from these primary variables as
$$\rho = \beta\sum_{k} \phi_k \rho^{o}_k$$
where $\tilde{\rho}_k = \beta\phi_k \rho^{o}_k$ and $\rho^{o}_k$ are the partial desnity and the reference density, respectively, of phase $k$. Also $\rho^{o} = \sum_{k} \phi_k \rho^{o}_k$ is the reference density of the mixture, which gives the mixture density of the initial volume element but with the present volume frection. 
The mass-averaged velocity $\mathbf{v}$ is defined as
$$\mathbf{v} = \frac{\sum_{k} \phi_k \rho^{o}_k \mathbf{v}_k}{\sum_{k} \phi_k \rho^{o}_k}$$
If we assume that all phases share the same speed of sound, then the pressure can be computed from the mixture density using the equation of state for the mixture, which is given by
$$p =  \rho^{o} c^2 (\beta - 1)$$
where $c$ is the artificial (reference) speed of sound. Furthermore, we assume all phase share the same pressure, which is a common assumption for multiphase flow models. Note that, under the weakly compressible fluid model, $c$ is chosen so that the density variation is less or about 1\%.  We assume that the inital mass of the volume element $m$ as
$$m = \rho_0 V_0$$
which is dependent of reference density and initial element volume only. As will be shown latter the mass of the fluid element is an invariant under the present average velocity moving frame.

For each phase, the density $\tilde{\rho}_k$ can be computed from the phase volume fraction and the reference density as $\tilde{\rho}_k = \beta\phi_k \rho^{o}_k$. The momentum of each phase can be computed as $\mathbf{m}_k = \tilde{\rho}_k \mathbf{v}_k = \beta\phi_k \rho^{o}_k \mathbf{v}_k$. Note that one can obtain the mixture velocity from the phase velocity as $\mathbf{v} = \frac{1}{\rho} \sum_{k} \beta\phi_k \rho$, which is consistent with the definition of the mass-averaged velocity due to the cancellation of $\beta$.

### Governing equations
For the mixture, one has the evolution of the compression ratio $\beta$, which is given as
$$\frac{d\beta}{dt} = -\beta \nabla \cdot \mathbf{v}$$
where the material derivative is defined as $\frac{d}{dt} = \frac{\partial}{\partial t} + \mathbf{v} \cdot \nabla$. 
Note the evelution of compression is purely kinematic and geoemtric and no physics involved.

The conservation equations of a Lagrangian volume element for the weakly compressible multiphase flow model,
two sets of equations can be derived, one is for each phase, the other is for the mixture. 

For each phase, first is the mass conservation equation, which is given as
$$\frac{d}{dt}\int_{V(t)}\tilde{\rho}_k\,dV + \int_{\partial V(t)}\mathbf{Q}_k\!\cdot\!\mathbf{n}\,dS = 0 $$
which describes the drift of each phase relative to the mixture.
Here, $\mathbf{Q}_k = \tilde{\rho}_k(\mathbf{v}_k-\mathbf{v})$ is the drift convection of phase $k$. 
The momentum conservation equation for each phase can be written as
$$\frac{d}{dt}\int_{V(t)}\tilde{\rho}_k\mathbf{v}_k\,dV + \int_{\partial V(t)} \mathbf{T}_k\!\cdot\!\mathbf{n}\,dS
= -\int_{V(t)}\phi_k\nabla p\,dV + \int_{V(t)}\mathbf{f}_k\,dV + \int_{\partial V(t)}(\phi_k\boldsymbol{\tau}_k)\!\cdot\!\mathbf{n}\,dS
$$
where $\mathbf{T}_k = \mathbf{v}_k \otimes \mathbf{Q}_k$ is the drift stress for phase $k$,  $\mathbf{f}_k$ and $\bm{\tau}_k$ is the body and shear forces acting on phase $k$.

For the mixture, the mass conservation equation can be obtained from the summation of the conservation eqation of each pahse and is
$$\frac{d}{dt}\int_{V(t)}\rho\,dV = 0 $$
which suggests that the mass within the moving volume element is invariant. Similarly, the momentum conservation equation becomes
$$\frac{d}{dt}\int_{V(t)}\rho\mathbf{v}\,dV + \int_{\partial V(t)}\mathbf{T}\!\cdot\!\mathbf{n}\,dS
= -\int_{V(t)}\nabla p\,dV + \int_{V(t)}\mathbf{f}\,dV + \int_{\partial V(t)}\boldsymbol{\tau}\!\cdot\!\mathbf{n}\,dS
$$
where $\mathbf{T} = \sum_k\mathbf{T}_k$ is the mixture drift stress, $\mathbf{f} = \sum_k \mathbf{f}_k$ is gravity and $\boldsymbol{\tau} = \sum_k \phi_k \boldsymbol{\tau}_k$ is the mixture shear stress. Note that, compared to single phase flow, the drift stress term is the extra contribution from the drift convection. Also note that, due to the canclation of mixture convection, one can also rewrite the mixture drfit stress as
$$
\mathbf{T} = \sum_k \mathbf{v}_k \otimes \mathbf{Q}_k = \tilde{\rho}_k\mathbf{u}_k \otimes \mathbf{u}_k
$$
where $\mathbf{u}_k = \mathbf{v}_k- \mathbf{v}$ is defined as phase slip velocity.  

### Boudary conditions of drift contributions

For a free surface moves with the mixture velocity, 
for each phase
$$
(\mathbf{v}_k - \mathbf{v})\cdot\mathbf{n} = 0
\quad\Longrightarrow\quad
\mathbf{Q}_k\cdot\mathbf{n} = 0,\quad
\mathbf{T}_k\cdot\mathbf{n} = 0 .
$$
The drift flux through the free surface is identically zero.
Similairly, as a solid wall is impermeable: $\mathbf{v}_k\cdot\mathbf{n} = \mathbf{v}\cdot\mathbf{n} = 0$, therefore
$$
(\mathbf{v}_k - \mathbf{v})\cdot\mathbf{n} = 0
\quad\Longrightarrow\quad
\mathbf{Q}_k\cdot\mathbf{n} = 0,\quad
\mathbf{T}_k\cdot\mathbf{n} = 0 .
$$
No drift flux crosses the wall either.

### SPH discretization of drift contributions

We consider to use a separated approach to handle the drift terms from others. Therefore, the mass conservation and the drift contribution from the momentum conservation will handeled seapartely.
The discretization of mass conservation equation for pahse $k$ on particle $i$ is gives as
$$
\frac{d}{dt}(m_{k, i}) = - \sum_j
\left(\mathbf{Q}_{k, i} + \mathbf{Q}_{k, j} \right) \nabla W_{ij}  V_iV_j$$
where $m_{k, i} = \rho^{o}_k V^{0}_{i}\phi_{k_i}$ is the mass of phase $k$. Note that, here, we assume that the change of mass for one phase in the particle is purely due to the change of volume fraction. Note that the mixture mass conservation becomes
$$\frac{d}{dt}(m_{i}) = 0$$
suggesting unchange of the mass of each particle and the net contribution from drfit convection is cancelled out. 
Similarly, the drift constribution of the momentum equation for phase $k$ is discretized as
$$
\left.\frac{d}{dt}(m_{k,i}\mathbf{v}_{k,i})\right|_{\mathrm{drift}}
= -\sum_j \bigl( \mathbf{T}_{k,i} + \mathbf{T}_{k,j} \bigr)\!\cdot\!\nabla_i W_{ij}\, V_i V_j .
$$
For the mixture moementum, the drift contibution becomes
$$\left.\frac{d}{dt}(m_{i} \mathbf{v}_{i})\right|_{\mathrm{drift}} = -\sum_j
\left(\mathbf{T}_{i} + \mathbf{T}_{j}\right) \nabla W_{ij} V_iV_j.$$
suggesting non-vanishing drfit contribution.

When the ambient phase (air or vacuum) is not modelled. Particles near the surface lack neighbours on the outside.
The missing outside particles would contribute $\mathbf{T}_{k,\mathrm{outside}}=0$, so the truncation introduces
no error.  No special boundary correction is needed; the conservative pairwise form automatically
guarantees that no momentum flows through the surface.
Similarly, if wall particles are used (for pressure and viscous terms), they are included in the drift summation with
$$
\mathbf{v}_{k,\mathrm{wall}} = \mathbf{v}_{\mathrm{wall}},\qquad
\mathbf{Q}_{k,\mathrm{wall}} = \mathbf{0},\qquad
\mathbf{T}_{k,\mathrm{wall}} = \mathbf{0}.
$$
This restores the full kernel support near the boundary without affecting conservation:
fluid–wall pairs contribute zero net momentum exchange.

### Time stepping

We first obtained the compression and mixture density. 


