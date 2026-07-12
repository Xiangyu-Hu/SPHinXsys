# Lagrangian incompressible multiphase flow

In this work, we consider an incompressible multiphase flow
using weakly compressible fluid model in average velocity frames.
The governing equations consist of the mass conservation
and momentum conservation equations for each phase,
as well as the equation of state for the mixture.
We first discuss the case with mass-averaged frame velocity
and then extend the formulations with volume-average
and more complex gradient-based average velocities.

## Primary and derived variables

We consider two sets of variables, one is for the mixture, the other is for each phase.
For the mixture, we use compression ratio $\beta$ as a primary variable,
which is defined as the ratio of the initial volume element
to the current volume element, that is $\beta = \frac{V^{o}}{V}$,
where $V^{o}$ is the initial volume and $V$ is the current volume.
For each phase $k$, we use the phase volume fraction $\phi_k$,
with $\sum_k \phi_k = 1$, and the phase velocity $\mathbf{v}_k$ as the primary variables.

The mixture density $\rho$ can be computed from these primary variables as
$$\rho = \beta\sum_{k} \phi_k \rho^{o}_k$$
where $\tilde{\rho}_k = \beta\phi_k \rho^{o}_k$ and $\rho^{o}_k$ are
the partial density and the reference density, respectively, of phase $k$.
Also $\rho^{o} = \sum_{k} \phi^{o}_k \rho^{o}_k$ is the reference density of the mixture,
which is given by the initial mixture density with
initial volume fractions $\phi^{o}_k$ in the initial volume element.
The mass-averaged velocity $\mathbf{v}^{m}$ is defined as
$$\mathbf{v}^{m} = \frac{\sum_{k} \phi_k \rho^{o}_k \mathbf{v}_k}{\sum_{k} \phi_k \rho^{o}_k}$$
If we assume that all phases share the same speed of sound,
then the pressure can be computed from the mixture density
using the equation of state for the mixture, which is given by
$$p =  \rho^{o} c^2 (\beta - 1)$$
where $c$ is the artificial (reference) speed of sound.
Furthermore, we assume all phase share the same pressure,
which is a common assumption for multiphase flow models.
Note that, under the weakly compressible fluid model,
$c$ is chosen so that the density variation is less or about 1\%.
We assume that the initial mass of the volume element $m$
as $$m = \rho^{o} V^{o}$$ which is dependent of reference density
and initial element volume only.
As will be shown latter the mass of the fluid element
is an invariant under the mass-averaged velocity moving frame.

For each phase, the partial density $\tilde{\rho}_k$ can be computed
from the phase volume fraction and the reference density
as $\tilde{\rho}_k = \beta\phi_k \rho^{o}_k$.
Note that one can obtain the mixture velocity from
the phase velocity as $\mathbf{v}^{m} = \frac{1}{\rho} \sum_{k} \beta\phi_k \rho$,
which is consistent with the definition of the mass-averaged velocity
due to the cancellation of $\beta$.

## Governing equations

For the mixture, one has the evolution of the compression ratio $\beta$, which is given as
$$\frac{d\beta}{dt} = -\beta \nabla \cdot \mathbf{v}^{m}$$
where the material derivative is defined as 
$\frac{d}{dt} = \frac{\partial}{\partial t} + \mathbf{v}^{m} \cdot \nabla$.
Note the evolution of compression ratio is purely kinematic and geometric and no physics involved.

There are two sets of conservation equations of a Lagrangian volume element
for the weakly compressible multiphase flow model:
one is for each phase, the other is for the mixture.

For each phase, first is the mass conservation equation, which is given as
$$\frac{d}{dt}\int_{V(t)}\tilde{\rho}_k\,dV + \int_{\partial V(t)}\mathbf{Q}_k\!\cdot\!\mathbf{n}\,dS = 0 $$
which describes the drift of each phase relative to the mixture.
Here, $\mathbf{Q}^{m}_k = \tilde{\rho}_k(\mathbf{v}_k-\mathbf{v}^{m})$ is the drift convection of phase $k$.
The momentum conservation equation for each phase can be written as

$$
\frac{d}{dt}\int_{V(t)}\tilde{\rho}_k\mathbf{v}_k\,dV + \int_{\partial V(t)} \mathbf{T}^{m}_k\!\cdot\!\mathbf{n}\,dS
= -\int_{V(t)}\phi_k\nabla p\,dV + \int_{V(t)}\mathbf{f}_k\,dV + \int_{\partial V(t)}(\phi_k\boldsymbol{\tau}_k)\!\cdot\!\mathbf{n}\,dS
$$

where $\mathbf{T}^{m}_k = \mathbf{v}_k \otimes \mathbf{Q}^{m}_k$ is the drift stress for phase $k$, $\mathbf{f}_k$ and $\bm{\tau}_k$ is the body and shear forces acting on phase $k$.

For the mixture, the mass conservation equation can be obtained from the summation of the conservation equation of each phase and is
$$\frac{d}{dt}\int_{V(t)}\rho\,dV = 0 $$
which suggests that the mass within the moving volume element is invariant. Similarly, the momentum conservation equation becomes

$$
\frac{d}{dt}\int_{V(t)}\rho\mathbf{v}^{m}\,dV + \int_{\partial V(t)}\mathbf{T}^{m}\!\cdot\!\mathbf{n}\,dS
= -\int_{V(t)}\nabla p\,dV + \int_{V(t)}\mathbf{f}\,dV + \int_{\partial V(t)}\boldsymbol{\tau}\!\cdot\!\mathbf{n}\,dS
$$

where $\mathbf{T}^{m} = \sum_k\mathbf{T}^{m}_k$ is the mixture drift stress, $\mathbf{f} = \sum_k \mathbf{f}_k$ is gravity and $\boldsymbol{\tau} = \sum_k \phi_k \boldsymbol{\tau}_k$ is the mixture shear stress. Note that, compared to single phase flow, the drift stress term is the extra contribution from the drift convection. Also note that, due to the canceling of mixture convection, one can also rewrite the mixture drift stress as

$$
\mathbf{T}^{m} = \sum_k \mathbf{v}_k \otimes \mathbf{Q}^{m}_k = \tilde{\rho}_k\mathbf{u}^{m}_k \otimes \mathbf{u}^{m}_k
$$

where $\mathbf{u}^{m}_k = \mathbf{v}_k- \mathbf{v}^{m}$ is defined as mass-averaged phase slip velocity.
Although the mixture drift stress is symmetric and depends only on the slip velocities,
the individual phase drift stress cannot be expressed solely in terms
of the slip velocity because the additional term $\mathbf{v}^{m} \otimes \mathbf{Q}^{m}_k$
cancels only after summation over all phases.

## Boundary conditions of drift contributions

For a free surface moves with the mixture velocity,
for each phase

$$
(\mathbf{v}_k - \mathbf{v})\cdot\mathbf{n} = 0
\quad\Longrightarrow\quad
\mathbf{Q}_k\cdot\mathbf{n} = 0,\quad
\mathbf{T}_k\cdot\mathbf{n} = 0 .
$$

The drift flux through the free surface is identically zero.
Similarly, as a solid wall is impermeable: $\mathbf{v}_k\cdot\mathbf{n} = \mathbf{v}^{m}\cdot\mathbf{n} = 0$, therefore

$$
(\mathbf{v}_k - \mathbf{v})\cdot\mathbf{n} = 0
\quad\Longrightarrow\quad
\mathbf{Q}_k\cdot\mathbf{n} = 0,\quad
\mathbf{T}_k\cdot\mathbf{n} = 0 .
$$

No drift flux crosses the wall either.
Note that, we omit the superscript $^{m}$ as these boundary conditions are still valid
when other average moving frame is chosen. 

## SPH discretization of drift contributions

We consider to use a separated approach to handle the drift terms from others. Therefore, the mass conservation and the drift contribution from the momentum conservation will handled separately.
The discretization of mass conservation equation for phase $k$ on particle $i$ is gives as

$$
\frac{d}{dt}(m_{k, i}) = - \sum_j
\left(\mathbf{Q}^{m}_{k, i} + \mathbf{Q}^{m}_{k, j} \right) \cdot \nabla W_{ij}  V_iV_j
$$

where $m_{k, i} = \rho^{o}_k V^{0}_{i}\phi_{k_i}$ is the mass of phase $k$. Note that, here, we assume that the change of mass for one phase in the particle is purely due to the change of volume fraction. Note that the mixture mass conservation becomes
$$\frac{d m_{i}}{dt} = 0$$
suggesting no change of the mass of each particle and the net contribution from the drift convection is cancelled out.
Similarly, the drift contribution of the momentum equation for phase $k$ is discretized as

$$
\left.\frac{d}{dt}(\mathbf{m}_{k,i})\right|_{\mathrm{drift}}
= -\sum_j \bigl( \mathbf{T}^{m}_{k,i} + \mathbf{T}^{m}_{k,j} \bigr)\!\cdot\!\nabla_i W_{ij}\, V_i V_j
$$
where $\mathbf{m}_{k,i} = m_{k,i}\mathbf{v}_{k,i}$ is the momentum of phase $k$ on particle $i$.
For the mixture momentum, the drift contribution becomes

$$
\left.\frac{d}{dt}(\mathbf{m}^{m}_{i})\right|_{\mathrm{drift}} = -\sum_j
\left(\mathbf{T}^{m}_{i} + \mathbf{T}^{m}_{j}\right) \cdot \nabla W_{ij} V_iV_j.
$$
where $\mathbf{m}^{m}_{i} = m_{i} \mathbf{v}^{m}_{i}$ is the mixture momentum on particle $i$,
suggesting non-vanishing drift contribution.

When the ambient phase (air or vacuum) is not modelled. Particles near the surface lack neighbors on the outside.
The missing outside particles would contribute $\mathbf{T}_{k,\mathrm{outside}}=0$, so the truncation introduces
no error. No special boundary correction is needed; the conservative pairwise form automatically
guarantees that no momentum flows through the surface.
Similarly, if wall particles are used (for pressure and viscous terms), they are included in the drift summation with

$$
\mathbf{v}_{k,\mathrm{wall}} = \mathbf{v}_{\mathrm{wall}},\qquad
\mathbf{Q}_{k,\mathrm{wall}} = \mathbf{0},\qquad
\mathbf{T}_{k,\mathrm{wall}} = \mathbf{0}.
$$

This restores the full kernel support near the boundary without affecting conservation:
fluid–wall pairs contribute zero net momentum exchange.

## Time stepping

We first obtained the compression ratio and mixture density.
Then, we update the phase mass for each phase from the mass conservation equations.
After this, we need a regularization approach so that
the volume fractions ensuring non-negative and unit sum of the volume fractions,
under the condition of invariant mixture mass for each particle.
To this end, we first clip the negative phase mass so that the phase mass $\tilde{m}_k$ is non-negative.
Then we assume the regularization factor $\delta_k = a \tilde{m}_k + b$, so the
regularized phase mass is given as $m_k = \tilde{m}_k (1 + \delta_k)$
where $a$ and $b$ are the coefficients to be determined.
Note that, the regularization factor ensures larger phase mass will have larger modification,
which is consistent with the physical meaning of the regularization.
The coefficients $a$ and $b$ are determined by the following two conditions:

$$
(1+ b)\sum_k \tilde{m}_k  + a \sum_k \tilde{m}^2_k  = m
$$

and

$$
(1+ b)\sum_k \frac{\tilde{m}_k}{\rho^{o}_k}  + a \sum_k\frac{\tilde{m}^{2}_k}{\rho^{o}_k}  = V^{o}
$$

After obtain the coefficients $a$ and $b$ by solving the above two equations,
we can update the phase volume fraction for each phase as $\phi_k = \frac{m_k}{\rho^{o}_k V^{o}}$.

Then, the momentum of each phase is updated from the momentum conservation equations.
Note that, if the phase mass is clipped to zero, the corresponding momentum is also set to zero,
and the clipped momentum $\mathbf{m}^{c}$ will be redistributed to other phases according the phase mass,
that is, the incremental momentum of phase $k$ is
$$\Delta_k = \frac{m_k \mathbf{m}^{c}}{\sum_{k} m_k} $$

Finally, we update the mixture and phase velocities as
$$  \mathbf{v}_i = \frac{\mathbf{m}^{m}_{i}}{m_{i}}, \quad 
\mathbf{v}_{k,i} = \frac{\mathbf{m}_{k,i} +\epsilon_1 \mathbf{m}^{m}_{i}}{m_{k,i} +\epsilon_1 m_{i}} $$
where $\epsilon_1$ is a small number to avoid division by zero 
when the phase mass is clipped to zero.

## Governing equations in the Lagrangian frame moving with volume-averages velocity

In this case, the mixture velocity defined as
the volume-averaged velocity $\mathbf{v}^{\phi}$ is defined as
$$\mathbf{v}^{\phi} = \sum_{k} \phi_k \mathbf{v}_k$$
which leads to no formal changes of the phase equations except now
the drift convection is defined as $\mathbf{Q}^{\phi}_k =  \tilde{\rho}_k(\mathbf{v}_k-\mathbf{v}^{\phi})$.
The most significant change is the mixture conservation equation now is

$$
\frac{d}{dt}\int_{V(t)} \rho \, dV
+ \int_{\partial V(t)} \mathbf{Q}^{\phi}\!\cdot\!\mathbf{n}\, dS = 0
$$

where

$$
\mathbf{Q}^{\phi} \equiv \sum_k \mathbf{Q}^{\phi}_k
= \sum_k \tilde{\rho}_k (\mathbf{v}_k - \mathbf{v}^{\phi}) \neq 0
$$

Note that one still can evaluate the fluid density in the same formulation 
used for mass-average frame velocity but the reference density now
$\rho^{o} = \sum_{k} \phi_k \rho^{o}_k$
is dependent of the present volume fractions,
by assuming volume fraction independent of compression,
other than the initial value.
Another slightly change is the form of the mixture drift stress

$$
\mathbf{T}^{\phi} = \sum_k \mathbf{v}_k \otimes \mathbf{Q}^{\phi}_k
$$

which is not a symmetric form anymore.
Note that, for SPH discretization of drift contributions,
the only considerable change is now the mixture mass evolution equation

$$
\frac{dm_{i}}{dt} = - \sum_j
\sum_k\left(\mathbf{Q}^{\phi}_{k, i} + \mathbf{Q}^{\phi}_{k, j} \right) \cdot \nabla W_{ij}  V_iV_j
$$

Similarly, the time stepping almost has no change except
the mixture density is updated after that of particle mass,
which later is used in the first condition of regularization.

## A Lagrangian frame moving with complex-average velocity

We first compute the gradient of each phase $\nabla \phi_k$.
Our assumption is that, when the gradient is significant,
the velocity along it will be mass averaged. Otherwise,
the flow interface is considered as smeared,
or single phase is assumed,
hence volume average will be used.
Therefore, one can define a complex average velocity $\mathbf{v}^{c}$ blended by

$$
\mathbf{v}^{c} = \mathbf{v}^{\phi}
+ \sum_k \phi_k\,
\frac{\bigl(\mathbf{v}^{m} - \mathbf{v}^{\phi}\bigr)\!\cdot\!\nabla\phi_k}
{|\nabla\phi_k|^2 + \epsilon_2/h^2}\,
\nabla\phi_k
$$

Note that, with the new frame velocity, the governing equation and SPH discretization
are formally the same as those using volume average.
Here, the small number can be choosing as $\epsilon_2 \ll h^2$,
where $h$ is smoothing length of SPH smoothing kernel.

## Consistency of these models for inviscid, immiscible and initially pure phase flows

For a flow initially composed of pure, immiscible inviscid fluids,
the drift terms vanish exactly at the particle level,
the phase masses remain constant, and the numerical method reproduces pure fluid dynamics with a sharp interface. 
The model is entirely consistent and does not suffer from artificial mixing in this scenario.