# Plate model

## Kinetics

Along with the Reissner–Mindlin type shell model, we introduce the base $\{\mathbf{e}_x,\mathbf{e}_y, \mathbf{e}_z\}$ for the global coordinate system, and the initial local coordinate $\bm{\xi} = (\xi, \eta, \zeta)$ for the reference configuration, associated with initial normal (fiber) vector $\mathbf{n}^{o} = (0, 0, 1)$.
Note that, for plate, the global coordinate system shares the same base with the initial local coordinate system

$$\mathbf{e}_x = \mathbf{e}_\xi, \quad \mathbf{e}_y = \mathbf{e}_\eta, \quad \mathbf{e}_z = \mathbf{e}_\zeta $$

Under such assumption, a initial material position $\mathbf{r}^{o}$ presented with the reference configuration can be expressed in terms of the local coordinate system as

$$\mathbf{r}^{o}(\xi,\eta,\zeta) = \mathbf{r}^{0}_m(\xi,\eta) + \zeta\,\mathbf{n}^{o}(\xi,\eta), \quad \zeta \in [-d/2, d/2]$$

where $(\xi,\eta) \equiv (\xi,\eta, 0)$ and $\mathbf{r}^{0}_m(\xi,\eta)$ is the material position vector of the mid-surface of the plate, and $d$ is the thickness of the plate.
The present material position $\mathbf{r}$ can be expressed,
with the reference configuration, as

$$
\mathbf{r}(\xi,\eta,\zeta) = \mathbf{r}_m(\xi,\eta) + \zeta\,\mathbf{n}(\xi,\eta),
\qquad \zeta \in [-d/2,\; d/2]
$$

where $\mathbf{r}_m(\xi,\eta)$ is the present material position vector of the mid-surface of the plate, and $\mathbf{n}(\xi,\eta)$ is
the pseudo-normal vector of the mid-surface of the plate,
which is obtained by the tensor $\bm{\Lambda}(\xi,\eta)\in \text{SO}(3)$ in the reference configuration from the initial normal vector

$$\mathbf{n}(\xi,\eta) = \bm{\Lambda}(\xi,\eta) \mathbf{n}^{o}(\xi,\eta)$$

and is not necessarily perpendicular to the mid-surface of the plate.

Under the total Lagrangian description, the deformation gradient $\mathbf{F}$ can be expressed as

$$\mathbf{F} = \nabla^{0} \mathbf{r} = \left(\frac{\partial\mathbf{r}}{\partial\xi},\frac{\partial\mathbf{r}}{\partial\eta}, \frac{\partial\mathbf{r}}{\partial\zeta}\right)$$

where

$$
\begin{aligned}
\frac{\partial\mathbf{r}}{\partial\xi} &= \mathbf{r}_{m,\xi} + \zeta\,\mathbf{n}_{,\xi}, \\[2mm]
\frac{\partial\mathbf{r}}{\partial\eta} &= \mathbf{r}_{m,\eta} + \zeta\,\mathbf{n}_{,\eta}, \\[2mm]
\frac{\partial\mathbf{r}}{\partial\zeta} &= \mathbf{n}(\xi,\eta).
\end{aligned}
$$

Note that, for the flat plate the reference mid-surface tangents are

$$\mathbf{r}^{0}_{m,\xi} = \mathbf{e}_x, \quad \mathbf{r}^{0}_{m,\eta} = \mathbf{e}_y$$
but the present mid-surface tangents are
$$\mathbf{r}_{m,\xi} = \mathbf{r}_{m,\xi}(\xi,\eta), \quad \mathbf{r}_{m,\eta} = \mathbf{r}_{m,\eta}(\xi,\eta)$$

From these tangents, we can define the present mid-surface normal vector $\mathbf{n}_m$ as

$$\mathbf{n}_m = \frac{\mathbf{r}_{m,\xi} \times \mathbf{r}_{m,\eta}}{\|\mathbf{r}_{m,\xi} \times \mathbf{r}_{m,\eta}\|}$$

Note that the present mid-surface normal vector $\mathbf{n}_m$ is not necessarily equal to the pseudo-normal vector $\mathbf{n}$, which is obtained from the rotation tensor $\bm{\Lambda}$.
With the present mid-surface normal and the material point at middle surface and the rotation tensor $\mathbf{Q} \in \text{SO}(3)$, one can define the present local coordinate system $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$ as

$$
\mathbf{e}_1 = \mathbf{Q} \mathbf{e}_\xi, \quad
\mathbf{e}_2 = \mathbf{Q}\,\mathbf{e}_\eta,\quad
\mathbf{e}_3 = \mathbf{Q}\,\mathbf{e}_\zeta = \mathbf{n}_m
$$

With the deformation gradient, one can define the Green-Lagrange strain tensor $\mathbf{E}$ as
$$\mathbf{E} = \frac{1}{2}(\mathbf{F}^{T}\mathbf{F} - \mathbf{I})$$

Note that, along the fiber direction,
one has the following strain component
$$E_{33} = \frac{1}{2}(\mathbf{n}\cdot\mathbf{n} - 1)$$
which is essentially zero.
The right Cauchy-Green deformation tensor $\mathbf{C}$ as

$$\mathbf{C} = \mathbf{F}^{T}\mathbf{F} = \mathbf{I} + 2\mathbf{E}$$

and the left Cauchy-Green deformation tensor $\mathbf{B}$ as

$$\mathbf{B} = \mathbf{F}\mathbf{F}^{T}$$

## Constitutive model

First, we transfer deformation gradient and the Green-Lagrangian strain tensor to the present local coordinate using the rotation tensor $\mathbf{Q}$ as

$$\bar{\mathbf{F}} = \mathbf{Q}^{T} \mathbf{F}, \quad \bar{\mathbf{E}} = \mathbf{Q}^{T} \mathbf{E} \mathbf{Q}$$

Then, we can obtain the second Piola–Kirchhoff stress tensor $\mathbf{S}$ in the present local coordinate system as
$$\bar{\mathbf{S}} = f(\bar{\mathbf{E}}), \quad \text{or} \quad \overline{\bar{S}} = f(\bar{\mathbf{C}}), \quad \bar{\mathbf{C}} = 2 \bar{\mathbf{E}} - \mathbf{I} $$

To impose the plane stress condition, we can obtain the Cauchy stress tensor $\bm{\sigma}$ in the present local coordinate system as
$$\bar{\bm{\sigma}} = \frac{1}{J} \bar{\mathbf{F}} \bar{\mathbf{S}} \bar{\mathbf{F}}^{T}, \quad J = \det(\bar{\mathbf{F}})$$

Then, we solve for the strain component $\bar{E}_{33}$ in terms of the other strain components so that the stress component $\bar{\bm{\sigma}}_{33} = 0$. This allows us to express the stress-strain relationship in a reduced form suitable for plate analysis.

Note that to satisfy another boundary condition of free in-plane shear stress at the top and bottom surfaces of the plate, one need to introduce a shear correction factor to modify the shear stress components in the constitutive model as
$$\bar{\bm{\sigma}}_{13} = \kappa \bar{\bm{\sigma}}_{13}, \quad \bar{\bm{\sigma}}_{23} = \kappa \bar{\bm{\sigma}}_{23}$$

where $\kappa$ is the shear correction factor, which can be determined based on the geometry and material properties of the plate.

After obtaining the reduced stress-strain relationship, we obtain corresponding second Piola–Kirchhoff stress tensor $\bar{\mathbf{S}}$ in the present local coordinate system as
$$\bar{\mathbf{S}} = f(\bar{\mathbf{E}})$$
and then we can transform it back to the global coordinate system using the rotation tensor $\mathbf{Q}$ as
$$\mathbf{S} = \mathbf{Q} \bar{\mathbf{S}} \mathbf{Q}^{T}$$

for the total Lagrangian description.

## Governing equations

The general governing equations in the total Lagrangian description is
$$\rho^{0}\ddot{\mathbf{r}} = \nabla^{0} \cdot \mathbf{P} + \mathbf{b}^{o}$$

where $\mathbf{P} = \mathbf{F}\mathbf{S}$ is the first Piola–Kirchhoff stress tensor and $\mathbf{b}^{o}$ is the body force per unit reference volume.
To obtain the governing equations for the plate model, we integrate the above equation through the thickness of the plate.

Define stress resultants (forces and moments per unit reference length):

$$
\begin{align*}
\mathbf{N}_\alpha(\xi,\eta) &= \int_{-d/2}^{d/2} \mathbf{P}_\alpha(\xi,\eta,\zeta)\,d\zeta,\qquad \alpha = \xi,\eta,\\[2mm]
\mathbf{Q}(\xi,\eta) &= \int_{-d/2}^{d/2} \mathbf{P}_\zeta(\xi,\eta,\zeta)\,d\zeta,\\[2mm]
\mathbf{M}_\alpha(\xi,\eta) &= \int_{-d/2}^{d/2} \zeta\,\mathbf{P}_\alpha(\xi,\eta,\zeta)\,d\zeta,\qquad \alpha = \xi,\eta.
\end{align*}
$$

Integrated body forces and inertia coefficients:

$$
\begin{align*}
\bar{\mathbf{b}}^{o} &= \int_{-d/2}^{d/2} \mathbf{b}^{o}\,d\zeta,\qquad
\bar{\mathbf{b}}^{o}_{\text{rot}} = \int_{-d/2}^{d/2} \zeta\,\mathbf{b}^{o}\,d\zeta,\\
I_0 &= \int_{-d/2}^{d/2} \rho_0\,d\zeta,\quad
I_1 = \int_{-d/2}^{d/2} \rho_0\,\zeta\,d\zeta,\quad
I_2 = \int_{-d/2}^{d/2} \rho_0\,\zeta^2\,d\zeta.
\end{align*}
$$

For a homogeneous plate with the mid‑surface chosen as the reference middle plane: $I_1=0$, $I_0 = \rho_0 d$, $I_2 = \rho_0 d^3/12$.
The integrations can be obtained numerically using Gaussian quadrature or other numerical integration methods, depending on the complexity of the stress distribution through the thickness.
Note that, in order to satisfy the plane stress condition and
the surface in-plane shear condition, one need to use the above reduced stress-strain relationship in the constitutive model to at integration points.

We obtain the translational equilibrium equation

$$
I_0\,\ddot{\mathbf{r}}_m + I_1\,\ddot{\mathbf{n}} =
\frac{\partial\mathbf{N}_\xi}{\partial\xi} + \frac{\partial\mathbf{N}_\eta}{\partial\eta}
- \mathbf{T}^{o} + \bar{\mathbf{b}}^{o},
$$

where $\mathbf{T}^{o} = \bigl[\mathbf{P}_\zeta\bigr]_{-d/2}^{d/2}$ is the total surface traction on top and bottom faces.

and Rotational (moment) equilibrium equation

$$
I_1\,\ddot{\mathbf{r}}_m + I_2\,\ddot{\mathbf{n}} =
\frac{\partial\mathbf{M}_\xi}{\partial\xi} + \frac{\partial\mathbf{M}\_\eta}{\partial\eta}
- \mathbf{Q}
+ \frac{d}{2}\bigl(\mathbf{P}_\zeta(\text{top}) + \mathbf{P}_\zeta(\text{bottom})\bigr)
+ \bar{\mathbf{b}}^{o}_{\text{rot}} .
$$

These two vector equations are the strong form governing equations for the mid‑surface of the Reissner–Mindlin plate in the total Lagrangian description.
