# Algorithm for one-sided FSI problem
The basic idea is to retrieve kernel completeness by using dummy shell contact particles. The algorithm is summarized below:

# Ghost particle integration
For a shell particle with index j in the neighborhood of a fluid particle i, we integrate the force dW_ijV_je_ij by considering ghost particles.

First we check if the shell particle itself is in the support of fluid particle i. If $r_{ij}^0$ is smaller than the cut-off radius of fluid, we calculate the force $\nabla W_{ij}^0 A_j^0=dW_{ij}^0 \mathbf{e}_{ij}^0 A_j^0$. Here $A_j^0$ is the area of shell particle.

Then we calculate the contribution of the k-th ghost particle. The distance between ghost particles is the reference spacing of shell $dp_s$. The position of the k-th particle is then:

$$\mathbf{r}_{ij}^k=\mathbf{r}_{ij}^0 + k\cdot dp_s \cdot \mathbf{n}_j=\mathbf{r}_{ij}^{k-1}+dp_s \cdot \mathbf{n}$$

where $\mathbf{n}_j$ is the corrected normal direction of shell particle j, which points from fluid to shell, i.e. opposite to $\mathbf{e}_{ij}$. If $e_{ij} \cdot n_j=-1$, $\mathbf{n}_{corrected} = \mathbf{n}_j$. If $e_{ij} \cdot n_j=1 \text{ or } 0$, $\mathbf{n}_{corrected} = -\mathbf{n}_j$.

If the distance $r_{ij}^k=\|\mathbf{r}_{ij}^k\|$ is less than the cut-off radius of fluid particle, we will sum up the force provided by this ghost particle $\nabla W_{ij}^nA_j$. The process is repeated until $r_{ij}^N$ becomes greater than the cut-off radius. Then the equavilent force can be calculated as:

$$\overline{dW_{ij}A_j}=\sum_{k=0}^N dW_{ij}^k A_j^k$$
$$\overline{\mathbf{e}_{ij}}=\frac{\sum_{k=0}^N dW_{ij}^k A_j^k\mathbf{e}_{ij}^k}{\sum_{k=0}^N dW_{ij}^k A_j^k}$$

To transform the area to volume, we have:
$$\overline{dW_{ij}V_j}=\overline{dW_{ij}A_j} \cdot dp_s$$

Similarly, we can integrate the weight $W_ij$:

$$\overline{W_{ij}}=\frac{1}{A_j^0}\sum_{k=0}^N W_{ij}^k A_j^k \cdot dp_s$$

This term is also multiplied by $dp_s$, since in `DensitySummationComplex`, the mass of contact shell `this->contact_particles_[k]->mass_` is actually the area mass. To avoid changing of codes, we turn the weight to a volumetric one here.

Unlike the force, which is integrated in the kernel support, the density ghost particles are summed up within shell thickness, i.e. until $(k+1)\cdot dp_s<t$, where $t$ is the thickness of shell. (Not modified yet)

The method is illustrated in the graph below:

![](./ghost_particles.png)

The calculation of ghost particle area $A_j^k$ will be explained in details in the next sections.

# Shell curvature
We can compute the total curvature H and Gaussian curvature K by:

$$H=\kappa_1 + \kappa_2 = \nabla \cdot \mathbf{n}_r$$

$$K=\kappa_1 \cdot \kappa_2 = \frac{1}{2}((\nabla \cdot \mathbf{n}_r)^2-\sum_i \sum_j (\frac{\partial n_j}{\partial x_i})^2)$$

where $\mathbf{n}_r$ is the normal direction of the surface, $\kappa_1$ and $\kappa_2$ are the two principle curvatures. 

For the initial configuration, the gradient of normal can be calculated as:

$$\nabla^0{\mathbf{n}^0}=-\sum_j{(\mathbf{n}^0_i-\mathbf{n}^0_j) \times \nabla^0W_{ij}V_j^0}$$

Here $\nabla^0$ is the gradient in the initial global coordinate.

Note that the bending deformation tensor in the current global configuration is defined as:

$$\mathbb{F}_2=-\sum_j[(\mathbf{n}_i-\mathbf{n}_i^0)-(\mathbf{n}_j-\mathbf{n}_j^0)] \times \nabla^0W_{ij}V_j^0
=\nabla^0{\mathbf{n}}-\nabla^0{\mathbf{n}^0}$$

The $\mathbf{n}$ used in the calculation of bending deformation is actually the pseudo normal, but here we assume that the changing of real and pseudo normal is very close to each other.

Thus $\nabla^0{\mathbf{n}}$ can be computed from the local bending deformation gradient $\mathbb{F}_2^L$:

$$\nabla^0{\mathbf{n}}=\nabla^0{\mathbf{n}^0}+(\mathbb{Q}^0)^T \cdot \mathbb{F}_2^L \cdot \mathbb{Q}^0$$

We can reuse $\mathbb{F}_2^L$ (`F_bending_[index_i]`).

Calculating in the initial global configuration instead of the current one will definitely bring error, but since the curvature is only used for approximate volume change, we will use this formulation at the moment.

The gradient in the current configuration can be computed by:

$$\nabla{\mathbf{n}}=(F_1^{-T}\nabla^0){\mathbf{n}}=-\sum_j(\mathbf{n}_i-\mathbf{n}_j)
 \times F_m^{-T} \cdot \mathbf{e}_{ij}^0dW_{ij}^0V_j^0=
-\sum_j(\mathbf{n}_i-\mathbf{n}_j)
 \cdot [dW_{ij}^0V_j^0(\mathbf{e}_{ij}^0)^T \cdot F_m^{-1}] $$
 
$$F_1^{-1}=(\mathbb{Q}^0)^T \cdot (\mathbb{F}_1^L)^{-1} \cdot \mathbb{Q}^0$$

where $\mathbb{F}_1^L$ is `F_[index_i]`.

# Update volume change
![](./volume_change.png)

The change in the area of so-called dummy particles are updated as a cone. For the 2D case, the length of the particle along the normal direction can be written as:

$$l^k=\alpha (R+k\cdot dp)$$

where superscript $n$ refers to the nth dummy particle, $dp$ is the particle reference resolution, curve radius $R=1/k$ and $\alpha$ is the angle of the cone.

The length of the real shell particle is:

$$l^0=\alpha R$$

Thus, the length of the k-th dummy particle takes the form:

$$l^k=l^0(1+k\cdot dp/R)=l^0(1+k\cdot dp \cdot \kappa)$$

where $\kappa$ is the curvature of the curve.

For 3D problems, we can calculate the length in two principle directions, but since this is only a rough estimation, we can estimate the length change by the mean curvature $H/2$:

$${(l^k)}^2 =(l^0)^2(1+k\cdot dp \cdot H/2)^2$$

The area of both 2D and 3D shells can be written as:

$$A^k =[A^0(1+k\cdot dp \cdot \frac{H}{D-1})]^{D-1}$$

# Reference
- Nitschke, Ingo, Axel Voigt, and Jörg Wensch. "A finite element approach to incompressible two-phase flow on manifolds." Journal of Fluid Mechanics 708 (2012): 418-438.