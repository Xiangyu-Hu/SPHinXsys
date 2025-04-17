# Solid-to-shell coupling
The documentation describes the algorithm of solid-to-shell coupling for SPH solvers. The data transfer sequence is presented in the first section, then the mapping algorithm for unmatching interfaces are discussed in the second section.

## Dirichlet-Neumann coupling method
The Dirichilet-Neumann strategy is adopted for the coupling, where the force is mapped from origin side to the destiny side and the displacement is mapped from the destiny side to the origin side. 

Unlike FEM, the boundary of a body is not clearly defined in SPH. Hence, instead of exchanging information on the interface, we extend the solid by one layer, and transfer data on the overlapping layer of solid and shell. The solid is chosen as the destiny side and the shell is chosen as the origin side. The weight of each particle will be discussed in the next section.

![geometry](./img/geometry.png)

The geometry of the problem is different from the model, as the reference resolution of the shell can be different from the thickness, which can be a source of error.

![data_transfer](./img/data_transfer.png)

The algorithm is summarized below:

1. Execute 1st half integration of the shell body
2. Distribute the internal forces `force_` on the extented particles to the external force 'force_prior_' of the coupled solid particles
3. Execute 1st and 2nd half integration of the solid body
4. Distribute the velocity of the solid particles to the coupled shell particles
5. Execute 2nd half of the shell body

## Data mapping for unmatching interfaces
Reference: Chapter 2.2, Lindner, Florian. "Data transfer in partitioned multi-physics simulations: interpolation & communication." (2019).

A mapping is required to transfer the data between the solid and shell particles. Following the formulations of the coupling problem in FEM, the continuous fields can be approximated by the shape functions and nodal values:

$$\mathbf{u_O(x)}\approx \mathbf{N_O(x)^T} \mathbf{U_O(x)}$$
$$\mathbf{t_O(x)}\approx \mathbf{N_O(x)^T} \mathbf{T_O(x)}$$
$$\mathbf{u_D(x)}\approx \mathbf{N_D(x)^T} \mathbf{U_D(x)}$$
$$\mathbf{t_D(x)}\approx \mathbf{N_D(x)^T} \mathbf{T_D(x)}$$

where $\mathbf{u}$ and $\mathbf{U}$ represent the continuous and nodal displacement respectively, while $\mathbf{t}$ and $\mathbf{T}$ represent the continuous and nodal traction respectively. The underscript O and D denote the origin side (solid) and destiny side (shell).

For displacement, usually a direct/consistent mapping is used, where the row sum of $\mathbf{H_{OD}}$ is equal to 1.

$$\mathbf{u_O} = \mathbf{H_{OD}} \mathbf{u_D}$$

For the force, a conservative mapping is usually required. The conservation of energy at the interface is expressed as

$$\int_{\Gamma} \mathbf{u_D(x)^T} \mathbf{t_D(x)}d\Gamma = \int_{\Gamma} \mathbf{u_O(x)^T} \mathbf{t_O(x)}d\Gamma$$

, which can be reformulated as

$$\mathbf{U_D(x)^T} \int_{\Gamma} \mathbf{N_D(x)}\mathbf{N_D(x)^T} \mathbf{T_D(x)} d\Gamma = 
\mathbf{U_O(x)^T} \int_{\Gamma} \mathbf{N_O(x)}\mathbf{N_O(x)^T} \mathbf{T_O(x)} d\Gamma
$$

With the nodal force defined as

$$\mathbf{F_D(x)} = \int_{\Gamma} \mathbf{N_D(x)}\mathbf{N_D(x)^T} \mathbf{T_D(x)} d\Gamma = \mathbf{M_{DD}}\mathbf{T_D(x)}$$
$$\mathbf{F_O(x)} = \int_{\Gamma} \mathbf{N_O(x)}\mathbf{N_O(x)^T} \mathbf{T_O(x)} d\Gamma = \mathbf{M_{OO}}\mathbf{T_O(x)}$$

, we have the discrete form of energy conservation:

$$\mathbf{U_D(x)^T}\mathbf{F_D(x)} = \mathbf{U_O(x)^T}\mathbf{F_O(x)}$$

From this equation, we can directly obtain the conservative mapping matrix of the nodal force:

$$\mathbf{F_D(x)} = \mathbf{H_{DO}} \mathbf{F_O(x)} = \mathbf{H_{OD}^T} \mathbf{F_O(x)}$$

## Direct/consistent mapping with SPH approximation
For SPH method, a natural choice of consistent mapping is the SPH interpolant with Shepard filter. The displacement of a shell particle (denoted by subscription S) i can be interpolated from all the volumetric solid particles (denoted by subscription V):

$$\mathbf{U_{S,i}} = \sum_{j\in V} \tilde{\omega_{ij}} \mathbf{U_{V,j}}$$
$$\tilde{\omega_{ij}} = \frac{W_{ij}(h_{ij})V_j}{\sum_{j\in V} W_{ij}(h_{ij})V_j}$$

For a particle j outside of the neighborhood of i, the weight $\tilde{\omega_{ij}}$ is equal to 0. The (i,j) entry of the mapping matrix $\mathbf{H_{SV}}$ is the weight $\tilde{\omega_{ij}}$. Apparently, the consistency condition $\sum_j \tilde{\omega_{ij}} = 1$ is satisfied.

For the force mapping, from the discrete force of energy conservation, a conservative mapping matrix is obtained from the transpose of $\mathbf{H_{SV}}$.

$$\mathbf{F_{V,i}} = \sum_{j\in S} H_{SV,ji} \mathbf{F_{S,j}}$$

where the (j,i) entry of $\mathbf{H_{SV}}$ is:

$$H_{SV,ji} = \tilde{\omega_{ji}} = \frac{W_{ji}(h_{ji})V_i}{\sum_{k\in V} W_{jk}(h_{jk})V_k}$$

I choose $h_{ij} = max(h_i, h_j)$ as the smoothing length of two particles will different resolutions, which guarantees the symmetry $W_{ij} = W_{ji}$. Hence, the equation can also be written as:

$$H_{SV,ji} = \frac{W_{ij}(h_{ij})V_i}{\sum_{k\in V} W_{jk}(h_{jk})V_k}$$

Since $W_{ij}$ is zero outside of the neighborhood of i, we only need to consider the contribution of shell particles inside the neighborhood of solid particle i.

To implement this conservative mapping for the solid, we need to first compute and record $\sum_{k\in V} W_{jk}(h_{jk})V_k$ for each couplied shell particle j. This value should be computed from the contact relation shell-solid, instead of solid-shell. 

## Example
An example is used to prove the feasibility of this coupling method, where a solid cube is placed on a shell stripe. The cube is subjected to a constant gravity force.

