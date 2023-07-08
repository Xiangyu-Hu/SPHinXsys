# ![](../../assets/logo.png) SPHinXsys

## Conservative SPH formulation with zero- and first-order consistency

In SPH, if we assume all particles have the same mass $m$, 
the density is obtained $\rho_i = m \sum_j W_{ij}$ and 
the particle volume is obtained from $V_i = m/ \rho_i$.
The typical anti-symmetric formulation of SPH approximation of gradient is

$$
\nabla \phi_i = - \sum_j 
(\phi_i + \phi_j)  \nabla W_{ij} V_j.
$$

To ensure zero-order consistency of the approximation, 
the particle position is corrected by relaxing the zero-order 
consistency error iteratively

$$
\Delta \mathbf{r}_i = 2 \alpha h^2 \sum_j  \nabla W_{ij} V_j,
$$

where the parameter $\alpha \approx 0.2$ to ensure numerical stability, 
until the residue is sufficiently small.

To achieve first-order consistency, one introduces the correction matrix

$$
\mathbf{B}_i = \left(\sum_j 
\mathbf{r}_{ij} \otimes \nabla W_{ij} V_j \right)^{-1},
$$

and relax the particle positions similarly according to 

$$
\Delta \mathbf{r}_i = \alpha h^{2} 
\sum_j (\mathbf{B}_i + \mathbf{B}_j) \nabla W_{ij} V_j.
$$

Then, the new anti-symmetric formulation of SPH approximation of gradient is given as

$$
\nabla \phi_i = - \sum_j 
(\phi_i \mathbf{B}_j + \phi_j \mathbf{B}_i)  \nabla W_{ij} V_j.
$$

Such approximation can be rewritten as 

$$
\nabla \phi_i = 
- \phi_i \sum_j (\mathbf{B}_i + \mathbf{B}_j)  \nabla W_{ij} V_j +
\sum_j (\phi_i - \phi_j ) \mathbf{B}_i \nabla W_{ij} V_j.
$$

Note that, while the first term of the rhs vanishes due to the particle relaxation
and suggests zero-order consistency, 
the second term achieves first-order consistency 
as suggested by the property of the correction matrix.
