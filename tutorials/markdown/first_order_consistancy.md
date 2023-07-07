# ![](../../assets/logo.png) SPHinXsys

## Conservative zero- and first-order consistency corrections

In SPH, if we assume all particles have the same mass $m$, 
the density is obtained $\rho_i = m \sum_j W_{ij}$ and 
the particle volume is obtained from $V_i = m/ \rho_i$.
The typical anti-symmetric formulation of SPH approximation of gradient is

$$
\nabla \phi_i = - \sum_j 
(\phi_i + \phi_j)  \nabla W_{ij} V_j.
$$

To ensure the zero-order consistency of the approximation, 
the particle position is corrected by relaxing the zero-order 
consistency error iteratively

$$
\Delta \boldsymbol{r}_i = 2 \alpha h^{2} \sum_j  \nabla W_{ij} V_j,
$$

where the parameter $\alpha \approx 0.2$ to ensure numerical stability.

$$
\boldsymbol{B}_i = \left(\sum_j 
\boldsymbol{r}_{ij} \otimes \nabla W_{ij} V_j \right)^{-1},
$$