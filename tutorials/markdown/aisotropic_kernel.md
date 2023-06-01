# ![](../../assets/logo.png) SPHinXsys

## Anisotropic SPH (ASPH) kernel

For anisotropic kernel, one introduces a mapping between physical to normalized
position by $(\boldsymbol{r} \rightarrow \boldsymbol{\eta})$ with a linear transformation $\boldsymbol{G}$. 
To compute the kernel function in SPH, 
compared with isotropic kernel, one has the following change on the normalization as

$$ \text{SPH:} ~\boldsymbol{\eta} = \boldsymbol{r} / h \rightarrow \text{ASPH:} ~\boldsymbol{\eta} = \boldsymbol{G}\boldsymbol{r}.
$$

Under such notation, one can obtain a ASPH kernel value by $W'(\boldsymbol{G}\boldsymbol{r})$ and the kernel should still have the normalization property of 

$$
1 = \int W(\boldsymbol{\eta}) d \boldsymbol{\eta}
= \int W'(\boldsymbol{G}\boldsymbol{r}) d \boldsymbol{r} 
= \int |\boldsymbol{G}| W(\boldsymbol{\eta}) d \boldsymbol{r}. 
$$

Here, $|\boldsymbol{G}| dr = d \boldsymbol{\eta}$, 
where $|\boldsymbol{G}|$ is the determinate, 
and $W(\boldsymbol{\eta}) \equiv W(\boldsymbol{G}\boldsymbol{r})$ is the value obtained
by the original isotropic kernel.
Based on this, one has the approximation of derivative of a field by

$$
(\nabla \phi)_i = \int \left(\phi(\boldsymbol{r}_i) - \phi(\boldsymbol{r})\right) \nabla W’(\boldsymbol{G}(\boldsymbol{r}_i - \boldsymbol{r})) 
d \boldsymbol{r}
=  |\boldsymbol{G}| \boldsymbol{G} 
\int \left(\phi(\boldsymbol{r}_i) - \phi(\boldsymbol{r})\right)
\nabla W(\boldsymbol{\eta}) 
d \boldsymbol{r}.
$$

Another way to obtain the approximation of derivative is using
the following identity

$$
\boldsymbol{I} = - \int \boldsymbol{\eta} \otimes \nabla W(\boldsymbol{\eta})  d \boldsymbol{\eta},
$$

where $\boldsymbol{I}$ is the identity matrix. 
With the transformation, one has

$$
\boldsymbol{I} = - |\boldsymbol{G}| \boldsymbol{G} \int \boldsymbol{r} \otimes \nabla W(\boldsymbol{\eta})  d \boldsymbol{r}.
$$

Note that the above equation can also be considered as the kernel approximation of 
the gradient  
$\nabla \boldsymbol{\eta} = \boldsymbol{I}$.
Therefore, one can substitute $\boldsymbol{r}$ in the above equation with 
the difference of the field to position $i$, that is 
$\boldsymbol{r} \rightarrow \phi(\boldsymbol{r}) - \phi(\boldsymbol{r}_i)$,
to obtain the same approximation of derivate.

To obtain the formulation of Laplacian, we first consider the relation

$$
\Delta \phi(\boldsymbol{r}) = 
\text{Tr}[  \nabla \boldsymbol{D} \nabla \phi(\boldsymbol{\eta})],
$$

where $\boldsymbol{D} = \boldsymbol{G} \boldsymbol{G}^T$.
From Espanol and Revenga (2003), one has 

$$
\nabla \nabla \phi(\boldsymbol{\eta}) =
\int \frac{\phi(\boldsymbol{0}) - \phi(\boldsymbol{\eta})}{\eta}
\left[\frac{\boldsymbol{\eta}}{\eta} 
\otimes \nabla W(\boldsymbol{\eta}) (d + 2)  
- \boldsymbol{I}\frac{\boldsymbol{\eta}}{\eta} \cdot \nabla W(\boldsymbol{\eta}) \right] 
d \boldsymbol{\eta}.
$$

The final formulation of the Laplacian is 

$$
(\Delta \phi)_i = |\boldsymbol{G}|  
\int \frac{\phi(\boldsymbol{r}_i) - \phi(\boldsymbol{r})}{\eta}
\left[ \boldsymbol{D}: 
\frac{\boldsymbol{\eta}}{\eta}\frac{\boldsymbol{\eta}}{\eta}
(d + 2)  - \text{Tr}(\boldsymbol{D}) 
\frac{\boldsymbol{\eta}}{\eta} \cdot \frac{\boldsymbol{\eta}}{\eta} \right] 
\frac{\partial W}{\partial \eta}  d \boldsymbol{r}.
$$

The discretized form is

$$
(\Delta \phi)_i = |\boldsymbol{G}|  
\sum_j \frac{\phi_{ij}}{\eta_{ij}}
\left[ \boldsymbol{D}: 
\boldsymbol{e}_{ij} \boldsymbol{e}_{ij}
(d + 2)  - \text{Tr}(\boldsymbol{D}) 
\boldsymbol{e}_{ij} \cdot \boldsymbol{e}_{ij} \right] 
\frac{\partial W_{ij}}{\partial \eta_{ij}} V_j,
$$

where $\boldsymbol{e}_{ij} = \boldsymbol{\eta}_{ij}/\eta_{ij}$.
To increase numerical accuracy, one can introduce the correction matrix

$$
\boldsymbol{B}_i = - \left(|\boldsymbol{G}|\sum_j 
\boldsymbol{\eta}_{ij} \boldsymbol{e}_{ij} 
\frac{\partial W_{ij}}{\partial \eta_{ij}} V_j \right)^{-1},
$$

so that 

$$
(\Delta \phi)_i = |\boldsymbol{G}|  
\sum_j \frac{\phi_{ij}}{\eta_{ij}}
\left[ \boldsymbol{D}: 
\boldsymbol{e}_{ij} \overline{\boldsymbol{B}}_{ij} \boldsymbol{e}_{ij}
(d + 2)  - \text{Tr}(\boldsymbol{D}) 
\boldsymbol{e}_{ij} \cdot \overline{\boldsymbol{B}}_{ij} \boldsymbol{e}_{ij} \right] 
\frac{\partial W_{ij}}{\partial \eta_{ij}} V_j,
$$

where $\overline{\boldsymbol{B}}_{ij} = (\boldsymbol{B}_{i} + \boldsymbol{B}_{j})/2$ is the inter-particle average of the correction matrix. 

