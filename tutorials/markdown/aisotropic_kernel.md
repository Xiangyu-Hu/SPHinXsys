# ![](../../assets/logo.png) SPHinXsys

## Anisotropic SPH (ASPH) kernel

For anisotropic kernel, one introduces a mapping between physical to normalized
position by $(\boldsymbol{r} \rightarrow \boldsymbol{\eta})$ with a linear transformation $\boldsymbol{G}$. 
To compute the kernel function in SPH, 
compared with isotropic kernel, one has the following change on the normalization as

$$ \text{SPH:} ~\boldsymbol{\eta} = \boldsymbol{r} / h \rightarrow \text{ASPH:} ~\boldsymbol{\eta} = \boldsymbol{G}\boldsymbol{r}.
$$ (1)

Under such notation, one can obtain a ASPH kernel value by $W'(\boldsymbol{G}\boldsymbol{r})$ and the kernel should still have the normalization property of 

$$
1 = \int W(\boldsymbol{\eta}) d \boldsymbol{\eta}
= \int W'(\boldsymbol{G}\boldsymbol{r}) d \boldsymbol{r} 
= \int |\boldsymbol{G}| W(\boldsymbol{\eta}) d \boldsymbol{r}. 
$$(2)

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

To obtain the formulation of Laplacian, we first consider the identity

$$
2 d = -2 \int \boldsymbol{\eta} \cdot \boldsymbol{\eta} F(\boldsymbol{\eta}) d \boldsymbol{\eta},
$$(3)

where $d$ is the dimension and 
$F(\boldsymbol{\eta}) = \frac{1}{\eta} \frac{\partial W}{\partial \eta}$.
Note that Eq. (3) can also be considered as the kernel approximation of 
the Laplacian 
$\nabla^2 (\boldsymbol{\eta} \cdot \boldsymbol{\eta}) = 2 d$.
By introducing the transformation, one has 

$$
2 d = -2 |\boldsymbol{G}| \boldsymbol{G} : \boldsymbol{G} 
\int \boldsymbol{r} \cdot \boldsymbol{r} F(\boldsymbol{\eta}) d \boldsymbol{r},
$$(4)

then the general approximation of Laplacian can be written as

$$
(\nabla^{2} \phi)_i = 2|\boldsymbol{G}| \boldsymbol{G} : \boldsymbol{G} 
\int \frac{\phi(\boldsymbol{r}_i) - \phi(\boldsymbol{r})}{\eta} 
\frac{\partial W}{\partial \eta} d \boldsymbol{r}.
$$(5)
