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

To obtain the formulation of Laplacian, similarly, 
we first consider the identity

$$
2 d = -2 \int \boldsymbol{\eta}^2 F(\boldsymbol{\eta}) d \boldsymbol{\eta},
$$

where $d$ is the dimension and 
$F(\boldsymbol{\eta}) = \frac{1}{\eta} \frac{\partial W}{\partial \eta}$.
Again, the above equation can also be considered as 
the kernel approximation of the Laplacian 
$\Delta \boldsymbol{\eta}^2 = 2 d$.
By introducing the transformation, one has 

$$
2 d = -2 |\boldsymbol{G}| \boldsymbol{G} : \boldsymbol{G} 
\int \boldsymbol{r}^2 F(\boldsymbol{\eta}) d \boldsymbol{r},
$$

then, with the similar substitution as for derivative
$\boldsymbol{r}^2 \rightarrow \phi(\boldsymbol{r}) - \phi(\boldsymbol{r}_i)$, 
the general approximation of Laplacian can be written as

$$
(\Delta \phi)_i = 2|\boldsymbol{G}| \boldsymbol{G} : \boldsymbol{G} 
\int \frac{\phi(\boldsymbol{r}_i) - \phi(\boldsymbol{r})}{\eta} 
\frac{\partial W}{\partial \eta} d \boldsymbol{r}.
$$

