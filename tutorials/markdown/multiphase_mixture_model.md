## Governing equations in the Lagrangian frame moving with average velocity

In this work, we consider a weakly compressible multiphase flow model in the average velocity frame. The governing equations consist of the mass conservation and momentum conservation equations for each phase, as well as the equation of state for the mixture.

### Primary and derived variables
We consider two different set of variables, one is for the mixture, the other is for each phase. 

For the mixture, we use compression ratio $\beta$ as a primary variable, which is defined as the ratio of the initial volume element to the current volume element, that is $\beta = \frac{V_0}{V}$, where $V_0$ is the initial volume and $V$ is the current volume. 
For each phase $k$, we use the phase volume fraction $\phi_k$ and the phase velocity $\mathbf{v}_k$ as the primary variables.

The mixture density $\rho$ and the mixture velocity $\mathbf{v}$ can be computed from these primary variables as follows:
$$\rho = \beta\sum_{k} \phi_k \rho^{o}_k$$
and the mass-averaged velocity is defined as
$$\mathbf{v} = \frac{\sum_{k} \phi_k \rho^{o}_k \mathbf{v}_k}{\sum_{k} \phi_k \rho^{o}_k}$$
where $\rho^{o}_k$ is the reference density of phase $k$. If we assume that all phases share the same speed of sound, then the pressure can be computed from the mixture density using the equation of state for the mixture, which is given by
$$p =  \rho_0 c_0^2 (\beta - 1)$$
where $c_0$ is the reference sound speed and $\rho_0$ is the reference density of the mixture, which can be computed as $\rho_0 = \sum_{k} \phi_k \rho^{o}_k$. Furthermore, we assume all phase share the same pressure, which is a common assumption for multiphase flow models.

For each phase, the density $\rho_k$ can be computed from the phase volume fraction and the reference density as $\rho_k = \beta\phi_k \rho^{o}_k$. The momentum of each phase can be computed as $\mathbf{m}_k = \rho_k \mathbf{v}_k = \beta\phi_k \rho \mathbf{v}_k$. Note that one can obtain the mixture velocity from the phase velocity as $\mathbf{v} = \frac{1}{\rho} \sum_{k} \mathbf{m}_k$, which is consistent with the definition of the mass-averaged velocity.

### Governing equations
Base on the above definitions, the governing equations of a Lagrangian volume element for the weakly compressible multiphase flow model,
two set of equations can be derived, one is for the mixture, the other is for each phase. 

For the mixture, one has the  the evolution of the compression ratio $\beta$, which is given as
$$\frac{d\beta}{dt} = -\beta \nabla \cdot \mathbf{v}$$
where the material derivative is defined as $\frac{d}{dt} = \frac{\partial}{\partial t} + \mathbf{v} \cdot \nabla$. 
For each phase, one first has the mass conservation equation, which is given as
$$\frac{d\phi_k}{dt} = - \nabla \cdot \left[\phi_k  ( \mathbf{v}_k - \mathbf{v})\right]$$
which describes the drift of each phase relative to the mixture. The momentum conservation equation for each phase can be written as
$$\frac{d\mathbf{m}_k}{dt} = - \nabla \cdot \left[\mathbf{m}_k \otimes (\mathbf{v}_k - \mathbf{v})\right] - \phi_k \nabla p + \mathbf{f}_k + \nabla \cdot (\phi_k \bm{\tau}_k) $$
where $\mathbf{f}_k$ and $\bm{\tau}_k$ is the body and shear forces acting on phase $k$.
