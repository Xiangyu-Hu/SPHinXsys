import numpy as np
from stirring_metrics import metrics,viscous_power
from normalize_source import normalize_on_particles
# Rigid rotation must not produce viscous heat; simple shear power = mu*gamma^2*volume.
g=np.array([[[0,-2,0],[2,0,0],[0,0,0]]]);assert abs(viscous_power([3],g,[.1]))<1e-12
g=np.array([[[0,2,0],[0,0,0],[0,0,0]]]);assert abs(viscous_power([3],g,[.1])-1.2)<1e-12
v=np.array([.01,.02,.03]);q=normalize_on_particles([0,1,2],v,50000);assert abs(sum(q*v)-50000)<1e-8
s=dict(rho_kg_m3=2800,impeller_diameter_m=.2,rpm=10,mu_reference_Pa_s=50,cp_J_kg_K=1500,k_W_m_K=5,P_Joule_W=56000)
s['P_viscous_W']=47*50*(10/60)**2*.2**3
r=metrics(s);assert abs(r['NpRe_using_prescribed_mu']-47)<1e-10;assert abs(r['apparent_mu_from_A47_Pa_s']-50)<1e-10
print('PASS: rigid rotation, simple shear dissipation, 3-D source normalization, SI nondimensionalization.')
