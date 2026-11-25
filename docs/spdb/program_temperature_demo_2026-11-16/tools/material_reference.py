"""Independent scalar reference only; not a SPHinXsys solid solver."""
import math,csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
C={'alpha_MPa_inv':[.01277,-.04834,.18608,-.36932,.37353,-.15155], 'n':[7.19785,-12.73951,57.97832,-145.77045,168.05383,-71.91728], 'Q_J_mol':[454181,-319488,-28978.6,1777650,-3257280,1745480], 'lnA':[40.18281,-29.51957,2.19563,134.56396,-247.06753,130.34468]}
def poly(c,e):
    v=0.
    for a in reversed(c):v=v*e+a
    return v
def values(e):return [poly(C[k],e) for k in C]
def rate(T,stress_MPa,e):
    if stress_MPa==0:return 0.
    a,n,q,l=values(e)
    return math.exp(l-q/(8.314*T)+n*math.log(math.sinh(a*stress_MPa)))
def stress(T,r,e):
    a,n,q,l=values(e)
    return math.asinh(math.exp((math.log(r)+q/(8.314*T)-l)/n))/a
if __name__=='__main__':
    rows=[]
    for e in [.05,.1,.5,.75]:
        sig=stress(1223.15,.01,e)
        assert abs(rate(1223.15,sig,e)/.01-1)<1e-10
        rows.append([e,*values(e),sig,rate(1223.15,5,e)])
    with (ROOT/'reference/material_checkpoints.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['strain','alpha_MPa_inv','n','Q_J_mol','lnA','flow_stress_MPa_at_rate_0.01_T1223.15','rate_at_5MPa_T1223.15']);w.writerows(rows)
    assert abs(rows[0][-2]-79.57775328757)<1e-8
    print('Scalar material roundtrip/checkpoint checks passed. No SPH simulation performed.')
