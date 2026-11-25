"""SI inputs; mechanical dissipation power is separate from Joule/heat-transfer power."""
import json,math,argparse

def viscous_power(mu, velocity_gradient, volume):
    import numpy as np
    g=np.asarray(velocity_gradient); D=.5*(g+g.swapaxes(-1,-2))
    return float(np.sum(np.asarray(volume)*np.asarray(mu)*2*np.sum(D*D,axis=(-2,-1))))

def metrics(s):
    rho,d,N,mu,cp,k=[float(s[x]) for x in ['rho_kg_m3','impeller_diameter_m','rpm','mu_reference_Pa_s','cp_J_kg_K','k_W_m_K']]
    N/=60
    if min(rho,d,N,mu,cp,k)<=0:raise ValueError('Positive SI inputs required; N is revolutions per second after conversion.')
    P=float(s['P_viscous_W']);J=float(s['P_Joule_W'])
    if P<=0 or J<=0:raise ValueError('Positive mechanical and Joule powers required.')
    npow=P/(rho*N**3*d**5);re=rho*N*d*d/mu;pr=mu*cp/k
    ree=47/npow;mue=rho*N*d*d/ree;pre=mue*cp/k
    out=dict(N_rev_s=N,Re_using_prescribed_mu=re,Pr_using_prescribed_mu=pr,Np=npow,NpRe_using_prescribed_mu=npow*re,reference_laminar_A=47,apparent_mu_from_A47_Pa_s=mue,Re_equivalent=ree,Pr_equivalent=pre,warning='Np*Re_equivalent=47 is a definition, not an independent validation; A47 requires matching impeller geometry.')
    if 'torque_resisting_Nm' in s:
        pc=2*math.pi*N*float(s['torque_resisting_Nm']);out.update(P_torque_W=pc,torque_vs_dissipation_relative_difference=(pc-P)/P)
    if 'P_heat_removed_stirrer_W' in s:
        q=float(s['P_heat_removed_stirrer_W']);ratio=q/J;pred=.005*ree**(1/3)*pre**(1/3)
        out.update(heat_transfer_fraction=ratio,Fig16_reference_fraction=pred,Fig16_relative_deviation=(ratio-pred)/pred)
    if all(x in s for x in ['cooled_area_m2','T_bulk_K','T_stirrer_K']):
        dt=float(s['T_bulk_K'])-float(s['T_stirrer_K']);a=float(s['cooled_area_m2'])
        if dt>0 and a>0 and 'P_heat_removed_stirrer_W' in s:
            h=float(s['P_heat_removed_stirrer_W'])/(a*dt);out.update(h_engineering_W_m2_K=h,Nu_engineering_using_d=h*d/k,Nu_definition_note='Engineering convention; paper does not specify enough detail to assert identical Nu normalization.')
    return out
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('json_file');a=ap.parse_args();print(json.dumps(metrics(json.load(open(a.json_file))),ensure_ascii=False,indent=2))
