#!/usr/bin/env python3
"""Bind user-verified time/iteration values; does not run SPHinXsys."""
import argparse,json,math,hashlib
from pathlib import Path

def bind(c):
 def n(k,default=None):
  v=c.get(k,default)
  if v is None:v=default
  if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
   raise ValueError(f'{k}: finite real value must be supplied')
  return float(v)
 def grid(a,b,count):
  if b<a:raise ValueError('reversed interval')
  if b==a:return [a]
  return [a+(b-a)*j/(count-1) for j in range(count)]
 def unique(xs):
  ys=[]
  for v in sorted(xs):
   if not ys or not math.isclose(v,ys[-1],rel_tol=1e-12,abs_tol=1e-12):ys.append(v)
  return ys
 def period():
  v=n('P_ref')
  if v<=0:raise ValueError('P_ref must be positive')
  return v
 coord=c.get('coordinate','time')
 policy=c['probe_policy']; fp=c.get('field_policy'); field_on=c.get('field_enabled',False)
 mode=c.get('mode')
 if mode=='benchmark' and field_on:raise ValueError('benchmark must disable VTP')
 aux=[];thermal=[]
 if coord=='iteration':
  raw=n('iterations')
  if raw<1 or not raw.is_integer():raise ValueError('iterations must be positive integer')
  total=int(raw)
  if policy not in ['iterations100','uniform100']:raise ValueError('iteration axis accepts iterations100/uniform100 only')
  probes=sorted({round(x) for x in grid(0,total,min(100,total+1))})
  if fp=='iterations8':fields=[0,1,.01*total,.05*total,.2*total,.5*total,.8*total,total]
  elif fp=='transient5':fields=[f*total for f in [.2,.4,.6,.8,1]]
  else:
   if field_on:raise ValueError('unsupported iteration field policy')
   fields=[]
  fields=sorted({min(total,max(0,round(x))) for x in fields}) if field_on else []
  sparse=[0,round(total/2),total]
  start=0;end=total
 elif coord=='time':
  total=n('T')
  if total<=0:raise ValueError('T must be positive')
  start=n('S',0)
  end=n('E',total)
  if policy=='cycles100' or (policy=='industrial51' and mode!='benchmark'):
   end=start+4*period()
   if c.get('E') is not None and not math.isclose(n('E'),end,rel_tol=1e-10):
    raise ValueError('cycles100/industrial showcase E must equal S+4P_ref')
   probes=grid(start,end,100)
  elif policy in ['uniform100','event100']:probes=grid(start,end,100)
  elif policy in ['window51','industrial51']:probes=grid(start,end,51)
  else:raise ValueError('unknown time probe_policy')
  if not 0<=start<end<=total:raise ValueError('require 0 <= S < E <= T')
  if policy=='event100':
   a=n('event_start');b=n('event_end')
   if not start<=a<b<=end:raise ValueError('event interval outside probe window')
   probes=unique(probes+grid(a,b,21))
  if c.get('thermal_start') is not None or c.get('thermal_end') is not None:
   h0=n('thermal_start');h1=n('thermal_end')
   if not 0<=h0<h1<=total:raise ValueError('invalid thermal interval')
   thermal=grid(h0,h1,100)
  if c.get('include_startup_health',False) and start>0:aux=grid(0,start,51)
  if not field_on:fields=[]
  elif fp=='transient5':fields=[f*total for f in [.2,.4,.6,.8,1]]
  elif fp in ['cycle8','industrial8']:
   p=period();fs=n('field_start',start+3*p)
   fields=[fs+k*p/8 for k in range(8)]
   if fp=='industrial8' and thermal:fields.append(thermal[-1])
  elif fp=='evolution8':
   fs=n('field_start',start);fe=n('field_end',end)
   if fe<=fs:raise ValueError('field_end must exceed field_start')
   fields=[fs+f*(fe-fs) for f in [0,.1,.25,.4,.55,.7,.85,1]]
  elif fp=='event8':
   a=n('event_start');b=n('event_end')
   if not 0<=a<b<=total:raise ValueError('invalid field event interval')
   fields=[0,a/2,a,a+(b-a)/3,a+2*(b-a)/3,b,b+(total-b)/2,total]
  else:raise ValueError('unknown time field_policy')
  fields=unique(fields);sparse=[start,(start+end)/2,end]
  if any(t<0 or t>total for t in fields):raise ValueError('field times outside [0,T]')
 else:raise ValueError('coordinate must be time or iteration')
 if len(fields)>10:raise ValueError('more than 10 field groups')
 canonical=json.dumps(c,sort_keys=True,separators=(',',':'),ensure_ascii=False)
 return dict(case_id=c['case_id'],mode=mode,record_kind='bound_configuration_not_measured',
  coordinate=coord,coordinate_unit='iteration' if coord=='iteration' else 's (caller must declare pseudo-time separately)',
  input_sha256=hashlib.sha256(canonical.encode()).hexdigest(),
  probe_policy=policy,probe_targets=probes,probe_target_count=len(probes),
  startup_health_targets=aux,thermal_probe_targets=thermal,
  sparse_expensive_diagnostic_targets=sparse,field_policy=fp,field_targets=fields,field_group_count=len(fields),
  telemetry_interval_wall_s=1,performance_flush_wall_s=None,
  sampling_note='Target lists, not computed states. Record actual time; no duplicated state for crossed targets. Do not alter solver dt.',
  field_note='One group may contain multiple body VTP files. Verify field mappings and storage estimate in preflight.')

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('input_json',type=Path);p.add_argument('output_json',type=Path)
 a=p.parse_args()
 try:r=bind(json.loads(a.input_json.read_text(encoding='utf-8')))
 except (ValueError,KeyError,TypeError) as e:p.error(str(e))
 # Exclusive create avoids replacing a schedule that may already back real runs.
 a.output_json.parent.mkdir(parents=True,exist_ok=True)
 with a.output_json.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
 print(json.dumps(dict(output=str(a.output_json),probe_target_count=r['probe_target_count'],field_group_count=r['field_group_count']),ensure_ascii=False))
if __name__=='__main__':main()
