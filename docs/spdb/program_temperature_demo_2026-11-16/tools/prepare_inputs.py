"""Read supplied ASCII STL in mm; create two aligned, rotated closed bodies in SI."""
from pathlib import Path
import json, math, hashlib
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
def load(path):
    vs=[tuple(float(x) for x in l.split()[1:]) for l in path.read_text().splitlines() if l.strip().startswith('vertex')]
    assert len(vs)%3==0
    return [vs[i:i+3] for i in range(0,len(vs),3)]
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def info(ts):
    vs=[v for t in ts for v in t]; e=Counter(); vol=0
    for a,b,c in ts:
        vol+=sum(x*y for x,y in zip(a,cross(b,c)))/6
        for v,w in [(a,b),(b,c),(c,a)]:e[tuple(sorted((v,w)))]+=1
    return {'bounds':[[min(v[k] for v in vs) for k in range(3)],[max(v[k] for v in vs) for k in range(3)]], 'signed_volume':vol,'triangles':len(ts),'all_edges_two_faces':all(n==2 for n in e.values())}
def write(path,ts):
    with path.open('w') as f:
        f.write('solid specimen\n')
        for a,b,c in ts:
            n=cross(sub(b,a),sub(c,a)); m=math.sqrt(sum(x*x for x in n)); n=tuple(x/m for x in n) if m else (0,0,0)
            f.write(' facet normal '+' '.join(f'{x:.12g}' for x in n)+'\n  outer loop\n')
            for v in (a,b,c):f.write('   vertex '+' '.join(f'{x:.12g}' for x in v)+'\n')
            f.write('  endloop\n endfacet\n')
        f.write('endsolid specimen\n')
src=ROOT/'geometry/specimen_original_mm.stl';ts=load(src)
# Original slot opens at z=20 mm. Translate lower to [-20,0] mm.
lower=[[(x*.001,y*.001,(z-20)*.001) for x,y,z in t] for t in ts]
# Rotate 180 degrees about x (det=+1), then translate original z=20 to 0.
upper=[[(x*.001,-y*.001,(20-z)*.001) for x,y,z in t] for t in ts]
write(ROOT/'geometry/lower_m.stl',lower);write(ROOT/'geometry/upper_m.stl',upper)
report={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'original_units':'mm_assumed_by_user','generated_units':'m','original':info(ts),'lower':info(lower),'upper':info(upper),'nominal_slot_each_mm':[10,10,10],'nominal_combined_void_mm':[10,10,20],'note':'Two separate bodies; interface at z=0; no welding. Exact source chamfers preserved.'}
assert report['lower']['all_edges_two_faces'] and report['upper']['all_edges_two_faces']
assert report['lower']['signed_volume']>0 and report['upper']['signed_volume']>0
(ROOT/'geometry/geometry_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
