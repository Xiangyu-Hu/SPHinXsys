"""Structural checks only; source-overlay visual review remains necessary."""
from pathlib import Path
import csv,json,hashlib,math
root=Path(__file__).resolve().parents[1]/'reference'
meta=json.loads((root/'digitization_provenance.json').read_text())
for key,sha in meta['image_sha256'].items():
    assert hashlib.sha256((root/(key+'_native.png')).read_bytes()).hexdigest()==sha
assert meta['absolute_Tmax_K'] is None
assert not meta['fig8']['absolute_F1_ready']
rows=list(csv.DictReader((root/'fig8_labelled_contour_points.csv').open()))
assert len(rows)==4928
assert set(float(r['T_over_Tmax']) for r in rows)=={.5,.6,.7,.8,.9}
for r in rows:
    assert 0<=float(r['r_over_R'])<=1 and 0<=float(r['z_over_H'])<=1
v=list(csv.DictReader((root/'fig7_direction_samples.csv').open()))
assert len(v)==20 and all(r['speed_m_s']=='' for r in v)
for r in v:
    assert abs(math.hypot(float(r['display_direction_dr']),float(r['display_direction_dz']))-1)<1e-12
g=list(csv.DictReader((root/'fig8_partial_linear_reconstruction.csv').open()))
assert len(g)==10201 and all(r['usable_for_absolute_F1']=='0' for r in g)
assert any(r['status']=='OUTSIDE_OBSERVATION_HULL' for r in g)
assert any(r['status']=='HOT_CORE_UNRESOLVED_ABOVE_0P9' for r in g)
print('PASS: source hashes, counts, coordinates, levels, explicit missing speed/Kelvin, reconstruction masks')
