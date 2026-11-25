"""Normalize a mapped nonnegative source on actual FLUID particles, not a uniform 2-D area."""
import numpy as np

def normalize_on_particles(shape,volume_m3,absorbed_power_W):
    q=np.asarray(shape,dtype=float);v=np.asarray(volume_m3,dtype=float)
    if q.shape!=v.shape or not np.isfinite(q).all() or not np.isfinite(v).all() or (v<=0).any() or (q<0).any() or absorbed_power_W<=0:raise ValueError('Invalid nonnegative source / positive volumes / power.')
    z=float(np.sum(q*v))
    if z<=0:raise ValueError('Zero integrated source.')
    return q*(absorbed_power_W/z)
