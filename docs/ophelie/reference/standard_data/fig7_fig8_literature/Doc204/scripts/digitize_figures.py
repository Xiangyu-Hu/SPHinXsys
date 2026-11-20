"""Reproduce the Fig.8 line extraction and Fig.7 manual direction annotations.

Input images are native 1-bit figure rasters extracted from the user-provided PDF.
No velocity magnitude or absolute temperature is recovered from these images.
Run from the deliverable root: python3 scripts/digitize_figures.py
Requires numpy, scipy, matplotlib, Pillow.
"""
from pathlib import Path
import csv, json, hashlib
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'reference'; QA=ROOT/'qa'
QA.mkdir(exist_ok=True)
CAL={
 'fig7':dict(x0=145,x1=2309,y_top=12,y_bottom=1606,width=2318,height=1882),
 'fig8':dict(x0=145,x1=2326,y_top=12,y_bottom=1618,width=2335,height=1907),
}
def writecsv(name,fields,rows):
    with (REF/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def coord(x,y,c):
    return ((x-c['x0'])/(c['x1']-c['x0']),
            (c['y_bottom']-y)/(c['y_bottom']-c['y_top']))
def bg(ax,key):
    c=CAL[key];im=Image.open(REF/(key+'_native.png')).convert('L')
    ax.imshow(im,cmap='gray',vmin=0,vmax=255,
      extent=[-c['x0']/(c['x1']-c['x0']),
              (c['width']-c['x0'])/(c['x1']-c['x0']),
              (c['y_bottom']-c['height'])/(c['y_bottom']-c['y_top']),
              c['y_bottom']/(c['y_bottom']-c['y_top'])],aspect='auto')
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='r/R',ylabel='z/H')

# Line components are manually assigned ONLY where the printed level is readable.
# Components touching the top frame are split by the explicit crop. Unlabelled
# outer curves are deliberately not assigned inferred temperatures.
a=np.array(Image.open(REF/'fig8_native.png').convert('L'))<128
a[:25]=False;a[1600:]=False;a[:,:160]=False;a[:,2310:]=False
lab,n=ndi.label(a,np.ones((3,3)))
component_levels={11:.9,5:.8,6:.8,4:.7,7:.7,26:.6,8:.6,30:.5,12:.5}
expected_sizes={11:9942,5:7732,6:4190,4:7904,7:5397,26:3196,8:7066,30:1509,12:7920}
rows=[];c=CAL['fig8']
for component,level in component_levels.items():
    yy,xx=np.where(lab==component)
    assert len(xx)==expected_sizes[component], (component,len(xx))
    bins={}
    for x,y in zip(xx,yy):bins.setdefault((int(x)//5,int(y)//5),[]).append((x,y))
    for k,pts in sorted(bins.items()):
        x,y=np.mean(pts,axis=0);r,z=coord(x,y,c)
        rows.append(dict(component_id=component,r_over_R=r,z_over_H=z,
          T_over_Tmax=level,x_px=x,y_px=y,
          estimated_coordinate_error_r=0.005,estimated_coordinate_error_z=0.007,
          provenance='DIGITIZED_LABELLED_CONTOUR'))
writecsv('fig8_labelled_contour_points.csv',list(rows[0]),rows)

# Intersections use observed pixels, not the reconstructed field. Nearby columns
# may contain two branches of the same isotherm; preserve both.
crossings=[]
for cutaxis,cut in [('r_over_R',.5),('z_over_H',.5)]:
    for lev in [.5,.6,.7,.8,.9]:
        pts=np.array([[p['r_over_R'],p['z_over_H']] for p in rows if p['T_over_Tmax']==lev])
        axis=0 if cutaxis=='r_over_R' else 1
        sel=pts[np.abs(pts[:,axis]-cut)<.0035]
        if not len(sel):continue
        vals=np.sort(sel[:,1-axis]);groups=np.split(vals,np.where(np.diff(vals)>.025)[0]+1)
        for g in groups:
            v=float(np.median(g))
            crossings.append(dict(cut_axis=cutaxis,cut_value=cut,
              r_over_R=cut if axis==0 else v,z_over_H=v if axis==0 else cut,
              T_over_Tmax=lev,coordinate_uncertainty=.01,
              source='OBSERVED_CONTOUR_INTERSECTION'))
writecsv('fig8_midline_contour_intersections.csv',list(crossings[0]),crossings)

# Fig.7 manual readings: positions are representative locations on visible arrows;
# orientation is rounded to a 15-degree step, not a measured velocity value.
# Angle is measured in the DISPLAYED normalized (r/R,z/H) plane. This is not
# physical atan2(uz,ur); physical directions require R/H correction below.
samples=[
 (.45,.94,180,25),(.55,.90,165,25),(.27,.88,225,25),
 (.16,.75,250,20),(.075,.56,270,20),(.14,.39,285,25),
 (.22,.28,315,25),(.32,.39,0,20),(.43,.51,15,20),
 (.51,.65,15,20),(.54,.74,15,25),(.67,.73,45,25),
 (.77,.74,75,20),(.82,.57,90,20),(.82,.88,90,20),
 (.91,.95,345,30),(.96,.85,300,30),(.96,.73,270,25),
 (.90,.65,225,30),(.87,.80,105,30),
]
vectors=[]
for i,(r,z,ang,unc) in enumerate(samples,1):
    theta=np.deg2rad(ang);dr,dz=np.cos(theta),np.sin(theta)
    physical=np.array([.25*dr,.185*dz]);physical/=np.linalg.norm(physical)
    vectors.append(dict(sample_id=i,r_over_R=r,z_over_H=z,
      display_direction_dr=dr,display_direction_dz=dz,
      direction_ur_assuming_R025_H0185=float(physical[0]),
      direction_uz_assuming_R025_H0185=float(physical[1]),
      display_angle_deg=ang,estimated_angle_uncertainty_deg=unc,
      estimated_position_uncertainty=.025,
      speed_m_s='',source='MANUAL_APPROXIMATE_ARROW_ORIENTATION'))
writecsv('fig7_direction_samples.csv',list(vectors[0]),vectors)
features=[
 dict(feature='main_cell_center_search_region',r_min=.38,r_max=.58,z_min=.72,z_max=.88,
      type='APPROXIMATE_VISUAL_ROI_NOT_EXACT_STAGNATION_POINT'),
 dict(feature='peripheral_cell_center_search_region',r_min=.85,r_max=.97,z_min=.72,z_max=.93,
      type='APPROXIMATE_VISUAL_ROI_NOT_EXACT_STAGNATION_POINT'),
 dict(feature='upper_upwelling_band',r_min=.76,r_max=.86,z_min=.60,z_max=.95,
      type='APPROXIMATE_VISUAL_REGION'),
 dict(feature='axis_side_downflow',r_min=.00,r_max=.16,z_min=.30,z_max=.85,
      type='APPROXIMATE_VISUAL_REGION'),
 dict(feature='peripheral_downflow',r_min=.94,r_max=1.0,z_min=.60,z_max=.90,
      type='APPROXIMATE_VISUAL_REGION'),
]
writecsv('fig7_topology_regions.csv',list(features[0]),features)

# A PARTIAL reconstruction: only interpolates in the convex hull of observed
# contour coordinates. No fake boundary values, Tmax anchor, or kelvin conversion.
pts=np.array([[p['r_over_R'],p['z_over_H']] for p in rows]);vals=np.array([p['T_over_Tmax'] for p in rows])
interp=LinearNDInterpolator(pts,vals,fill_value=np.nan)
rgrid=np.linspace(0,1,101);zgrid=np.linspace(0,1,101)
rr,zz=np.meshgrid(rgrid,zgrid);tt=np.asarray(interp(rr,zz))
grid=[]
for r,z,t in zip(rr.ravel(),zz.ravel(),tt.ravel()):
    # Regions inside the 0.9 isotherm are unconstrained ABOVE 0.9. Do not report
    # the piecewise linear plateau as a uniquely recovered hot-core temperature.
    hot=bool(np.isfinite(t) and t>=.899999)
    grid.append(dict(r_over_R=r,z_over_H=z,T_over_Tmax='' if not np.isfinite(t) else float(t),
      status='OUTSIDE_OBSERVATION_HULL' if not np.isfinite(t) else
             ('HOT_CORE_UNRESOLVED_ABOVE_0P9' if hot else 'INTERPOLATED_NOT_OBSERVED'),
      usable_for_absolute_F1=0))
writecsv('fig8_partial_linear_reconstruction.csv',list(grid[0]),grid)

fig,ax=plt.subplots(figsize=(11,8));bg(ax,'fig8')
for lev in [.5,.6,.7,.8,.9]:
    q=np.array([[p['r_over_R'],p['z_over_H']] for p in rows if p['T_over_Tmax']==lev])
    ax.scatter(q[:,0],q[:,1],s=2,label=f'T/Tmax={lev:.1f}')
ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=9);ax.set_title('Fig.8: extracted labelled contour pixels (no fitted curves)')
fig.tight_layout();fig.savefig(QA/'fig8_contour_overlay.png',dpi=170,bbox_inches='tight');plt.close(fig)

fig,ax=plt.subplots(figsize=(11,8));bg(ax,'fig7')
ax.quiver([s[0] for s in samples],[s[1] for s in samples],
          [np.cos(np.deg2rad(s[2])) for s in samples],
          [np.sin(np.deg2rad(s[2])) for s in samples],
          color='#de3428',angles='xy',scale_units='xy',scale=18,width=.003)
for i,s in enumerate(samples,1):ax.text(s[0]+.006,s[1]+.008,str(i),fontsize=7,color='#c00000')
from matplotlib.patches import Rectangle
for f in features[:2]:
    ax.add_patch(Rectangle((f['r_min'],f['z_min']),f['r_max']-f['r_min'],f['z_max']-f['z_min'],
                          fill=False,color='#1369d4',linewidth=1.4))
ax.set_title('Fig.7: approximate arrow directions; boxes are center SEARCH regions')
fig.tight_layout();fig.savefig(QA/'fig7_direction_overlay.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(10,7));pm=ax.pcolormesh(rr,zz,np.ma.masked_invalid(tt),shading='auto',vmin=.5,vmax=.9)
fig.colorbar(pm,ax=ax,label='T/Tmax (linear interpolation of observed contours only)')
ax.set(xlabel='r/R',ylabel='z/H',title='PARTIAL reconstruction: blank exterior; hot core above 0.9 unresolved')
fig.tight_layout();fig.savefig(QA/'fig8_partial_reconstruction.png',dpi=150);plt.close(fig)

meta=dict(source_title='Numerical modelling of natural convection in molten glass heated by induction',
 source_pdf_sha256='6dddf36ebe7078dd9d5dcf90f852e99032fe1cdb1aa07df30606105bf519c92b',
 doi='10.1016/j.cep.2007.01.013',pdf_page_1based=5,printed_page=453,
 normalization='Tstar=T/Tmax; NOT (T-Tcold)/(Thot-Tcold)',
 absolute_Tmax_K=None,case_side_h_W_m2_K=None,case_match='PARTIAL_UNRESOLVED',
 R_m_assumed=.25,H_m_assumed=.185,geometry_provenance='R user-confirmed; H thesis supplement, not read from image',
 calibration=CAL,coordinate_origin='r=0 symmetry axis; z=0 bottom; full plot height assumed H',
 image_sha256={k:hashlib.sha256((REF/(k+'_native.png')).read_bytes()).hexdigest() for k in CAL},
 fig8=dict(status='DIGITIZED_LABELLED_CONTOURS',point_count=len(rows),
   levels=[.5,.6,.7,.8,.9],estimated_reading_error_r=.005,estimated_reading_error_z=.007,
   uncertainty_note='Conservative manual coordinate allowance, not statistical confidence interval. Level text read exactly; contour spacing does not bound interpolation error.',
   reconstruction_status='PARTIAL_INTERPOLATION_NOT_AUTHOR_DATA',absolute_F1_ready=False),
 fig7=dict(status='DIGITIZED_APPROXIMATE_DIRECTIONS_AND_REGIONS',direction_sample_count=len(vectors),
   velocity_magnitude_available=False,streamfunction_available=False,
   center_ROIs_are_not_exact_vortex_centers=True,
   orientation_note='Manual visual estimates rounded to 15 degrees in normalized plot coordinates; no arrow-length speed inference.'),
 limitations=['Do not generate exact velocity/streamfunction from Fig.7.',
 'Do not feed partial reconstruction directly to a Kelvin-only input.',
 'Do not use a reconstructed T field both to force the simulation and independently validate its temperature.',
 'No temperature or velocity accuracy claim outside observed image constraints.'])
(REF/'digitization_provenance.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'contour_points':len(rows),'directions':len(vectors),'crossings':len(crossings),
 'grid_rows':len(grid),'finite_grid':int(np.isfinite(tt).sum())},indent=2))
