from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from matplotlib.path import Path as MPath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
# Manually read contour vertices in rendered reference-image pixel coordinates.
# These are low-density visual traces, not original author numerical data.
traces={
'STIR_Fig8a':[
(.9,[(96,73),(76,76),(59,96),(49,121),(47,164),(50,206),(66,235),(95,251),(147,252),(187,241),(206,226),(211,207),(204,188),(191,174),(162,166),(168,149),(148,148),(127,152),(115,145),(112,127),(115,103),(113,89),(106,79)]),
(.9,[(596,74),(582,77),(573,90),(570,107),(571,128),(567,150),(548,153),(518,145),(494,136),(477,125),(473,133),(478,154),(474,175),(464,192),(464,214),(476,236),(500,246),(531,252),(566,253),(601,246),(621,227),(632,202),(636,170),(634,136),(628,108),(615,84)]),
(.8,[(193,275),(141,275),(95,275),(66,268),(49,252),(39,227),(35,184),(36,136),(42,94),(55,67),(72,50),(97,45),(132,50),(168,60),(201,71),(228,88),(241,109),(245,130),(241,146),(230,158),(238,175),(249,198),(261,223),(262,241),(249,256),(228,267)]),
(.8,[(494,275),(540,275),(586,275),(613,267),(634,246),(644,217),(646,175),(644,128),(636,86),(624,61),(603,47),(580,46),(550,54),(516,63),(479,74),(453,86),(440,103),(434,124),(435,154),(428,183),(420,205),(415,229),(418,248),(439,262),(465,271)]),
(.7,[(47,282),(29,253),(25,219),(24,175),(26,128),(30,87),(39,57),(52,41),(80,35),(122,36),(173,39),(219,46),(251,61),(263,83),(268,116),(267,155),(270,189),(280,220),(295,244)]),
(.7,[(391,243),(403,214),(404,187),(408,153),(412,120),(414,92),(423,69),(440,53),(466,44),(510,38),(552,36),(591,36),(620,43),(640,63),(651,94),(655,139),(656,183),(655,227),(645,258),(629,281)]),
(.6,[(293,31),(302,48),(307,75),(307,104),(301,127),(301,155),(306,184),(314,198),(335,205)]),
(.6,[(375,31),(366,49),(364,77),(363,105),(370,129)]),
(.6,[(249,288),(285,279),(318,274),(338,271),(363,272),(389,278),(420,288)])],
'STIR_Fig8c':[
(.9,[(251,154),(222,157),(194,174),(174,197),(163,224),(144,233),(119,228),(110,209),(115,184),(117,166),(109,167),(94,194),(81,224),(73,247)]),
(.9,[(295,155),(325,164),(351,179),(379,201),(403,228),(421,258),(435,281),(439,300),(434,325),(419,352),(402,367),(385,363),(376,368),(356,365),(352,351),(330,355),(309,361),(284,363),(258,358),(230,352),(204,340),(187,326),(177,307),(173,289)]),
(.8,[(252,197),(224,207),(204,218),(182,231),(157,243),(125,249)]),
(.8,[(300,201),(322,211),(344,223),(366,237),(387,249)]),
(.8,[(174,288),(191,307),(218,320),(246,327),(273,330),(302,322),(328,313),(350,307),(368,309)]),
(.7,[(320,232),(300,224),(279,217),(254,217),(233,222),(211,232),(188,241),(165,248)]),
(.7,[(181,286),(210,295),(234,307),(255,313),(280,314)]),
(.7,[(319,301),(340,294),(358,291),(369,293)]),
(.9,[(267,50),(350,69),(416,119),(462,192),(480,276),(461,360),(414,429),(340,470),(257,481),(175,463),(106,416),(63,344),(49,261),(66,178),(115,109),(187,64)])],
'STIR_Fig8d':[
(.4,[(243,245),(257,242),(271,247),(290,282),(277,295),(257,300),(239,292),(226,281)]),
(.45,[(215,250),(220,241),(235,225),(252,219),(273,222),(291,234),(305,251),(311,281),(299,305),(280,318),(256,323),(235,319),(214,306),(202,285)]),
(.5,[(173,249),(181,225),(200,204),(225,189),(256,182),(287,189),(316,207),(336,232),(344,254),(344,282),(333,310),(313,335),(282,351),(256,354),(225,350),(195,333),(177,310),(169,282)]),
(.55,[(119,250),(126,215),(148,180),(178,153),(216,135),(256,129),(295,137),(332,154),(363,186),(383,219),(396,254),(396,288),(386,326),(363,361),(332,389),(295,403),(257,407),(218,401),(182,384),(150,354),(129,320),(117,284)]),
(.55,[(255,57),(315,66),(369,91),(414,132),(447,184),(465,243),(469,301),(450,362),(413,412),(364,450),(312,469),(255,473),(193,462),(137,435),(91,393),(59,338),(46,277),(49,216),(70,157),(110,108),(163,73),(217,59)])],
'NC_Fig6':[
(.9,[(684,78),(707,83),(723,98),(728,123),(728,148),(713,159),(692,155),(678,146),(653,145),(643,128),(644,107),(631,99),(628,86),(647,79),(667,77)]),
(.5,[(666,47),(698,45),(723,53),(735,65),(751,80),(757,98),(771,116),(781,144),(786,170),(784,201),(772,224),(750,227)]),
(.5,[(659,224),(639,213),(623,201),(611,183),(591,172),(578,154),(575,132),(568,116),(567,97),(578,84),(592,73),(611,59),(632,52)]),
(.3,[(639,288),(627,280),(615,267),(606,257),(610,244),(597,232),(579,226),(562,215),(562,194),(548,180),(535,172),(522,172),(510,160),(500,146),(498,126),(502,104),(513,87),(531,76),(550,65),(560,64),(561,51),(579,44),(593,34),(621,28),(648,26),(689,26),(727,28),(742,42),(762,46),(780,61),(792,84),(796,114),(809,119),(808,151),(799,172),(798,209),(787,235),(779,254),(757,259),(744,276),(715,287)]),
(.1,[(660,381),(631,379),(607,372),(580,369),(561,355),(560,337),(548,335),(533,325),(519,325),(504,315),(486,298),(475,298),(457,281),(453,270),(430,255),(421,252),(416,239),(400,230),(386,229),(382,217),(381,195),(376,178),(374,160),(378,138),(378,123),(391,101),(406,89),(423,77),(440,70),(455,65),(456,59),(475,45),(491,33),(518,28),(546,24),(581,20),(583,5)]),
(.1,[(736,5),(736,20),(773,24),(790,32),(804,50),(811,71),(817,119),(838,123)]),
(.1,[(839,172),(818,172),(816,226),(806,271),(801,307),(787,336),(784,356),(777,367),(743,373)])]
}
traces['NC_Fig6']=[item for item in traces['NC_Fig6'] if item[0] != .5] # omit incompletely verified intermediate contour
bounds={'STIR_Fig8a':(15,27,667,306),'STIR_Fig8c':(22,20,510,509),'STIR_Fig8d':(14,16,500,509),'NC_Fig6':(47,0,840,585)}
polys={'STIR_Fig8a':[(320,27),(359,27),(359,93),(379,77),(400,78),(381,115),(361,132),(361,180),(400,216),(374,221),(342,196),(318,196)],'STIR_Fig8c':[(71,250),(399,250),(405,317),(370,318),(369,285),(111,285),(111,306),(77,306)],'STIR_Fig8d':[(76,253),(383,253),(388,316),(355,316),(351,282),(112,282),(112,305),(81,305)]}
summary=[]
for name,lines in traces.items():
 im=np.array(Image.open(ROOT/'figures'/f'{name}_reference.png'))
 l,t,r,b=bounds[name]; rows=[]
 for ci,(val,pts) in enumerate(lines):
  for order,(px,py) in enumerate(pts):
   rows.append(dict(curve_id=ci,point_order=order,pixel_x=px,pixel_y=py,x_unit=(px-l)/(r-l),z_or_y_unit=(b-py)/(b-t),value_normalized=val,T_K=val*1620 if name.startswith('STIR') else np.nan,source_type='manual_visual_contour_trace',coordinate_uncertainty_pixels=4,value_uncertainty_estimate=.025 if name.startswith('STIR') else .05))
 df=pd.DataFrame(rows)
 if name in polys:
  keep=~MPath(polys[name]).contains_points(df[['pixel_x','pixel_y']].values);df=df.loc[keep].copy()
 df.to_csv(ROOT/'data'/f'{name}_contour_points.csv',index=False)
 # unconnected observations for interpolation; no invented Dirichlet wall values
 xx,yy=np.meshgrid(np.linspace(0,1,161),np.linspace(0,1,121)); pts=df[['x_unit','z_or_y_unit']].values; v=df.value_normalized.values
 lin=griddata(pts,v,(xx,yy),method='linear'); extrap=np.isnan(lin); near=griddata(pts,v,(xx,yy),method='nearest'); field=np.where(extrap,near,lin)
 pixels=np.column_stack((l+xx.ravel()*(r-l),b-yy.ravel()*(b-t)))
 mask=np.ones(xx.size,dtype=bool)
 if name in polys:mask &= ~MPath(polys[name]).contains_points(pixels)
 if name in ('STIR_Fig8c','STIR_Fig8d'):mask &= ((xx.ravel()-.5)**2+(yy.ravel()-.5)**2<.25)
 dist=cKDTree(pts).query(np.column_stack((xx.ravel(),yy.ravel())))[0]
 if name=='NC_Fig6':
  # Low-source region below the 0.1 contour is not measured. An explicit decay prior is used ONLY for initialization.
  field=np.where(extrap,field*np.exp(-dist.reshape(xx.shape)/.12),field)
 out=pd.DataFrame(dict(x_unit=xx.ravel(),z_or_y_unit=yy.ravel(),value_normalized=field.ravel(),reference_visible_mask=mask.astype(int),nearest_extrapolation=extrap.ravel().astype(int),distance_to_trace_unit=dist))
 if name.startswith('STIR'):out['T_K']=out.value_normalized*1620
 else:out['Q_shape']=out.value_normalized # provisional only; full low-Q field is unobserved
 out.to_csv(ROOT/'data'/f'{name}_INITIAL_GUESS_grid.csv',index=False)
 fig,ax=plt.subplots(figsize=(8,5 if name.endswith('8a') else 7)); ax.imshow(im)
 for val, ps in lines:
  p=np.array(ps);ax.plot(p[:,0],p[:,1],'.-',lw=.65,ms=2,label=f'{val:g}')
 ax.set_title(name+' — manually traced contours (overlay)');fig.tight_layout();fig.savefig(ROOT/'figures'/f'{name}_trace_QA.png',dpi=160);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,4 if name.endswith('8a') else 6));f=field.copy(); f.ravel()[~mask]=np.nan
 z=ax.pcolormesh(xx,yy,f,vmin=0 if name=='NC_Fig6' else .3,vmax=1,cmap='inferno',shading='auto');fig.colorbar(z,ax=ax,label='Q/Qmax (seed)' if name=='NC_Fig6' else 'T/1620 K (seed)');ax.set(xlabel='normalized horizontal coordinate',ylabel='normalized vertical / second horizontal coordinate',title=name+' — reconstructed initial guess');fig.tight_layout();fig.savefig(ROOT/'figures'/f'{name}_INITIAL_GUESS.png',dpi=150);plt.close(fig)
 summary.append(dict(figure=name,points=len(df),curves=len(lines),grid_points=len(out),extrapolated_fraction=float(extrap.ravel()[mask].mean())))
pd.DataFrame(summary).to_csv(ROOT/'data'/'digitization_summary.csv',index=False)
re=np.geomspace(.1,10,31)
pd.DataFrame({'Re':re,'Np_reference':47/re,'source':'STIR Fig11 published correlation (not digitized markers)'}).to_csv(ROOT/'data'/'STIR_Fig11_analytic_reference.csv',index=False)
pd.DataFrame({'Re':re,'heat_fraction_over_Pr_one_third':.005*re**(1/3),'source':'STIR Fig16 published correlation (not digitized markers)'}).to_csv(ROOT/'data'/'STIR_Fig16_analytic_reference.csv',index=False)
print(pd.DataFrame(summary).to_string(index=False))
