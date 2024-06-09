#%%
import sys
sys.path.append('../')
from waverider_generator.generator import waverider as wr
from waverider_generator.cad_export import to_CAD

M_inf=5
beta=15
height=1.34
width=3
dp=[0.11,0.63,0,0.46]
n_planes=20
n_streamwise=30
delta_streamwise=0.1
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise,delta_streamwise=delta_streamwise)

# export cad
waverider_cad=to_CAD(waverider=waverider,sides='right',export=True,filename='waverider.step')
# visualise in jupyter notebook
waverider_cad