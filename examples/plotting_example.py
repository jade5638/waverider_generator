#%%
import sys
sys.path.append('../')
from waverider_generator.generator import waverider as wr
from waverider_generator.plotting_tools import Plot_Base_Plane,Plot_Leading_Edge
import matplotlib.pyplot as plt

M_inf=5
beta=15
height=1.34
width=3
dp=[0.11,0.63,0,0.46]
n_planes=40
n_streamwise=20
delta_streamwise=0.01
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise,delta_streamwise=delta_streamwise)
#%%
'''
PLOT BASE PLANE
'''
base_plane=Plot_Base_Plane(waverider=waverider)
#%%
'''
PLOT LEADING EDGE
'''
leading_edge=Plot_Leading_Edge(waverider=waverider)
plt.show()