#%%
from waverider_generator.generator import waverider as wr
from waverider_generator.plotting_tools import Plot_Base_Plane,Plot_Leading_Edge,Plot_CAD
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import cadquery as cq
from cadquery import exporters 

M_inf=5
beta=15
height=1.34
width=3
dp=[0.11,0.63,0,0.46]
n_planes=20
n_streamwise=10
delta_streamwise=0.1
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise,delta_streamise=delta_streamwise)
le=waverider.leading_edge
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
#%%
'''
PLOT CAD'''
CAD_plot=Plot_CAD(waverider=waverider)
CAD_plot.plot()

plt.plot()