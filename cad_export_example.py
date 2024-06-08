#%%
from waverider_generator.generator import waverider as wr
from waverider_generator.cad_export import to_CAD
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import cadquery as cq
from cadquery import exporters 

M_inf=5
beta=15
height=1.34
width=3
dp=[0.04,0.98,0.87,0.64]
n_planes=50
n_streamwise=10
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise)
#%%
waverider_cad=to_CAD(waverider=waverider)
waverider_cad