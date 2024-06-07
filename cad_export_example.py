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
dp=[0.21,0.31,0.21,0.33]
n_planes=50
n_streamwise=20
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise)

to_CAD(waverider=waverider)