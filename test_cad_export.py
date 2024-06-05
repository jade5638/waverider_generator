#%%
from waverider_generator.generator import waverider as wr
import matplotlib.pyplot as plt
import numpy as np
import bezier.curve as bcurve
from mpl_toolkits.mplot3d import Axes3D
import cadquery as cq
from cadquery import exporters 

M_inf=5
beta=15
height=1.34
width=3
dp=[0.33,0,0.06,0.37]
n_planes=20
n_streamwise=11

waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise)

streams=waverider.upper_surface_streams
le=waverider.leading_edge
surface_points = []
for i in range(len(streams)):
    for j in range(1, streams[i].shape[0]-1):
        surface_points.append(tuple(streams[i][j]))

#%%
print(streams[0])