#%%
from waverider_generator.generator import waverider as wr
import matplotlib.pyplot as plt
import numpy as np
import bezier.curve as bcurve

M_inf=8
beta=15
height=1.34
width=3
dp=[0.25,0.5,0.5,0.5]

waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=1000,n_shockwave=1000,n_planes=15,n_streamwise=11)


