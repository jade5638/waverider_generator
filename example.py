#%%
from waverider_generator.generator import waverider as wr
import matplotlib.pyplot as plt
import numpy as np
import bezier.curve as bcurve

M_inf=5
beta=15
height=1.34
width=3
dp=[0.04,0.26,0.97,0.95]

waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=30,n_streamwise=11)

inters=waverider.local_intersections_us
lowers=np.column_stack([waverider.z_local_shockwave,waverider.y_bar_shockwave])


plt.figure()
for point1, point2 in zip(inters, lowers):
    x_values = [point1[0], point2[0]]
    z_values = [point1[1], point2[1]]
    plt.plot(x_values, z_values, 'bo-')
shockwave=np.vstack([np.array([0,0]),lowers,waverider.s_P4])
upper_surface=np.vstack([np.array([0,waverider.height]),inters,waverider.us_P3])
plt.plot(upper_surface[:,0],upper_surface[:,1],'r-')
plt.plot(shockwave[:, 0], shockwave[:, 1],'b-')
X2=waverider.us_P3
plt.plot(X2[0],X2[1],'bo')
plt.plot([0, 0],[0, waverider.height],'bo-')
plt.xlabel('Z-axis')
plt.ylabel('Y-axis')
plt.title('Osculating planes')
plt.gca().set_aspect('equal')

le=waverider.leading_edge
plt.figure()
plt.plot(le[:,2],-le[:,0],'b-')
plt.plot(-le[:,2],-le[:,0],'b-')
plt.gca().set_aspect('equal')

cone_centers=waverider.cone_centers
