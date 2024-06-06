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
dp=[0.17,0.98,0.87,0.64]

waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=20,n_streamwise=11)
le=waverider.leading_edge
#%%
ls_streams=waverider.lower_surface_streams
i=11
# x=np.linspace(0,width,20)

print(ls_streams[i])
print("\n")
print(le[i,:])
# print(x[i])
#%%

inters=waverider.local_intersections_us
lowers=np.column_stack([waverider.z_local_shockwave,waverider.y_local_shockwave])

plt.figure()
for point1, point2 in zip(inters, lowers):
    x_values = [point1[0], point2[0]]
    z_values = [point1[1], point2[1]]
    plt.plot(x_values, z_values, 'bo-')
shockwave=np.vstack([np.array([0,0]),lowers,waverider.s_P4])
upper_surface=np.vstack([np.array([0,waverider.height]),inters,waverider.us_P3])
ls_streams=waverider.lower_surface_streams
lower_surface = np.vstack([stream[-1,:] for stream in ls_streams])
z_ls=lower_surface[:,2]
y_ls = lower_surface[:,1]+height
plt.plot(z_ls,y_ls,'--ok')
plt.plot(upper_surface[:,0],upper_surface[:,1],'r-')
plt.plot(shockwave[:, 0], shockwave[:, 1],'b-')
X2=waverider.us_P3
plt.plot(X2[0],X2[1],'bo')
plt.plot([0, 0],[0, waverider.height],'bo-')
plt.xlabel('Z-axis')
plt.ylabel('Y-axis')
plt.title('Osculating planes')
plt.gca().set_aspect('equal')
#%%
le=waverider.leading_edge
plt.figure()
plt.plot(le[:,2],-le[:,0],'b-')
plt.plot(-le[:,2],-le[:,0],'b-')
plt.gca().set_aspect('equal')

x=waverider.upper_surface_x
y=waverider.upper_surface_y
z=waverider.upper_surface_z
# x_tip_row = np.full((1, 11), waverider.length)
# y_tip_row=np.full((1,11),height*dp[1]-height)
# z_tip_row=np.full((1,11),width)
# x = np.vstack([x, x_tip_row])
# y = np.vstack([y, y_tip_row])
# z = np.vstack([z, z_tip_row])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# for i in range(21):
#     ax.plot(upper_surface_x[i,:], upper_surface_y[i,:], upper_surface_z[i,:])
ax.plot_surface(x,y,z,color='blue')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.gca().set_aspect('equal')


# streams=waverider.streams
# x_ls= np.array([line[:, 0] for line in streams])
# y_ls = np.array([line[:, 1] for line in streams])
# z_ls = np.array([line[:, 2] for line in streams])
# ax.plot_surface(x_ls, y_ls, z_ls,color='red')
plt.show()

#%%
from waverider_generator.flowfield import shock_angle,cone_angle
import numpy as np
M=5
gamma=1.4
# beta=shock_angle(M,theta,gamma)
beta=20
theta_comp=cone_angle(M,beta,gamma)
# %%
le=waverider.leading_edge
i=15
x=np.linspace(0,width,20)

print(ls_streams[8])
print("\n")
print(le[i,:])
print(x[i])
