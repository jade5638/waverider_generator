#%%
from waverider_generator.generator import waverider as wr
from waverider_generator.plotting_tools import Plot_Base_Plane,Plot_Leading_Edge
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
waverider=wr(M_inf=M_inf,beta=beta,height=height,width=width,dp=dp,n_upper_surface=10000,n_shockwave=10000,n_planes=n_planes,n_streamwise=n_streamwise)
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
plt.show()
#%%
# x=waverider.upper_surface_x
# y=waverider.upper_surface_y
# z=waverider.upper_surface_z
# # x_tip_row = np.full((1, 11), waverider.length)
# # y_tip_row=np.full((1,11),height*dp[1]-height)
# # z_tip_row=np.full((1,11),width)
# # x = np.vstack([x, x_tip_row])
# # y = np.vstack([y, y_tip_row])
# # z = np.vstack([z, z_tip_row])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # for i in range(21):
# #     ax.plot(upper_surface_x[i,:], upper_surface_y[i,:], upper_surface_z[i,:])
# ax.plot_surface(x,y,z,color='blue')
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# plt.gca().set_aspect('equal')


# # streams=waverider.streams
# # x_ls= np.array([line[:, 0] for line in streams])
# # y_ls = np.array([line[:, 1] for line in streams])
# # z_ls = np.array([line[:, 2] for line in streams])
# # ax.plot_surface(x_ls, y_ls, z_ls,color='red')
# plt.show()

# #%%
# from waverider_generator.flowfield import shock_angle,cone_angle
# import numpy as np
# M=5
# gamma=1.4
# # beta=shock_angle(M,theta,gamma)
# beta=20
# theta_comp=cone_angle(M,beta,gamma)
# # %%
# le=waverider.leading_edge
# i=15
# x=np.linspace(0,width,20)

# print(ls_streams[8])
# print("\n")
# print(le[i,:])
# print(x[i])
