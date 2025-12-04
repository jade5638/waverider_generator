#%%

import sys
# Assumes that cwd is the repo folder
sys.path.append('waverider_generator')

import matplotlib.pyplot as plt 
from waverider_generator import *


try :
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 11})
except :
    pass
    
if __name__ == "__main__":


    geo_params = GeometricParameters(X1 = 0.2, 
                                    X2 = 0.6, 
                                    X3 = 0.2, 
                                    X4 = 0.1, 
                                    w = 4.2, 
                                    h = 1.876)

    flow_params = FlowParameters(M_design = 5.561, 
                                 betaDeg = 15, 
                                 gamma = 1.4)

    waverider = Waverider(geo_params=geo_params, flow_params=flow_params)

    streams, curves = waverider.Build(n_planes=50, n_streamline_pts=20, dx_LS_streamline=0.05)

    solid = waverider.to_CAD(streams=streams, 
                             sides = 'both', 
                             export_filename='waverider_example.step')

    sc = curves["SC"]
    usc = curves["USC"]
    lsc = curves["LSC"]
    le = curves['LE']

    fig, axs = plt.subplots(1, 3)

    # plot base plane
    axs[0].plot(sc.z, sc.y, 'b-o')
    axs[0].plot(usc.z, usc.y, 'r-o')
    axs[0].plot(lsc.z, lsc.y, 'k-o')

    for i in range(len(sc)):
        axs[0].plot([sc[i].z, usc[i].z], [sc[i].y, usc[i].y], 'k--', alpha=0.7)
                
    axs[0].set_xlabel('z [m]')
    axs[0].set_ylabel('y [m]')

    axs[0].set_title('Waverider Base Plane')

    axs[0].set_aspect('equal')

    leg = axs[0].legend(['Shockwave Curve', 'Upper Surface Curve', 'Lower Surface Curve', 'Osculating Planes'], loc = 'upper right')
    leg.set_draggable(True)

    # plot waverider from the top
    axs[1].plot(le.z, le.x, 'b-')
    axs[1].plot([0, 0], [le.x[0], usc.x[0]], 'k-')
    axs[1].plot(usc.z, usc.x, 'r-')

    axs[1].invert_yaxis()

    axs[1].set_xlabel('z [m]')
    axs[1].set_ylabel('x [m]')

    axs[1].set_title('Waverider Top View')

    axs[1].set_aspect('equal')

    leg = axs[1].legend(['Leading Edge', 'Symmetry Plane', 'Base Plane'], loc = 'upper right')
    leg.set_draggable(True)

    # plot upper and lower surface streams in 3D
    pos = axs[2].get_position()
    axs[2].remove()

    ax3d = fig.add_axes(pos, projection='3d')
    for us_stream, ls_stream in zip(streams.US_Streams, streams.LS_Streams):
        ax3d.plot(us_stream.z, us_stream.x, us_stream.y, 'r-')
        ax3d.plot(ls_stream.z, ls_stream.x, ls_stream.y, 'k-')

    leg = ax3d.legend(['Upper Surface Streams', 'Lower Surface Streams'], loc = 'upper right')
    leg.set_draggable(True)
        
    ax3d.set_xlabel('z [m]')
    ax3d.set_ylabel('x [m]')
    ax3d.set_zlabel('y [m]')

    ax3d.invert_xaxis()

    ax3d.set_title('Waverider Streams 3D View')

    ax3d.set_aspect('equal')

    plt.show()
