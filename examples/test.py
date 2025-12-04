#%%

import sys
# Assumes that cwd is the repo folder
sys.path.append('waverider_generator')

import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec 
from waverider_generator import *
    
if __name__ == "__main__":

    # Define geometric parameters
    geo_params = GeometricParameters(X1 = 0.2, 
                                    X2 = 0.6, 
                                    X3 = 0.2, 
                                    X4 = 0.1, 
                                    w = 4.2, 
                                    h = 1.876)

    # Define flow parameters
    flow_params = FlowParameters(M_design = 5.561, 
                                 betaDeg = 15, 
                                 gamma = 1.4)

    # Instantiate the waverider
    waverider = Waverider(geo_params=geo_params, flow_params=flow_params)

    # Build the waverider
    streams, curves = waverider.Build(n_planes=50, n_streamline_pts=20, dx_LS_streamline=0.05)

    # Use streams to export to CAD
    solid = waverider.to_CAD(streams=streams, 
                             sides = 'both', 
                             export_filename='waverider_example.step')

    # Extract various curves
    sc = curves["SC"]
    usc = curves["USC"]
    lsc = curves["LSC"]
    le = curves['LE']

    fig = plt.figure()
    gs = GridSpec(2, 2)

    ax1 = fig.add_subplot(gs[0,0])

    # plot base plane
    ax1.plot(sc.z, sc.y, 'b-o')
    ax1.plot(usc.z, usc.y, 'r-o')
    ax1.plot(lsc.z, lsc.y, 'k-o')

    # osculating planes
    for i in range(len(sc)):
        ax1.plot([sc[i].z, usc[i].z], [sc[i].y, usc[i].y], 'k--', alpha=0.7)
                
    ax1.set_xlabel('z [m]')
    ax1.set_ylabel('y [m]')

    ax1.set_title('Waverider Base Plane')

    ax1.set_aspect('equal')
    
    leg = ax1.legend(['Shockwave Curve', 'Upper Surface Curve', 'Lower Surface Curve', 'Osculating Planes'], loc = 'upper right')
    leg.set_draggable(True)

    # plot waverider from the top
    ax2 = fig.add_subplot(gs[:,1])
    
    ax2.plot(le.z, le.x, 'b-')
    ax2.plot([0, 0], [le.x[0], usc.x[0]], 'k-')
    ax2.plot(usc.z, usc.x, 'r-')

    ax2.invert_yaxis()

    ax2.set_xlabel('z [m]')
    ax2.set_ylabel('x [m]')

    ax2.set_title('Waverider Top View')

    ax2.set_aspect('equal')

    leg = ax2.legend(['Leading Edge', 'Symmetry Plane', 'Base Plane'], loc = 'upper right')
    leg.set_draggable(True)

    # plot upper and lower surface streams in 3D
    ax3d = fig.add_subplot(gs[1,0], projection = '3d')
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