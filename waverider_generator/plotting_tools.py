import matplotlib.pyplot as plt
import numpy as np
from waverider_generator.generator import waverider


def Plot_Base_Plane(waverider: waverider,latex: bool):
    
    if latex==True:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    elif latex==False:
        pass
    else: 
        raise ValueError("input 'latex' must be a boolean")

    fig, ax=plt.subplots()

    # intersections with upper surface
    inters=waverider.local_intersections_us
    inters=np.vstack([np.array([0,waverider.height]),inters,waverider.us_P3])

    # shockwave points
    shockwave=np.column_stack([waverider.z_local_shockwave,waverider.y_local_shockwave])
    shockwave=np.vstack([np.array([0,0]),shockwave,waverider.s_P4])

    # lower surface at base plane
    lower_surface=waverider.lower_surface_streams
    lower_surface = np.vstack([stream[-1,:] for stream in lower_surface])
    z_ls=lower_surface[:,2]
    y_ls = lower_surface[:,1]+waverider.height

    ax.plot([0, 0], [0, waverider.height], 'b-')
    # osculating planes
    for i, (point1, point2) in enumerate(zip(inters, shockwave)):
        x_values = [point1[0], point2[0]]
        z_values = [point1[1], point2[1]]
        label = 'Osculating Planes' if i == 0 else None  # label only the first line for the legend
        ax.plot(x_values, z_values, 'b-', label=label)

    # lower surface, upper surface and shockwave
    ax.plot(shockwave[:, 0], shockwave[:, 1], 'go--',label="Shockwave")
    ax.plot(inters[:, 0], inters[:, 1], 'r-o',label="Upper Surface")
    ax.plot(z_ls, y_ls, '-ok',label="Lower Surface")

    
    
    X2 = waverider.us_P3
    ax.plot(X2[0], X2[1], 'bo',label="Tip")

    
    # set labels and title
    ax.set_xlabel('z')
    ax.set_ylabel('y')
    ax.set_title(f'Base Plane [X1,X2,X3,X4]=[{waverider.X1},{waverider.X2},{waverider.X3},{waverider.X4}]')
    ax.set_aspect('equal')

    ax.legend()
    # return the figure object
    return fig

def Plot_Leading_Edge(waverider: waverider,latex: bool):

    if latex==True:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
    elif latex==False:
        pass
    else: 
        raise ValueError("input 'latex' must be a boolean")
    fig, ax=plt.subplots()

    le=waverider.leading_edge
    base_point=[5,0]
    ax.plot(le[:,2],le[:,0],'b-',label='Leading Edge')

    ax.plot([le[0,2],base_point[1]],[le[0,0],base_point[0]],'--k',label='Symmetry Plane')
    ax.plot([le[-1,2],base_point[1]],[le[-1,0],base_point[0]],'--r',label='Base Plane')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f"Leading Edge Shape [X1,X2,X3,X4]=[{waverider.X1},{waverider.X2},{waverider.X3},{waverider.X4}] (Top View)")
    ax.set_xlabel('z')
    ax.set_ylabel('x')
    ax.legend()
    return fig