from waverider_generator.generator import waverider
import cadquery as cq
from cadquery import exporters
import numpy as np

def to_CAD(waverider:waverider):

    # extract streams from waverider object
    us_streams=waverider.upper_surface_streams
    ls_streams=waverider.lower_surface_streams

    # compute LE 
    le_right  = np.vstack([x[0] for x in us_streams])
    le_right=le_right
    le_left = np.flip(le_right.copy(), 0)
    le_left[:, 2] = -le_left[:, 2]
    le = np.vstack([le_left[:-1, :], le_right]) 

    # compute TE upper surface
    te_right_upper_surface=np.vstack([x[-1] for x in us_streams])
    te_left_upper_surface=np.flip(te_right_upper_surface.copy(), 0)
    te_left_upper_surface[:,2] = -te_left_upper_surface[:,2]
    te_upper_surface=np.vstack([te_left_upper_surface[:-1,:],te_right_upper_surface])

    # compute TE lower surface
    te_right_lower_surface=np.vstack([x[-1] for x in ls_streams])
    te_left_lower_surface=np.flip(te_right_lower_surface.copy(), 0)
    te_left_lower_surface[:,2] = -te_left_lower_surface[:,2]
    te_lower_surface=np.vstack([te_left_lower_surface[:-1,:],te_right_lower_surface])

    # add interior points for upper surface
    us_points=[]
    for i in range(len(us_streams)):
        for j in range(1, us_streams[i].shape[0]-1):
            us_points.append(tuple(us_streams[i][j]))
    
    # add interior points for lower surface
    ls_points=[]
    for i in range(len(ls_streams)):
        for j in range(1, ls_streams[i].shape[0]-1):
            ls_points.append(tuple(ls_streams[i][j]))

    # create boundaries
    edge_wire_te_upper_surface = cq.Workplane("XY").spline([tuple(x) for x in le])
    edge_wire_te_lower_surface = cq.Workplane("XY").spline([tuple(x) for x in le])
    edge_wire_te_upper_surface = edge_wire_te_upper_surface.add(cq.Workplane("XY").spline([tuple(x) for x in te_upper_surface]))
    edge_wire_te_lower_surface = edge_wire_te_lower_surface.add(cq.Workplane("XY").spline([tuple(x) for x in te_lower_surface]))

    # create upper surface
    upper_surface= cq.Workplane("XY").interpPlate(edge_wire_te_upper_surface, us_points, 0)

    # create lower surface
    lower_surface= cq.Workplane("XY").interpPlate(edge_wire_te_lower_surface, ls_points, 0)

    # add back as a plane
    e1 =cq.Edge.makeSpline([cq.Vector(tuple(x)) for x in te_lower_surface])
    e2=cq.Edge.makeSpline([cq.Vector(tuple(x)) for x in te_upper_surface])
    back = cq.Face.makeFromWires(cq.Wire.assembleEdges([e1, e2]))

    # Create solid
    # waverider_solid =cq.Shell.makeShell([upper_surface.objects[0],lower_surface.objects[0],back])
    waverider_solid = cq.Solid.makeSolid(cq.Shell.makeShell([upper_surface.objects[0], lower_surface.objects[0], back]))
    cq.exporters.export(waverider_solid, f'waverider.step')
    return edge_wire_te_upper_surface
    # return waverider_solid
    