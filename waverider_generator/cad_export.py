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
    '''
    le_right=le_right
    le_left = np.flip(le_right.copy(), 0)
    le_left[:, 2] = -le_left[:, 2]
    le = np.vstack([le_left[:-1, :], le_right]) 
    '''
    le=le_right
    # compute TE upper surface
    te_right_upper_surface=np.vstack([x[-1] for x in us_streams])
    '''
    te_left_upper_surface=np.flip(te_right_upper_surface.copy(), 0)
    te_left_upper_surface[:,2] = -te_left_upper_surface[:,2]
    te_upper_surface=np.vstack([te_left_upper_surface[:-1,:],te_right_upper_surface])
    '''
    te_upper_surface=te_right_upper_surface
    # compute TE lower surface
    te_right_lower_surface=np.vstack([x[-1] for x in ls_streams])
    '''
    te_left_lower_surface=np.flip(te_right_lower_surface.copy(), 0)
    te_left_lower_surface[:,2] = -te_left_lower_surface[:,2]
    te_lower_surface=np.vstack([te_left_lower_surface[:-1,:],te_right_lower_surface])
    '''
    te_lower_surface=te_right_lower_surface
    # add interior points for upper surface
    us_points=[]
    for i in range(len(us_streams)):
        for j in range(1, us_streams[i].shape[0]-1):
            us_points.append(tuple(us_streams[i][j]))

    # for i in range(1, len(us_streams)):
    #     for j in range(1, us_streams[i].shape[0]-1):
    #         pt = us_streams[i][j]
    #         pt[2] = -pt[2]
    #         us_points.append(tuple(pt))
    
    # add interior points for lower surface
    ls_points=[]
    for i in range(len(ls_streams)):
        for j in range(1, ls_streams[i].shape[0]-1):
            ls_points.append(tuple(ls_streams[i][j]))

    # for i in range(1, len(ls_streams)):
    #     for j in range(1, ls_streams[i].shape[0]-1):
    #         pt = ls_streams[i][j]
    #         pt[2] = -pt[2]
    #         ls_points.append(tuple(pt))

    # create boundaries
    # Define points
    points_upper_surface = [(0, 0, 0), (waverider.length, 0, 0)]
    points_lower_surface=[(0,0,0),(waverider.length,ls_streams[0][-1,1],0)]
    # Create a workplane and draw lines between points
    workplane = cq.Workplane("XY")
    edge_wire_te_upper_surface = workplane.moveTo(points_upper_surface[0][0], points_upper_surface[0][1])
    edge_wire_te_lower_surface=workplane.moveTo(points_lower_surface[0][0], points_lower_surface[0][1])

    for point in points_upper_surface[1:]:
        edge_wire_te_upper_surface = edge_wire_te_upper_surface.lineTo(point[0], point[1])
    for point in points_lower_surface[1:]:
        edge_wire_te_lower_surface = edge_wire_te_lower_surface.lineTo(point[0], point[1])

    edge_wire_te_upper_surface = edge_wire_te_upper_surface.add(cq.Workplane("XY").spline([tuple(x) for x in le]))
    edge_wire_te_lower_surface = edge_wire_te_lower_surface.add(cq.Workplane("XY").spline([tuple(x) for x in le]))
    edge_wire_te_upper_surface = edge_wire_te_upper_surface.add(cq.Workplane("XY").spline([tuple(x) for x in te_upper_surface]))
    edge_wire_te_lower_surface = edge_wire_te_lower_surface.add(cq.Workplane("XY").spline([tuple(x) for x in te_lower_surface]))
    
    # create upper surface
    upper_surface= cq.Workplane("XY").interpPlate(edge_wire_te_upper_surface, us_points, 0)

    # # create lower surface
    lower_surface= cq.Workplane("XY").interpPlate(edge_wire_te_lower_surface, ls_points, 0)
    
    # add back as a plane
    e1 =cq.Edge.makeSpline([cq.Vector(tuple(x)) for x in te_lower_surface])
    e2=cq.Edge.makeSpline([cq.Vector(tuple(x)) for x in te_upper_surface])
    sym_edge=np.vstack(((waverider.length, 0, 0),(waverider.length,ls_streams[0][-1,1],0)))
    v1 = cq.Vector(*sym_edge[0])
    v2 = cq.Vector(*sym_edge[1])
    e3 = cq.Edge.makeLine(v1, v2)
    back = cq.Face.makeFromWires(cq.Wire.assembleEdges([e1, e2,e3]))
    
    # Create solid
    # waverider_solid =cq.Shell.makeShell([upper_surface.objects[0],lower_surface.objects[0],back])
    right_side = cq.Solid.makeSolid(cq.Shell.makeShell([upper_surface.objects[0], lower_surface.objects[0], back]))
    # left_side= right_side.mirror(mirrorPlane='XY')
    # # waverider_solid = cq.Workplane("XY").add(right_side).add(left_side).combineSolids()
    # waverider_solid = (
    # cq.Workplane("XY")
    # .newObject([right_side])
    # .union(left_side)
    # )

    waverider_solid=right_side
    # cq.exporters.export(waverider_solid, f'waverider.step')
    
    return waverider_solid
    