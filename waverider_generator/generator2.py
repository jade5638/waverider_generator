import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from  scipy.integrate import solve_ivp
from waverider_generator.flowfield import cone_angle,cone_field
from typing import Union, Optional
from waverider_generator.input_validation import *
from abc import ABC, abstractmethod
from waverider_generator.curves import *
import math

class IWaverider(ABC) :

    __slots__ = ("_geo_params", "_flow_params")
    
    @property
    @abstractmethod
    def geo_params(self) -> GeometricParameters : ...

    @property
    @abstractmethod
    def flow_params(self) -> FlowParameters : ...

class WaveriderCore(IWaverider) : 

    __slots__ = ()

    def __init__(self, geo_params   : GeometricParameters, 
                       flow_params  : FlowParameters) :
        
        if not isinstance(geo_params,   GeometricParameters) : raise TypeError("geo_params must be a 'GeometricParameters' instance")
        if not isinstance(flow_params,  FlowParameters)      : raise TypeError("flow_params must be a 'FlowParameters' instance")

        self._geo_params     = geo_params
        self._flow_params    = flow_params

    @property
    def geo_params(self) -> GeometricParameters:
        return self._geo_params

    @property
    def flow_params(self) -> FlowParameters:
        return self._flow_params
    
class Waverider(WaveriderCore) :

    __slots__ = ()

    @property
    def thetaDeg(self) -> float :
        return calculate_deflection_angle(self._flow_params)
    
    @property
    def thetaRad(self) -> float :
        return np.radians(self.thetaDeg)
    
    @property
    def length(self) -> float :
        return self.h / np.tan(self.betaRad)
    
    @property
    def cone_angle_deg(self) -> float :
        return cone_angle(self.M_design, self.betaDeg, self.gamma)
    
    @property
    def cone_angle_rad(self) -> float :
        return np.radians(self.cone_angle_deg)

    @property
    def M_design(self) -> float :
        return self.flow_params.M_design
    
    @property
    def betaDeg(self) -> float :
        return self.flow_params.betaDeg
    
    @property
    def betaRad(self) -> float : 
        return np.radians(self.betaDeg)
    
    @property
    def gamma(self) -> float :
        return self.flow_params.gamma
    
    @property
    def X1(self) -> float :
        return self.geo_params.X1
    
    @property
    def X2(self) -> float :
        return self.geo_params.X2
    
    @property
    def X3(self) -> float :
        return self.geo_params.X3
    
    @property
    def X4(self) -> float :
        return self.geo_params.X4
    
    @property
    def w(self) -> float :
        return self.geo_params.w
    
    @property
    def h(self) -> float :
        return self.geo_params.h
    
    @property
    def SC_ControlPoints(self) -> Curve :
        return Curve(
            x = self.length * np.ones(5),
            y = np.concatenate((np.zeros(4), np.array([self.h * self.X2]))) - self.h,
            z = np.linspace(self.X1 * self.w, self.w, 5)
            )
    
    @property
    def USC_ControlPoints(self) -> Curve :
        curve = Curve(
            x = self.length * np.ones(3),
            y = np.array([self.h, self.h - (1 - self.X2) * self.X3 * self.h, self.h - (1 - self.X2) * self.X4 * self.h]) - self.h,
            z = np.linspace(0, self.w, 4)[:3]
            )
        
        curve.add_point(self.SC_ControlPoints[-1])

        return curve
    
    @property
    def SC_Bezier(self) -> BezierCurve :
        return BezierCurve(self.SC_ControlPoints)

    @property
    def USC_Bezier(self) -> BezierCurve :
        return BezierCurve(self.USC_ControlPoints)
    
    @property
    def z_crit(self) -> float :
        return self.w * self.X1
    
    @property
    def lateral_tip_point(self) -> Point :
        return self.SC_ControlPoints[-1]

    def Build(self, 
            z_base                  : Union[Iterable, None] = None,
            n_planes                : int   = 50,
            n_streamline_pts        : int   = 50,  
            dx_LS_streamline        : float = 0.1) :
        
        # n_planes input validation 
        if int(n_planes) != n_planes : raise ValueError("n_planes must be an integer")
        if n_planes < 3 : raise ValueError("n_planes must be >= 3")

        # n_streamline_pts input validation
        if n_streamline_pts < 10 : raise ValueError("n_streamline_pts must be >= 10")

        if not (0 < dx_LS_streamline <= 0.2) : raise ValueError("dx_LS_streamlines_pts must be in (0, 0.2]")
        
        if z_base is None : z_base = self.w * exp_spacing(n_planes) 
        z_base = np.sort(np.array(z_base))
        if not (z_base[0]  == 0)         : raise ValueError('Minimum value in z_base must be 0')
        if not (z_base[-1] == self.w)    : raise ValueError('Maximum value in z_base must be the half-width of the waverider "w"')

        SC_Bezier           = self.SC_Bezier
        USC_Bezier          = self.USC_Bezier
        
        # initialise, all Curve instances in global coordinate system
        SC_base             = Curve()
        USC_intersections   = Curve()
        cone_centers        = Curve()
        leading_edge        = Curve()
        radii_curvature     = []
        dy_dz_slope         = []

        for z in z_base :

            if not (0 <= z <= self.w): raise ValueError('z must be in [0, w]')

            # case 1 : flat part of the SC
            if z <= self.z_crit or self.X2 == 0 :

                # get the SC point
                SC_point = Point(x = self.length, y = - self.h, z = z)

                # get the t_USC and get the USC point
                t_USC = find_t(USC_Bezier, z_target=z)
                USC_point = USC_Bezier.Evaluate(t_USC) 

                # get the 'cone' center point
                cone_point = Point(x = self.length - (USC_point.y - SC_point.y) / np.tan(self.betaRad), 
                                   y = USC_point.y,
                                   z = z)
                 
                # get the leading edge point
                # in the flat part this is equivalent to the 'cone center' 
                leading_edge_point = cone_point

                # infinite radius of curvature for the flat part
                radius_curvature = math.inf
                
                # slope undefined in flat part
                dy_dz = math.inf

            # case 2 : lateral tip of the waverider
            elif z == self.w :

                SC_point            = self.lateral_tip_point
                USC_point           = self.lateral_tip_point
                cone_point          = self.lateral_tip_point
                leading_edge_point  = self.lateral_tip_point
                radius_curvature    = calculate_radius_curvature(dC_dt= SC_Bezier.Evaluate_FirstDerivative(1), 
                                                                 dC2_dt2=SC_Bezier.Evaluate_SecondDerivative(1))
                
                dSC_dt = SC_Bezier.Evaluate_FirstDerivative(1) 
                dy_dz = dSC_dt.y / dSC_dt.z

            # case 3 : curved section of the SC
            else :

                # get t_SC and get the SC point
                t_SC = find_t(SC_Bezier, z)
                SC_point = SC_Bezier.Evaluate(t_SC)
                
                # get first derivative of SC bezier curve at t_SC
                dSC_dt = SC_Bezier.Evaluate_FirstDerivative(t_SC) 

                # calculate dydz
                dy_dz = dSC_dt.y / dSC_dt.z

                # calculate slope of the line going through SC_point
                # and intersecting with the USC
                slope = - 1 / dy_dz
                
                # get the intersection with the USC
                def f(z) :
                    t_USC = find_t(USC_Bezier,z_target=z) 
                    return USC_Bezier.Evaluate(t_USC).y - base_plane_line_equation(z, slope, SC_point)
                
                intersection=root_scalar(f,bracket=[0, self.w])
                z_USC = intersection.root
                t_USC = find_t(USC_Bezier, z_target=z_USC)
                USC_point = USC_Bezier.Evaluate(t_USC) 

                # calculate second derivative to get radius of curvature
                dSC2_dt2 = SC_Bezier.Evaluate_SecondDerivative(t_SC)  

                radius_curvature = calculate_radius_curvature(dSC_dt, dSC2_dt2)
                alpha = np.arctan(dy_dz)
                z_cone = SC_point.z - np.sin(alpha) * radius_curvature
                y_cone = SC_point.y + np.cos(alpha) * radius_curvature
                # distance between SC Point and cone center projection on base plane
                d = SC_point.distanceTo(Point(x = self.length, y = y_cone, z=z_cone))
                x_cone = self.length - d / np.tan(self.betaRad)

                cone_point = Point(x = x_cone,
                                   y = y_cone,
                                   z = z_cone)
                
                leading_edge_point = calculate_LE_point(cone_point, SC_point, USC_point.y)

            SC_base.add_point(SC_point)
            USC_intersections.add_point(USC_point)
            cone_centers.add_point(cone_point)
            leading_edge.add_point(leading_edge_point)
            radii_curvature.append(radius_curvature)
            dy_dz_slope.append(dy_dz)

        # calculate upper surface streams
        US_Streams = []
        for leading_edge_point, USC_point in zip(leading_edge, USC_intersections):

            US_Streams.append(Curve(x = np.linspace(leading_edge_point.x, USC_point.x, n_streamline_pts),
                                    y = np.full(n_streamline_pts, leading_edge_point.y),
                                    z = np.full(n_streamline_pts, leading_edge_point.z)))
            
        # trace streamlines

        # propagate the streamlines
        Vr, Vt = cone_field(self.M_design, self.cone_angle_rad, self.betaRad, self.gamma)

        # ODE which propagates the streamlines
        def stode(t, x, y_max):

            th = np.arctan(x[1] / x[0])
            
            dxdt = np.zeros(2)
            
            dxdt[0] = Vr(th) * np.cos(th) - np.sin(th) * Vt(th)
            dxdt[1] = Vr(th) * np.sin(th) + np.cos(th) * Vt(th)
            
            return dxdt
        
        def back(t, y, y_max):
            return y[0] - y_max
        
        back.terminal = True

        LS_Streams = []
        for (z,
            cone_point, 
            leading_edge_point,
            USC_point, 
            SC_point, 
            radius_curvature,
            dy_dz) in zip(z_base, 
                                    cone_centers, 
                                    leading_edge, 
                                    USC_intersections,
                                    SC_base,
                                    radii_curvature, 
                                    dy_dz_slope):

            # case 1 : flat part of the SC
            if z <= self.z_crit or self.X2 == 0 :

                # trigonometry with deflection angle
                bottom_surface_y= leading_edge_point.y - np.tan(self.thetaRad) * (self.length - leading_edge_point.x)

                # store the x,y and z in a streams
                x_LS =np.linspace(leading_edge_point.x, self.length, n_streamline_pts)
                y_LS =np.linspace(leading_edge_point.y, bottom_surface_y, n_streamline_pts)
                z_LS =np.full(n_streamline_pts, leading_edge_point.z)

                LS_Streams.append(Curve(x = x_LS,
                                        y = y_LS,
                                        z = z_LS))

            # case 2 : lateral tip of the waverider
            elif z == self.w :

                LS_Streams.append(Curve([self.lateral_tip_point, self.lateral_tip_point]))

            # case 3 : curved section of the SC
            else :

                # need calculate R minus height of osculating plane
                eta_le = radius_curvature - SC_point.distanceTo(USC_point)

                # calculate the angle to rotate the streamlines by
                alpha=np.arctan(dy_dz)

                x_le=(eta_le) / np.tan(self.betaRad) 

                sol = solve_ivp(stode, (0, 1000), [x_le, eta_le], events=back, args=(radius_curvature / np.tan(self.betaRad),), max_step=dx_LS_streamline*self.length)
                stream = np.vstack([sol.y[0], -sol.y[1] * np.cos(alpha), sol.y[1] * np.sin(alpha)]).T

                # transform from cone center coordinate system to global
                x_LS = stream[:,0] + cone_point.x
                y_LS = stream[:,1] + cone_point.y
                z_LS = stream[:,2] + cone_point.z

                LS_Streams.append(Curve(x = x_LS,
                                        y = y_LS,
                                        z = z_LS))
                                        



        

        
        # return SC_base, USC_intersections, cone_centers, leading_edge, radii
        return US_Streams, LS_Streams, {
                "SC_base": SC_base,
                "USC_intersections": USC_intersections,
                "cone_centers" : cone_centers,
                "leading_edge" : leading_edge,
                "radii_curvature" : np.array(radii_curvature),
                "dy_dz" : dy_dz_slope
               }
        

def exp_spacing(n : int, k : int = 2) -> np.ndarray:
    u = np.linspace(0, 1, n)
    return 1 - (np.exp(k*u) - 1) / (np.exp(k) - 1)

def calculate_LE_point(cone_point : Point, SC_point : Point, y_USC : float):

    k=(y_USC-SC_point.y)/(cone_point.y-SC_point.y)

    x_I=SC_point.x + k*(cone_point.x-SC_point.x)
    y_I=y_USC
    z_I=SC_point.z + k*(cone_point.z-SC_point.z)

    return Point(x = x_I, y = y_I, z = z_I)

def calculate_radius_curvature(dC_dt : Point, dC2_dt2 : Point) -> float:

    dz_dt = dC_dt.z
    dz2_dt2 = dC2_dt2.z

    dy_dt = dC_dt.y
    dy2_dt2 = dC2_dt2.y

    radius= 1/(abs((dz_dt*dy2_dt2-dy_dt*dz2_dt2))/((dz_dt**2+dy_dt**2)**(3/2)))

    return radius

def base_plane_line_equation(z : float, m : float, known_point : Point) -> float:
    c = known_point.y - m * known_point.z 
    y = m * z + c
    return y

def find_t(bezier_curve : BezierCurve, z_target : float) -> float :

    def f(t):
        return bezier_curve.Evaluate(t).z - z_target

    intersection=root_scalar(f,bracket=[0,1])

    return intersection.root

def calculate_deflection_angle(flow_params : FlowParameters) -> float :

    betaRad = np.radians(flow_params.betaDeg)
    M_inf   = flow_params.M_design
    gamma   = flow_params.gamma

    tanTheta = 2*cot(betaRad)*(M_inf**2*np.sin(betaRad)**2-1)/(M_inf**2*(gamma+np.cos(2*betaRad))+2)

    thetaDeg = np.degrees(np.arctan(tanTheta))

    return thetaDeg

# cotangent
def cot(angle : float) -> float:
    return 1/np.tan(angle)
















