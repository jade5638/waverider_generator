import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from  scipy.integrate import solve_ivp
from waverider_generator.flowfield import cone_angle,cone_field
from typing import Union, Optional
from waverider_generator.input_validation import *
from abc import ABC, abstractmethod
from waverider_generator.curves import *

class IWaverider(ABC) :

    __slots__ = ("_geo_params", "_flow_params", "_opts",
                 "_theta", "_length")
    
    @property
    @abstractmethod
    def geo_params(self) -> GeometricParameters : ...

    @property
    @abstractmethod
    def flow_params(self) -> FlowParameters :...

    @property
    @abstractmethod
    def opts(self) -> OptionalParameters : ...

    @property
    @abstractmethod
    def theta(self) -> float : ...

    @property
    @abstractmethod
    def length(self) -> float : ...

class WaveriderCore(IWaverider) : 

    __slots__ = ()

    def __init__(self, geo_params   : GeometricParameters, 
                       flow_params  : FlowParameters,
                       opts         : Optional[OptionalParameters] = None) :
        
        if not isinstance(geo_params,   GeometricParameters) : raise TypeError("geo_params must be a 'GeometricParameters' instance")
        if not isinstance(flow_params,  FlowParameters)      : raise TypeError("flow_params must be a 'FlowParameters' instance")
        if opts is None:
            opts = OptionalParameters()
        elif not isinstance(opts, OptionalParameters):
            raise TypeError("opts must be an 'OptionalParameters' instance or None")
        
        self._geo_params     = geo_params
        self._flow_params    = flow_params
        self._opts           = opts

        # extract parameters
        beta     = self._flow_params.beta
        M_design = self._flow_params.M_design
        gamma    = self._flow_params.gamma
        h        = self._geo_params.h

        self._theta          = calculate_deflection_angle(beta, M_design, gamma)
        self._length         = h/np.tan(np.radians(beta))
         

    @property
    def geo_params(self) -> GeometricParameters:
        return self._geo_params

    @property
    def flow_params(self) -> FlowParameters:
        return self._flow_params

    @property
    def opts(self) -> OptionalParameters:
        return self._opts
    
    @property
    def theta(self) -> float:
        return self._theta
    
    @property
    def length(self) -> float:
        return self._length
    

class Waverider(WaveriderCore) :

    __slots__ = ()
    
    def Build(self) :

        X1 = self._geo_params.X1
        X2 = self._geo_params.X2
        X3 = self._geo_params.X3
        X4 = self._geo_params.X4
        w = self._geo_params.w
        h = self._geo_params.h
        length = self._length
        
        # SC Control Points in base plane coordinates
        SC_ControlPoints = Curve(x=self.length * np.ones(5),
                                y=np.concatenate((np.zeros(4), np.array([h * X2]))),
                                z=np.linspace(X1 * w, w, 5))
        SC_Bezier = BezierCurve(SC_ControlPoints)
        
        # USC Control Points in base plane coordinates
        USC_ControlPoints = Curve(x=self.length * np.ones(3),
                                  y=np.array([h, h - (1 - X2) * X3 * h, h - (1 - X2) * X4 * h]),
                                  z=np.linspace(0, w, 4)[:3])
        
        USC_ControlPoints.add_point(SC_ControlPoints[-1])
        USC_Bezier = BezierCurve(USC_ControlPoints)

        # z_base contains the z-coordinates to iterate through the 
        # width of the waverider in the base plane
        z_base = np.linspace(0, w, self._opts.n_planes)
        
        # z coordinate after which the SC is defined by the Bezier curve
        z_crit = w * X1
        
        # get the SC base of points and corresponding USC intersections
        SC_base = Curve()
        USC_intersections = Curve()
        for i, z in enumerate(z_base) :
            if z <= z_crit or X2 == 0 :

                SC_point = Point(x = length, y = 0, z = z)
                t_USC = find_t(USC_Bezier, z_target=z)
                USC_point : Point = USC_Bezier.Evaluate(t_USC) # type: ignore

            elif z == w :

                SC_point  = SC_ControlPoints[-1]
                USC_point = SC_ControlPoints[-1]

            else :

                t_SC = find_t(SC_Bezier, z)
                SC_point : Point = SC_Bezier.Evaluate(t_SC) # type: ignore
                
                dCdt : Point = SC_Bezier.Evaluate_FirstDerivative(t_SC) # type: ignore
                dydz = dCdt.y / dCdt.z
                slope = - 1 / dydz
                
                def f(z) :
                    t_USC = find_t(USC_Bezier,z_target=z) 
                    return USC_Bezier.Evaluate(t_USC).y - base_plane_line_equation(z, slope, SC_point)
                
                intersection=root_scalar(f,bracket=[0, w])
                z_USC = intersection.root
                t_USC = find_t(USC_Bezier, z_target=z_USC)
                USC_point : Point = USC_Bezier.Evaluate(t_USC) # type: ignore
            
            SC_base.add_point(SC_point)
            USC_intersections.add_point(USC_point)
        
        # to do : get the leading edge coordinates and cone centres
        


        return SC_base, USC_intersections
    
    def _local_to_global(self, point : Point) -> Point : 
        
        if not isinstance(point, Point) : return TypeError('point must be of type Point')

        return Point(x = point.x,
                     y = point.y - self._geo_params.h,
                     z = point.z)


def base_plane_line_equation(z : float, m : float, known_point : Point) -> float:
    c = known_point.y - m * known_point.z 
    y = m * z + c
    return y

def find_t(bezier_curve : BezierCurve, z_target : float) -> float :

    def f(t):
        return bezier_curve.Evaluate(t).z - z_target

    intersection=root_scalar(f,bracket=[0,1])

    return intersection.root

def calculate_deflection_angle(betaDeg : float, M_inf : float, gamma : float = 1.4) -> float :

    betaRad = np.radians(betaDeg)

    tanTheta = 2*cot(betaRad)*(M_inf**2*np.sin(betaRad)**2-1)/(M_inf**2*(gamma+np.cos(2*betaRad))+2)

    thetaDeg = np.degrees(np.arctan(tanTheta))

    return thetaDeg

# cotangent
def cot(angle : float) -> float:
    return 1/np.tan(angle)
















