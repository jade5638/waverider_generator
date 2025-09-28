from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
from typing import List, Literal, Union

@dataclass(frozen=True, slots=True) 
class GeometricParameters :

    X1      : float
    X2      : float
    X3      : float
    X4      : float

    w       : float # Half Width [m]
    h       : float # Desired Shockwave Height [m]

    def __init__(self, X1 : float, X2 : float, X3 : float, X4 : float, w : float, h : float) :

        X1 = float(X1)
        X2 = float(X2)
        X3 = float(X3)
        X4 = float(X4)
        w  = float(w)
        h  = float(h)
        
        if h <=0: raise ValueError("height must be a positive number")
        
        if w <=0: raise ValueError("width must be a positive number")

        if (0 <= X1 < 1) and (0 <= X2 <= 1):
            if not ((X2/((1-X1)**4))<(7/64)*(w/h)**4):
                raise ValueError("Condition for inverse design not respected, check value of design parameters X1 and X2")
        else:
            raise ValueError("X1 and/or X2 are not in the required range")
        
        if not (0 <= X3 <= 1) : raise ValueError("X3 must be between 0 and 1")
        if not (0 <= X4 <= 1) : raise ValueError("X4 must be between 0 and 1")

        object.__setattr__(self, 'X1', X1)
        object.__setattr__(self, 'X2', X2)
        object.__setattr__(self, 'X3', X3)
        object.__setattr__(self, 'X4', X4)
        object.__setattr__(self, 'w' , w)
        object.__setattr__(self, 'h' , h)

@dataclass(frozen=True, slots=True)
class FlowParameters :

    M_design    : float # Design Mach Number
    beta        : float # Shock Angle

    def __init__(self, M_design : float, beta : float) :
        
        M_design    = float(M_design)
        beta        = float(beta)

        if M_design <= 0        : raise ValueError('Design Mach Number must be a positive number')
        if not (0 < beta < 90)  : raise ValueError('Shock Angle must be between 0 and 90 degrees (exclusive)')

        object.__setattr__(self, 'M_design' , M_design)
        object.__setattr__(self, 'beta'     , beta)


@dataclass(frozen=True, slots=True)
class OptionalParameters:

    n_USC_pts               : int   = 1000  # Number of Upper Surface Curve points used for Interpolation
    n_SC_pts                : int   = 1000  # Number of Shockwave Curve points used for Interpolation
    n_planes                : int   = 50    # Number of Osculating Planes
    n_US_streamlines_pts    : int   = 50    # Number of points discretising the Upper Surface in the streamwise direction per Osculating Plane
    dx_LS_streamlines_pts   : float = 0.1   # Maximum step to be used in the tracing of the streamlines for the lower surface. Entered as a percentage of the total waverider length 


    def __post_init__(self) :

            # USC Points
            if int(self.n_USC_pts) != self.n_USC_pts :
                raise ValueError("Optional Parameter n_USC_pts must be an integer value")
            
            if self.n_USC_pts < 10 :
                    raise ValueError("Optional Parameter n_USC_pts must be greater or equal to 10")
            
            object.__setattr__(self, 'n_USC_pts', int(self.n_USC_pts))

            # SC Points
            if int(self.n_SC_pts) != self.n_SC_pts :
                raise ValueError("Optional Parameter n_SC_pts must be an integer value")
            
            if self.n_USC_pts < 10 :
                raise ValueError("Optional Parameter n_SC_pts must be greater or equal to 10")
            
            object.__setattr__(self, 'n_SC_pts', int(self.n_SC_pts))

            # Osculating Planes
            if int(self.n_planes) != self.n_planes :
                raise ValueError("Optional Parameter n_planes must be an integer value")
            
            if self.n_planes < 10 :
                raise ValueError("Optional Parameter n_planes must be greater or equal to 10")
            
            object.__setattr__(self, 'n_planes', int(self.n_planes))

            # US Streamline Points
            if int(self.n_US_streamlines_pts) != self.n_US_streamlines_pts :
                raise ValueError("Optional Parameter n_US_streamlines_pts must be an integer value")
            
            if self.n_US_streamlines_pts < 10 :
                raise ValueError("Optional Parameter n_US_streamlines_pts must be greater or equal to 10")
            
            object.__setattr__(self, 'n_US_streamlines_pts', int(self.n_US_streamlines_pts))

            # Step for streamline tracing
            if not (0 < self.dx_LS_streamlines_pts <= 0.2) :
                raise ValueError("Optional Parameter dx_LS_streamlines_pts must be a percentage value less than or equal to 20%")
            
            object.__setattr__(self, 'dx_LS_streamlines_pts', float(self.dx_LS_streamlines_pts))

            




                



        













    



        
