from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
from typing import List, Literal, Union

@dataclass(frozen=True, slots=True) 
class GeometricParameters :
    """
    Data class for the Geometric Parameters of the Waverider. \n
    For reference, see [GitHub README](https://github.com/jade5638/waverider_generator).\n
    These parameters are best visualised in the Base Plane.

    Properties
    ----------
    X1 :
        Parameter determining the length of the flat part of the Shockwave Curve.

    X2 :
        Parameter determining the height of the lateral tip of the Waverider
        relative to the height of the Shockwave Curve.

    X3 :
        Parameter determining the height of the first control point of the Upper
        Surface Curve.

    X4 :
        Parameter determining the height of the second control point of the Upper
        Surface Curve.

    w  :
        Waverider half width, in meters.

    h  :
        Desired Shockwave height, in meters.
    ----------

    """

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
    """
    Data class for the Flow Parameters of the Waverider.\n
    For reference, see [GitHub README](https://github.com/jade5638/waverider_generator).\n

    Properties 
    ----------
    M_design :
        Waverider Design Mach Number.

    betaDeg :
        Target Shock Angle at the Symmetry Plane, in degrees.

    gamma :
        Ratio of Specific Heats (default 1.4 for air).
    ----------

    """

    M_design    : float
    betaDeg     : float
    gamma       : float 

    def __init__(self, M_design : float, betaDeg : float, gamma : float = 1.4) :
        
        M_design    = float(M_design)
        betaDeg     = float(betaDeg)
        gamma       = float(gamma)

        if M_design <= 0        : raise ValueError('Design Mach Number must be a positive number')
        if not (0 < betaDeg < 90)  : raise ValueError('Shock Angle must be between 0 and 90 degrees (exclusive)')
        if not gamma > 1        : raise ValueError('Ratio of Specific Heats must be greater than 1')

        object.__setattr__(self, 'M_design' , M_design)
        object.__setattr__(self, 'betaDeg'  , betaDeg)
        object.__setattr__(self, 'gamma'    , gamma)