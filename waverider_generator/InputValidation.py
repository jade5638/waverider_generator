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









def isInRange(value : Union[float, int], lb : Union[float, int], ub : Union[float, int], 
              bInclusiveLB = True, 
              bInclusiveUB = True) -> bool:

    # lb check
    if bInclusiveLB:
        lower_ok = value >= lb
    else:
        lower_ok = value > lb

    # ub check
    if bInclusiveUB:
        upper_ok = value <= ub
    else:
        upper_ok = value < ub

    return lower_ok and upper_ok
    


        













    



        
