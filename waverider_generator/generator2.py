import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from  scipy.integrate import solve_ivp
from waverider_generator.flowfield import cone_angle,cone_field
from typing import Union, Optional
from waverider_generator.InputValidation import *
from abc import ABC, abstractmethod


class IWaverider(ABC) :

    __slots__ = ("_geo_params", "_flow_params", "_opts")
    
    @property
    @abstractmethod
    def geo_params(self) -> GeometricParameters : ...
    @property
    @abstractmethod
    def flow_params(self) -> FlowParameters :...
    @property
    @abstractmethod
    def opts(self) -> OptionalParameters : ...


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

    @property
    def geo_params(self):
        return self._geo_params

    @property
    def flow_params(self):
        return self._flow_params

    @property
    def opts(self) -> OptionalParameters:
        return self._opts


class Waverider(WaveriderCore) :

    __slots__ = ()
    
    pass














