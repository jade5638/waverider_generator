from dataclasses import dataclass, replace
import numpy as np
from typing import List, Literal, Union, Optional, Iterable
from math import factorial as factorial
@dataclass(frozen=True, slots=True)
class Point : 

    x : float
    y : float
    z : float

    def __init__(self, x : float, y : float, z : float) :
        
        x = float(x)
        y = float(y)
        z = float(z)

        object.__setattr__(self, 'x', x)
        object.__setattr__(self, 'y', y)
        object.__setattr__(self, 'z', z)

    def __add__(self, otherPoint) : 

        if not isinstance(otherPoint, Point) :
            return NotImplemented
        return Point(self.x + otherPoint.x, self.y + otherPoint.y, self.z + otherPoint.z)
    
    def __sub__(self, otherPoint):
        if not isinstance(otherPoint, Point):
            return NotImplemented
        return Point(self.x - otherPoint.x, self.y - otherPoint.y, self.z - otherPoint.z)
    
    def __mul__(self, scalar):
        if not (np.isreal(scalar)) :
            return NotImplemented
        scalar = float(scalar)
        return Point(self.x * scalar, self.y * scalar, self.z * scalar)
    
    def __rmul__(self, scalar) :
        return self.__mul__(scalar)
    
    def __truediv__(self, scalar):
        if not (np.isreal(scalar)) :
            return NotImplemented
        scalar = float(scalar)
        return Point(self.x / scalar, self.y / scalar, self.z / scalar)

    def __rtruediv__(self, scalar):
        if not (np.isreal(scalar)) :
            return NotImplemented
        scalar = float(scalar)
        return Point(scalar / self.x, scalar / self.y, scalar / self.z)
    
    def distanceTo(self, otherPoint) -> float:
        dx, dy, dz = self.x - otherPoint.x, self.y - otherPoint.y, self.z - otherPoint.z
        return np.sqrt(dx*dx + dy*dy + dz*dz)

class Curve:
    def __init__(self,
                 points: Optional[Iterable[Point]] = None,
                 x: Optional[np.ndarray] = None,
                 y: Optional[np.ndarray] = None,
                 z: Optional[np.ndarray] = None):
        """
        create curve from:
          - a sequence of Point objects
          - or three np arrays x, y, z
        """

        # from points
        if points is not None:
            if not isinstance(points, (list, tuple, np.ndarray)):
                raise TypeError("points must be a list, tuple, or numpy ndarray of Point objects.")
            for p in points:
                if not isinstance(p, Point):
                    raise TypeError(f"All elements in points must be Point objects.")
            self.points: list[Point] = list(points)
            return

        # from x, y z vectors
        if x is not None or y is not None or z is not None:
            # check if all three provided
            if x is None or y is None or z is None:
                raise ValueError("If providing coordinates, x, y, and z must all be provided.")
            
            # check if numpy vectors
            for name, arr in zip(('x', 'y', 'z'), (x, y, z)):
                if not isinstance(arr, (list, tuple, np.ndarray)):
                    raise TypeError(f"{name} must be a list, tuple, or numpy ndarray.")

            # check length
            if not (x.shape == y.shape == z.shape):
                raise ValueError("x, y, and z arrays must have the same shape.")

            # make into a list of points
            self.points: list[Point] = [Point(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
            return

        # empty curve
        self.points: list[Point] = []


    def add_point(self, point: Point) -> None:
        '''add point from a Point object '''
        if not isinstance(point, Point):
            raise TypeError(f"Can only add Point objects.")
        self.points.append(point)

    
    def add_coords(self, x: float, y: float, z: float) -> None:
        '''add point from set of x, y, z vectors '''
        self.add_point(Point(x, y, z))

    def __len__(self):
        return len(self.points)

    def __getitem__(self, index):
        return self.points[index]

    def __iter__(self):
        return iter(self.points)

    def __repr__(self):
        return f"Curve({len(self.points)} points)"

    @property
    def x(self) -> np.ndarray:
        return np.array([p.x for p in self.points])

    @property
    def y(self) -> np.ndarray:
        return np.array([p.y for p in self.points])

    @property
    def z(self) -> np.ndarray:
        return np.array([p.z for p in self.points])


class BezierCurve():

    def __init__(self, controlPoints : Curve) :
        
        if not isinstance(controlPoints, Curve) : 
            return TypeError('Control points must be a Curve instance')
        
        self.controlPoints = controlPoints
        self.nOrder = len(controlPoints) - 1

    def Evaluate(self, t : Union[float,Iterable]) -> Union[Curve, Point]:
        
        # handle scalar input
        if isinstance(t, (int, float)):
            tVec = [float(t)]
        else:
            tVec = list(np.sort(np.array(t)))

        curve = Curve()
        for t in tVec :
            if not (0 <= t <= 1) : # type: ignore
                raise ValueError('t must be between zero and one')
            point = sum([BernsteinPolynomial(self.nOrder, i, t) * self.controlPoints[i] for i in range(self.nOrder + 1)], 
                    Point(0, 0, 0))
            curve.add_point(point)

        if len(curve) > 1 :
            return curve
        else :
            return curve[0]
    
    def Evaluate_FirstDerivative(self, t : Union[float,Iterable]) -> Union[Curve, Point]:

        n = self.nOrder

        # handle scalar input
        if isinstance(t, (int, float)):
            tVec = [float(t)]
        else:
           tVec = list(np.sort(np.array(t)))

        curve = Curve()
        for t in tVec :
            if not (0 <= t <= 1) : # type: ignore
                raise ValueError('t must be between zero and one')
            c_prime = sum([
                BernsteinPolynomial(n - 1, i, t)
                * n 
                * (self.controlPoints[i + 1] - self.controlPoints[i]) for i in range(n)
                ],
                Point(0, 0, 0))
            curve.add_point(c_prime)
        
        if len(curve) > 1 :
            return curve
        else :
            return curve[0]
    
    def Evaluate_SecondDerivative(self, t : Union[float,Iterable]) -> Union[Curve, Point]:
        
        n = self.nOrder
        
        # handle scalar input
        if isinstance(t, (int, float)):
            tVec = [float(t)]
        else:
           tVec = list(np.sort(np.array(t)))

        curve = Curve()
        for t in tVec :
            if not (0 <= t <= 1) : # type: ignore
                raise ValueError('t must be between zero and one')
            c_prime_prime = sum([
                BernsteinPolynomial(n - 2, i, t)
                * n * (n - 1)
                * (self.controlPoints[i + 2] - 2 * self.controlPoints[i + 1] + self.controlPoints[i])
                for i in range(n - 1)
            ],
            Point(0, 0, 0))
            curve.add_point(c_prime_prime)
        
        if len(curve) > 1 :
            return curve
        else :
            return curve[0]
    
    
 
def BernsteinPolynomial(n, i, t) -> float:

    return float((factorial(n) / (factorial(i) * factorial(n - i))) * (t ** i) * ((1 - t) ** (n - i))) 

    
        
        
