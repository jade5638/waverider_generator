from dataclasses import dataclass, replace
import numpy as np
from typing import List, Literal, Union, Optional, Iterable, overload
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
    
    def __eq__(self, otherPoint):
        if not isinstance(otherPoint, Point):
            return NotImplemented
        return self.x == otherPoint.x and self.y == otherPoint.y and self.z == otherPoint.z
    
    def distanceTo(self, otherPoint) -> float:
        dx, dy, dz = self.x - otherPoint.x, self.y - otherPoint.y, self.z - otherPoint.z
        return np.sqrt(dx*dx + dy*dy + dz*dz)
    
    def toTuple(self) -> tuple :
        return (self.x, self.y, self.z)

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


# --------------------------
# Bézier Curve Class
# --------------------------
class BezierCurve:

    @overload
    def Evaluate(self, t: float) -> Point: ...

    @overload
    def Evaluate(self, t: Iterable) -> Curve: ...

    @overload
    def Evaluate_FirstDerivative(self, t: float) -> Point: ...

    @overload
    def Evaluate_FirstDerivative(self, t: Iterable) -> Curve: ...

    @overload
    def Evaluate_SecondDerivative(self, t: float) -> Point: ...

    @overload
    def Evaluate_SecondDerivative(self, t: Iterable) -> Curve: ...

    def __init__(self, controlPoints):


        if len(controlPoints) < 2:
            raise ValueError("Bezier curve requires at least 2 control points.")

        self.controlPoints = controlPoints
        self.nOrder = len(controlPoints) - 1

    def Evaluate(self, t: Union[float, Iterable]) -> Union[Point, Curve]:

        '''
        Evaluate Bezier curve at t
        '''

        tVec = [float(t)] if isinstance(t, (int, float)) else [float(x) for x in t]

        # initialise Curve
        curve = Curve()

        n = self.nOrder

        for ti in tVec:
            if not (0 <= ti <= 1):
                raise ValueError("t must be in [0, 1].")

            p = Point(0, 0, 0)

            # sum
            for i in range(n + 1):
                p += BernsteinPolynomial(n, i, ti) * self.controlPoints[i]

            # add point to curve
            curve.add_point(p)

        if len(curve) == 1 :
            point : Point = curve[0]
            return point 
        else :
            return curve

    def Evaluate_FirstDerivative(self, t: Union[float, Iterable]) -> Union[Point, Curve]:

        tVec = [float(t)] if isinstance(t, (int, float)) else [float(x) for x in t]

        # initialise Curve
        curve = Curve()
        
        n = self.nOrder

        for ti in tVec:
            if not (0 <= ti <= 1):
                raise ValueError("t must be in [0, 1].")

            p = Point(0, 0, 0)

            # sum
            for i in range(n):
                p += n * BernsteinPolynomial(n - 1, i, ti) * (self.controlPoints[i + 1] - self.controlPoints[i])

            # add point to curve
            curve.add_point(p)

        if len(curve) == 1 :
            point : Point = curve[0]
            return point 
        else :
            return curve

    def Evaluate_SecondDerivative(self, t: Union[float, Iterable]) -> Union[Point, Curve]:

        n = self.nOrder
        if n < 2:
            raise ValueError("Second derivative requires at least three control points (Bezier curve of order > 2).")

        tVec = [float(t)] if isinstance(t, (int, float)) else [float(x) for x in t]

        # initialise Curve
        curve = Curve()

        for ti in tVec:
            if not (0 <= ti <= 1):
                raise ValueError("t must be in [0, 1].")

            p = Point(0, 0, 0)

            # sum
            for i in range(n - 1):
                p += n*(n - 1) * BernsteinPolynomial(n - 2, i, ti) * (
                    self.controlPoints[i + 2]
                    - 2 * self.controlPoints[i + 1]
                    + self.controlPoints[i]
                )

            # add point to curve
            curve.add_point(p)

        if len(curve) == 1 :
            point : Point = curve[0]
            return point 
        else :
            return curve
    
def BernsteinPolynomial(n: int, i: int, t: float) -> float:
    """
    Bernstein polynomial B_{i,n}(t).
    """
    if i < 0 or i > n:
        return 0.0
    
    return (factorial(n) / (factorial(i) * factorial(n - i))) * (t**i) * ((1 - t)**(n - i))

    
        
        
