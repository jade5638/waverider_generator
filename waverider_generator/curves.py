from dataclasses import dataclass, replace
import numpy as np
from typing import List, Literal, Union, Optional, Iterable

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
    
    def distanceTo(self, otherPoint) -> float:
        dx, dy, dz = self.x - otherPoint.x, self.y - otherPoint.y, self.z - otherPoint.z
        return np.sqrt(dx*dx + dy*dy + dz*dz)
    
import numpy as np
from typing import Optional, Iterable
from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class Point:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float, z: float):
        object.__setattr__(self, 'x', float(x))
        object.__setattr__(self, 'y', float(y))
        object.__setattr__(self, 'z', float(z))


class Curve:
    def __init__(self,
                 points: Optional[Iterable[Point]] = None,
                 x: Optional[np.ndarray] = None,
                 y: Optional[np.ndarray] = None,
                 z: Optional[np.ndarray] = None):
        """
        Create a Curve either from:
          - a sequence of Point objects, OR
          - three NumPy arrays x, y, z of equal length.
        """

        # Case 1: Explicit list of Point objects
        if points is not None:
            if not isinstance(points, (list, tuple, np.ndarray)):
                raise TypeError("points must be a list, tuple, or numpy array of Point objects.")
            for p in points:
                if not isinstance(p, Point):
                    raise TypeError(f"All elements in points must be Point objects, got {type(p)}.")
            self.points: list[Point] = list(points)
            return

        # Case 2: NumPy coordinate vectors
        if x is not None or y is not None or z is not None:
            # Require all three
            if x is None or y is None or z is None:
                raise ValueError("If providing coordinates, x, y, and z must all be provided.")
            
            # Require numpy arrays
            for name, arr in zip(('x', 'y', 'z'), (x, y, z)):
                if not isinstance(arr, np.ndarray):
                    raise TypeError(f"{name} must be a NumPy ndarray, got {type(arr)}")

            # Require matching shape
            if not (x.shape == y.shape == z.shape):
                raise ValueError("x, y, and z arrays must have the same shape.")

            self.points: list[Point] = [Point(xi, yi, zi) for xi, yi, zi in zip(x, y, z)]
            return

        # Case 3: Empty curve
        self.points: list[Point] = []

    def add_point(self, point: Point) -> None:
        """Add a new Point to the end of the curve."""
        if not isinstance(point, Point):
            raise TypeError(f"Can only add Point objects, got {type(point)}.")
        self.points.append(point)

    def add_coords(self, x: float, y: float, z: float) -> None:
        """Create and add a Point from numeric coordinates."""
        self.add_point(Point(x, y, z))

    # ---- Magic methods ----
    def __len__(self):
        return len(self.points)

    def __getitem__(self, index):
        return self.points[index]

    def __iter__(self):
        return iter(self.points)

    def __repr__(self):
        return f"Curve({len(self.points)} points)"

    # ---- Convenience accessors ----
    @property
    def x(self) -> np.ndarray:
        return np.array([p.x for p in self.points])

    @property
    def y(self) -> np.ndarray:
        return np.array([p.y for p in self.points])

    @property
    def z(self) -> np.ndarray:
        return np.array([p.z for p in self.points])


        
        
