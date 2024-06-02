import numpy as np
from scipy.interpolate import interp1d
'''
+---------------------------+
| Created by Jade Nassif    |
|                           |
| jade.nassif2002@gmail.com |
+---------------------------+
'''

'''
The purpose of this code is to generate waverider geometries based on the osculating
cone inverse design method. The user inputs the following:
- Freestream Mach number 'M_inf'
- Shock angle 'beta'
- Height of the waverider at the base plane 'height"
- Width of the waverider at the base plane 'width"
- Design parameters 'X1', 'X2', 'X3', 'X4'
- Number of osculating planes 'n_planes'
- Number of points in the streamwise direction 'n_streamwise'

The parametrisation is based on the work by Son et. al [1]

[1] Jiwon Son, Chankyu Son, and Kwanjung Yee. 
'A Novel Direct Optimization Framework for Hypersonic Waverider Inverse Design Methods'.
In: Aerospace 9.7 (June 2022), p. 348. issn: 2226-4310. doi: 10.3390/aerospace9070348.

The output is a CAD geometry of the waverider defined by surfaces.

The code structure is based on the class "waverider" 

Note that the following convention is used in this code:
x --> streamwise direction
y --> transverse direction
z --> spanwise direction
with origin at the waverider tip

A local 2D coordinate system with origin at the shockwave symmetry plane also exists with 
y_bar=y-height, z_bar=z and x=length 
'''

class waverider():
    
    # constructor
    # expected input for dp (design parameters) is a list [X1,X2,X3,X4]
    # STATUS OF FUNCTION : STABLE
    def __init__(self,M_inf,beta,height,width,dp,n_upper_surface,**kwargs):

        #initialise class attributes below
        self.M_inf=M_inf
        self.beta=beta
        self.height=height
        self.width=width

        self.X1=dp[0]
        self.X2=dp[1]
        self.X3=dp[2]
        self.X4=dp[3]

        # check that condition for inverse design is respected
        if not ((self.X2/((1-self.X1)**4))<(7/64)*(self.width/self.height)**4):
            raise ValueError("Condition for inverse design not respected, check design parameters X1 and X2")
    
        # check optional input "n_planes"
        if "n_planes" in kwargs:
            n_planes = kwargs["n_planes"]
            if not (isinstance(n_planes, int) and n_planes >= 10):
                raise TypeError("The number of planes must be an integer and at least 10")
            self.n_planes = n_planes

        # check optional input "n_streamwise"
        if "n_streamwise" in kwargs:
            n_streamwise = kwargs["n_streamwise"]
            if not (isinstance(n_streamwise, int) and n_streamwise >= 10):
                raise TypeError("The number of streamwise points must be an integer and at least 10")
            self.n_streamwise = n_streamwise

        # obtain length of waverider from tip to base plane
        self.length=height/np.tan(self.beta*np.pi/180)

        ''''define the shockwave based on the control points 
        -----------------------------------------------------'''
        self.s_cp=np.zeros((5,2))

        # all five points are evenly distributed in z
        self.s_cp[:,0]=np.transpose(np.linspace(self.X1*self.width,self.width,5))

        # assign the y_bar of the last point
        self.s_cp[-1,1]=self.X2*self.height

        # express the control points as individual points
        # column 1 is z and column 2 is y_bar
        self.s_P0=self.s_cp[0,:]
        self.s_P1=self.s_cp[1,:]
        self.s_P2=self.s_cp[2,:]
        self.s_P3=self.s_cp[3,:]
        self.s_P4=self.s_cp[4,:]

        ''' define the upper surface curve
        ---------------------------------------------------'''
        self.us_cp=np.zeros((4,2))

        # assign z coordinates of all points defining upper surface, equally spaced
        self.us_cp[:,0]=np.transpose(np.linspace(0,self.width,4))

        #assign y_bar coordinates of all points on upper surface
        self.us_cp[0,1]=self.height
        self.us_cp[1,1]=self.height-(1-self.X2)*self.X3
        self.us_cp[2,1]=self.height-(1-self.X2)*self.X4

        #assign last point using the P4 computed for the shockwave
        self.us_cp[3,:]=self.s_P4

        #define control points individually
        self.us_P0=self.us_cp[0,:]
        self.us_P1=self.us_cp[1,:]
        self.us_P2=self.us_cp[2,:]
        self.us_P3=self.us_cp[3,:]

        # create an interpolation object for upper surface
        self.Create_Interpolated_Upper_Surface(n=n_upper_surface)



    def Create_Interpolated_Upper_Surface(self,n):

        # values of t for the bezier curve
        t_values=np.linspace(0,1,n)

        points=np.zeros((n,2))

        # get points along the bezier curve representing the upper surface
        for i, t in enumerate(t_values):
            points[i, :] = self.Bezier_Upper_Surface(t)

        # store interp1d objected as an attribute
        self.Interpolate_Upper_Surface=interp1d(points[:,0],points[:,1],kind='linear')







    """
    AUXILIARY FUNCTIONS    
    """    
    def Bezier_Shockwave(self,t):

        point=(1-t)**4*self.s_P0+4*(1-t)**3*t*self.s_P1+6*(1-t)**2*t**2*self.s_P2+4*(1-t)*t**3*self.s_P3+t**4*self.s_P4

        return point
    
    # returns slope, 
    def First_Derivative(self, t):

        first_derivative = 4 * (1-t)**3 * (self.s_P1 - self.s_P0) + 12 * (1-t)**2 * t * (self.s_P2 - self.s_P1) + 12 * (1-t) * t**2 * (self.s_P3 - self.s_P2) + 4 * t**3 * (self.s_P4 - self.s_P3)
    
        return first_derivative

    def Bezier_Upper_Surface(self, t):

        point = (1 - t)**3 * self.us_P0 + 3 * (1 - t)**2 * t * self.us_P1 + 3 * (1 - t) * t**2 * self.us_P2 + t**3 * self.us_P3

        return point
    
    def Local_to_Global(self,y):
    # convert local coordinates to global coordinates
        y=y-self.height

        return y
        

# Auxiliary Functions
def Euclidean_Distance(x1,y1,x2,y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

def Equation_of_Line(x,m,c):
    return m*x+c