import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from  scipy.integrate import solve_ivp
from waverider_generator.flowfield import cone_angle,cone_field
from typing import Union

'''
+-----------------------------------+
| Created by Jade Nassif            |
| Github : jade5638                 |
| Email : jade.nassif2002@gmail.com |
+-----------------------------------+
'''
'''
Documentation on the inputs is provided in the README file

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
    def __init__(self,M_inf: Union[float,int],beta: Union[float,int],height: Union[float,int],width: Union[float,int],dp:list,n_upper_surface:int,n_shockwave:int,**kwargs):

        ''''
        +-------------------+
        | initialise inputs |
        +-------------------+
        '''

        if not isinstance(M_inf, (float, int)) or M_inf <= 0:
            raise ValueError("Mach number must be a positive number")
        else:
            self.M_inf=float(M_inf)

        if not isinstance(beta, (float, int)) or not (0 < beta < 90):
            raise ValueError("beta must be a float or integer between 0 and 90 degrees.")
        else:
            self.beta=float(beta)

        if not isinstance(height,(float,int)) or height<=0:
            raise ValueError("height must be a positive number")
        else:
            self.height=float(height)
        
        if not isinstance(width,(float,int)) or width<=0:
            raise ValueError("width must be a positive number")
        else:
            self.width=float(width)
        
        error=""
        if not isinstance(dp,list):
            error=error+'dp must be a list\n'
        
        if len(dp)!=4:
            error=error+'dp must contain 4 elements'

        if error!="":
            raise ValueError(error)
        
        for parameter in dp:
            if not isinstance(parameter,(float,int)):
                raise ValueError('please enter a valid float or int for the design parameters')
            
        # extract the design parameters                
        self.X1=dp[0]
        self.X2=dp[1]
        self.X3=dp[2]
        self.X4=dp[3]

        # check that condition for inverse design is respected
        if (self.X1<1 and self.X1>=0) and (self.X2<=1 and self.X2>=0):
            if not ((self.X2/((1-self.X1)**4))<(7/64)*(self.width/self.height)**4):
                raise ValueError("Condition for inverse design not respected, check value of design parameters X1 and X2")
        else:
            raise ValueError("X1 and/or X2 are not in the required range")

        if not (self.X3>=0 and self.X3<=1):
            raise ValueError("X3 must be between 0 and 1")
        
        if not (self.X4>=0 and self.X4<=1):
            raise ValueError("X4 must be between 0 and 1")
        
        if not isinstance(n_upper_surface,int) or n_upper_surface<10:
            raise ValueError('number of points on the upper surface for interpolation must be an integer greater than or equal to 10')
        
        if not isinstance(n_shockwave,int) or n_shockwave<10:
            raise ValueError('number of points on the shockwave for interpolation must be an integer greater than or equal to 10')
                
        # check optional input "n_planes"
        if "n_planes" in kwargs:
            n_planes = kwargs["n_planes"]
            if not (isinstance(n_planes, int) and n_planes >= 10):
                raise TypeError("The number of planes must be an integer and at least 10")
            self.n_planes = n_planes
        else:
            self.n_planes = 10

        # check optional input "n_streamwise"
        if "n_streamwise" in kwargs:
            n_streamwise = kwargs["n_streamwise"]
            if not (isinstance(n_streamwise, int) and n_streamwise >= 10):
                raise TypeError("The number of streamwise upper surface points must be an integer and at least 10")
            self.n_streamwise = n_streamwise
        else:
            self.n_streamwise = 10

        # check optional input "delta_streamwise"
        if "delta_streamwise" in kwargs:
            delta_streamwise = kwargs["delta_streamwise"]
            if isinstance(delta_streamwise, float) and (delta_streamwise<=0.2 and delta_streamwise>0):
                self.delta_streamwise = delta_streamwise
            else:
                raise ValueError("delta_streamwise must be a percentage between 0 and 20 percent of the waverider length")
        else:
            self.delta_streamwise = 0.05


        #ratio of specific heats
        self.gamma=1.4

        #computes self.theta, the deflection angle corresponding to a shock angle in oblique shock relations
        self.Compute_Deflection_Angle()

        # obtain length of waverider from tip to base plane
        self.length=height/np.tan(self.beta*np.pi/180)

        ''''
        +--------------------------------------------------+
        | define the shockwave based on the control points |
        +--------------------------------------------------+
        '''
        # stores coordinates of all the control points in format z_bar,y_bar
        # four of the points are y_bar=0 already
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

        ''''
        +---------------------------------------------+
        | define the upper surface via control points |
        +---------------------------------------------+
        '''
        # same procedure as with shockwave curve
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

        ''''
        +-------------------------------------------------------------------+
        | create interpolation objects for shockwave and upper surface curve|
        +-------------------------------------------------------------------+
        '''
        # create an interpolation object for upper surface curve and curved part of the shockwave:
        # self.Interpolate_Upper_Surface
        # self.Interpolate_Shockwave
        self.Create_Interpolated_Upper_Surface(n=n_upper_surface)
        self.Create_Interpolated_Shockwave(n=n_shockwave)

        ''''
        +------------------------------------------------------------------+
        | find intersections of osculating planes with upper surface curve |
        +------------------------------------------------------------------+
        '''
        # next step is to calculate intersections with upper surface
        # start by defining the sample of points for the shockwave in local coordinates z and y
        self.z_local_shockwave=np.linspace(0,self.width,self.n_planes+2)
        self.z_local_shockwave=self.z_local_shockwave[1:-1] 

        #obtain the y_bar values for the z sample in self.y_local_shockwave
        self.y_local_shockwave=np.zeros((self.n_planes,1))
        self.Get_Shockwave_Curve()

        #obtain the intersection with the upper surface in self.local_intersections_us
        self.local_intersections_us=np.zeros((self.n_planes,2))
        self.Find_Intersections_With_Upper_Surface()
        ''''
        +-------------------------------------------+
        | compute the leading edge and cone centers |
        +-------------------------------------------+
        '''
        # next step is to obtain the LE points in the global coordinate system
        # initialise LE object
        self.leading_edge=np.zeros((self.n_planes+2,3))

        # tip is already at 0,0,0
        # set the last point (tip at base plane)
        self.leading_edge[-1,:]=np.array([self.length,self.Local_to_Global(self.X2*self.height),self.width])

        # initialise an object for the cone centers, note in the flat these are not "cones" but still better to have a single array
        self.cone_centers=np.zeros((self.n_planes,3))

        # osculate through the planes and populate self.cone_centers and self.leading_edge
        self.Compute_Leading_Edge_And_Cone_Centers()

        ''''
        +---------------------------+
        | compute the upper surface |
        +---------------------------+
        '''
        # next step is to compute the upper surface
        # stored in this format for easy visualisation in matplotlib
        self.upper_surface_x=np.zeros((self.n_planes+1,self.n_streamwise))
        self.upper_surface_y=np.zeros((self.n_planes+1,self.n_streamwise))
        self.upper_surface_z=np.zeros((self.n_planes+1,self.n_streamwise))

        # add the symmetry plane
        self.upper_surface_x[0,:]=np.linspace(0,self.length,self.n_streamwise)
        self.upper_surface_y[0,:]=0
        self.upper_surface_z[0,:]=0

        self.Compute_Upper_Surface()

        # add the tip point
        x_tip = np.full((1, self.n_streamwise), self.length)
        y_tip =np.full((1,self.n_streamwise),self.height*self.X2-self.height)
        z_tip =np.full((1,self.n_streamwise),self.width)

        self.upper_surface_x= np.vstack([self.upper_surface_x, x_tip])
        self.upper_surface_y= np.vstack([self.upper_surface_y, y_tip])
        self.upper_surface_z= np.vstack([self.upper_surface_z, z_tip])

        # store in a streams format 
        self.upper_surface_streams=[]
        self.Streams_Format()
        
        '''
        +---------------------------+
        | compute the lower surface |
        +---------------------------+
        '''
        
        # next step is to trace the streamlines

        # compute the cone angle in degrees
        self.cone_angle=cone_angle(self.M_inf,self.beta,self.gamma)

        self.lower_surface_streams=[]
        
        # populate the self.lower_surface_streams list by tracing the streamlines
        self.Streamline_Tracing()

    # convert the upper surface streams to the desired format
    def Streams_Format(self):

        for i in range(self.n_planes+2):

            x=self.upper_surface_x[i,:]
            y=self.upper_surface_y[i,:]
            z=self.upper_surface_z[i,:]
            self.upper_surface_streams.append(np.vstack([x,y,z]).T)

            # keep only twice the same point for the tip
            if i==self.n_planes+1:
                self.upper_surface_streams[i]=self.upper_surface_streams[i][0:2,:]

    def Streamline_Tracing(self):

        # propagate the streamlines
        Vr, Vt = cone_field(self.M_inf,self.cone_angle*np.pi/180,self.beta*np.pi/180,self.gamma)

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
        
        # make the following arrays of size n_planes+2
        leading_edge=self.leading_edge
        y_local_shockwave=self.y_local_shockwave
        y_local_shockwave=np.vstack((np.array([[0]]),y_local_shockwave,np.array([[self.X2*self.height]])))

        z_local_shockwave=self.z_local_shockwave[:,None]
        z_local_shockwave=np.vstack((np.array([[0]]),z_local_shockwave,np.array([[self.width]])))

        cone_centers=self.cone_centers
        cone_centers=np.vstack((np.array([[0,0,0]]),cone_centers,np.array([[self.length,self.Local_to_Global(self.X2*self.height),self.width]])))

        local_intersections_us=self.local_intersections_us
        local_intersections_us=np.vstack((np.array([[0,self.height]]),local_intersections_us,np.array([[self.width,self.X2*self.height]])))

        for i,le_point in enumerate(leading_edge):
            # tip
            if i==len(leading_edge)-1:
                stream=np.vstack((le_point,le_point))
                self.lower_surface_streams.append(stream)

            # flat region
            elif z_local_shockwave[i,0]<=self.X1*self.width or self.X2==0:

                # trigonometry with deflection angle
                bottom_surface_y=le_point[1]-np.tan(self.theta*np.pi/180)*(self.length-le_point[0])

                # store the x,y and z in a streams
                x=np.linspace(le_point[0],self.length,self.n_streamwise)[:,None]
                y=np.linspace(le_point[1],bottom_surface_y,self.n_streamwise)[:,None]
                z=np.full((y.shape),le_point[2])

                self.lower_surface_streams.append(np.column_stack([x,y,z]))

            # curved region
            else:

                # need calculate R minus height of osculating plane
                eta_le=Euclidean_Distance(
                    local_intersections_us[i,0],
                    self.Local_to_Global(local_intersections_us[i,1]),
                    cone_centers[i,2],
                    cone_centers[i,1]
                ) 
                r=Euclidean_Distance(
                    z_local_shockwave[i,0],
                    self.Local_to_Global(y_local_shockwave[i,0]),
                    cone_centers[i,2],
                    cone_centers[i,1]
                ) 

                # calculate the angle to rotate the streamlines by
                m,_,_=self.Get_First_Derivative(z_local_shockwave[i,0])
                alpha=np.arctan(m)

                x_le=(eta_le)/ np.tan(self.beta*np.pi/180) 

                sol = solve_ivp(stode, (0, 1000), [x_le, eta_le], events=back, args=(r / np.tan(self.beta*np.pi/180),), max_step=self.delta_streamwise*self.length)
                stream = np.vstack([sol.y[0], -sol.y[1] * np.cos(alpha), sol.y[1] * np.sin(alpha)]).T

                # transform from cone center coordinate system to global
                stream[:,0]=stream[:,0]+cone_centers[i,0]
                stream[:,1]=stream[:,1]+cone_centers[i,1]
                stream[:,2]=stream[:,2]+cone_centers[i,2]

                # append
                self.lower_surface_streams.append(stream)

    def Compute_Upper_Surface(self):
        
        for i in range(0,self.n_planes):
            self.upper_surface_x[i+1,:]=np.linspace(self.leading_edge[i+1,0],self.length,self.n_streamwise)
            self.upper_surface_y[i+1,:]=np.linspace(self.leading_edge[i+1,1],self.Local_to_Global(self.local_intersections_us[i,1]),self.n_streamwise)
            self.upper_surface_z[i+1,:]=np.linspace(self.leading_edge[i+1,2],self.local_intersections_us[i,0],self.n_streamwise)
        
    def Compute_Leading_Edge_And_Cone_Centers(self):

        for i,z in enumerate(self.z_local_shockwave):
             
            if z<=self.X1*self.width or self.X2==0:
                self.cone_centers[i,0]=self.length-((self.local_intersections_us[i,1]-self.y_local_shockwave[i,0])/np.tan(self.beta*np.pi/180))
                self.cone_centers[i,1]=self.Local_to_Global(self.local_intersections_us[i,1])
                self.cone_centers[i,2]=float(z)

                self.leading_edge[i+1,:]=self.cone_centers[i,:]

            else:

                #calculate corresponding t value
                t=self.Find_t_Value(z)
                # first derivative and radius
                first_derivative,_,_=self.First_Derivative(t)
                radius=self.Calculate_Radius_Curvature(t)

                self.leading_edge[i+1,:]=self.cone_centers[i,:]
                # get angle theta
                theta=np.arctan(first_derivative)
                
                # get x value for cone center
                self.cone_centers[i,0]=float(self.length-radius/np.tan(self.beta*np.pi/180))

                # get y value for cone center
                self.cone_centers[i,1]=float(self.Local_to_Global(self.y_local_shockwave[i,0])+np.cos(theta)*radius) 

                # get z value for cone center
                self.cone_centers[i,2]=float(z-radius*np.sin(theta))

                # get the location of the intersection
                self.leading_edge[i+1,:]=self.Intersection_With_Freestream_Plane(self.cone_centers[i,0],
                                                                                self.cone_centers[i,1],
                                                                                self.cone_centers[i,2],
                                                                                self.length,
                                                                                self.Local_to_Global(self.y_local_shockwave[i,0]),
                                                                                z,
                                                                                self.Local_to_Global(self.local_intersections_us[i,1]))
    
    # find all intersections with the upper surface
    def Find_Intersections_With_Upper_Surface(self):

        for i,z in enumerate(self.z_local_shockwave):

            if z<=self.X1*self.width or self.X2==0:
                self.local_intersections_us[i,0]=z
                self.local_intersections_us[i,1]=self.Interpolate_Upper_Surface(z)
            else:
                first_derivative,_,_=self.Get_First_Derivative(z)
                # print(first_derivative)
                intersection=self.Intersection_With_Upper_Surface(first_derivative=first_derivative,z_s=float(z),y_s=float(self.y_local_shockwave[i,:]))
                self.local_intersections_us[i,:]=intersection

    # create an interp1d object for the shockwave curve 
    def Create_Interpolated_Shockwave(self,n):
        
        t_values=np.linspace(0,1,n)
        points=np.zeros((n,2))

        for i,t in enumerate(t_values):
            points[i,:]=self.Bezier_Shockwave(t)

        self.Interpolate_Shockwave=interp1d(points[:,0],points[:,1],kind='linear')

    # creates an interp1d object for the upper surface, used to find intersection easily
    # with root_scalar
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
    # function which determines intersection of an osculating plane with the upper surface
    # curve
    # inputs:
    # - first derivative at the shockwave point (dy/dz)
    # - z_s and y_s which are the local coordinates of the point
    # outputs the coordinates of the intersection in local coordinates
    def Intersection_With_Upper_Surface(self,first_derivative,z_s,y_s):

        # get the constant c and slope for the line between the two points
        c=y_s+(1/first_derivative)*z_s
        m= -1/first_derivative

        # define the function used to find the intersection between the two curves
        def f(z):
            return Equation_of_Line(z,m,c) - self.Interpolate_Upper_Surface(z)
        
        # use root_scalar method to get the root
        intersection = root_scalar(f, bracket=[0, self.width])

        # extract local coordinates of the root
        z=intersection.root
        y=Equation_of_Line(z,m,c)

        return np.array([z,y])

    # get the first derivative for a z value along the shockwave
    def Get_First_Derivative(self,z):

        # get the corresponding t value
        t=self.Find_t_Value(z)

        first_derivative,dzdt,dydt=self.First_Derivative(t)

        return first_derivative,dzdt,dydt

    # get the y_bar coordinates of all points along the shockwave curve by means of 
    # interpolation
    def Get_Shockwave_Curve(self):

        for i,z in enumerate(self.z_local_shockwave):
            if z<=self.width*self.X1:
                self.y_local_shockwave[i,0]=0
            else:
                self.y_local_shockwave[i,0]=float(self.Interpolate_Shockwave(float(z)))
    
    # bezier curve defining the shockwave
    # returns np.array([z,y])
    def Bezier_Shockwave(self,t):

        point=(1-t)**4*self.s_P0+4*(1-t)**3*t*self.s_P1+6*(1-t)**2*t**2*self.s_P2+4*(1-t)*t**3*self.s_P3+t**4*self.s_P4

        return point
    
    # returns slope m, dz/dt and dy/dt
    def First_Derivative(self, t):

        first_derivative = 4 * (1-t)**3 * (self.s_P1 - self.s_P0) + 12 * (1-t)**2 * t * (self.s_P2 - self.s_P1) + 12 * (1-t) * t**2 * (self.s_P3 - self.s_P2) + 4 * t**3 * (self.s_P4 - self.s_P3)
    
        return first_derivative[1]/first_derivative[0],first_derivative[0],first_derivative[1]
    
    # returns components of second derivative of point along shockwave curve with respect to t
    # z and y respectively
    def Second_Derivative(self, t):

        second_derivative = 12 * (1-t)**2 * (self.s_P2 - 2 * self.s_P1 + self.s_P0) + 24 * (1-t) * t * (self.s_P3 - 2 * self.s_P2 + self.s_P1) + 12 * t**2 * (self.s_P4 - 2 * self.s_P3 + self.s_P2)
        return second_derivative[0],second_derivative[1]
    
    # Bezier curve of upper surface
    # output is an np.array([z,y]) in local coordinates
    def Bezier_Upper_Surface(self, t):

        point = (1 - t)**3 * self.us_P0 + 3 * (1 - t)**2 * t * self.us_P1 + 3 * (1 - t) * t**2 * self.us_P2 + t**3 * self.us_P3

        return point
    
    #convert from local to global y coordinate
    def Local_to_Global(self,y):

        y=y-self.height

        return y
    
    # get the deflection angle resulting from the shock angle and flow conditions
    def Compute_Deflection_Angle(self):

        tanTheta=2*cot(self.beta*np.pi/180)*(self.M_inf**2*np.sin(self.beta*np.pi/180)**2-1)/(self.M_inf**2*(self.gamma+np.cos(2*self.beta*np.pi/180))+2)

        self.theta=np.arctan(tanTheta)*180/np.pi
    
    # finds the intersection with the local freestream plane by finding the point where y=y_target=y_upper_surface
    def Intersection_With_Freestream_Plane(self,x_C,y_C,z_C,x_S,y_S,z_S,y_target):

        #  ALL COORDINATES IN GLOBAL SYSTEM
        # x_C,y_C,z_C are coordinates of cone center
        # x_S,y_S,z_S are coordinates of shock location in osculating plane

        # need to find where y=y_target
        # parametric curve
        k=(y_target-y_S)/(y_C-y_S)

        x_I=x_S+k*(x_C-x_S)
        y_I=y_target
        z_I=z_S+k*(z_C-z_S)

        return np.array([x_I,y_I,z_I])


    # calculate the radius of curvature for a given t along the bezier curve
    def Calculate_Radius_Curvature(self,t):
        
        _,dzdt,dydt=self.First_Derivative(float(t))
        dzdt2,dydt2=self.Second_Derivative(float(t))

        radius= 1/(abs((dzdt*dydt2-dydt*dzdt2))/((dzdt**2+dydt**2)**(3/2)))

        return radius

    # find the t value which corresponds to a z value
    def Find_t_Value(self,z):

        def f(t):
            return self.Bezier_Shockwave(t)[0]-z
        
        intersection=root_scalar(f,bracket=[0,1])

        return intersection.root

'''EXTERNAL AUXILIARY FUNCTIONS'''

# calculates the euclidean distance between two points in 2D
def Euclidean_Distance(x1,y1,x2,y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

# Equation of a straight line
def Equation_of_Line(z,m,c):
    return m*z+c

# cotangent
def cot(angle):
    return 1/np.tan(angle)
