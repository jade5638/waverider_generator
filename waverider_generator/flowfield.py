
#%%

'''
Note: The functions found here are largely based on the conical_flow file
found in https://github.com/j-herrera/waverider by Javier Herrera Montojo

The author of this package acknowledges the contributions of Javier Herrera Montojo 
'''

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.interpolate import UnivariateSpline

# Taylor-Maccoll ODE which describes conical flow

def Taylor_Maccoll(t, x, gamma):

    # calculate A constant
    A = (gamma - 1.0) / 2.0 * (1.0 - x[0]**2 - x[1]**2)

    dxdt = np.zeros(2)

    dxdt[0] = x[1]
    dxdt[1] = (x[1] * x[0] * x[1] - A * (2.0 * x[0] + x[1] / np.tan(t))) / (A - x[1] * x[1])

    return dxdt

def cone_field(Mach, theta_rad, beta_rad, gamma):
    
    # deflection angle
    d = np.arctan(2.0 / np.tan(beta_rad) * (Mach**2 * np.sin(beta_rad)**2 - 1.0) / (Mach**2 * (gamma + np.cos(2 * beta_rad)) + 2.0))

    # post shock mach
    Ma2 = 1.0 / np.sin(beta_rad - d) * np.sqrt((1.0 + (gamma - 1.0) / 2.0 * Mach**2 * np.sin(beta_rad)**2) / (gamma * Mach**2 * np.sin(beta_rad)**2 - (gamma - 1.0) / 2.0))

    # velocity magnitude
    V = 1.0 / np.sqrt(2.0 / ((gamma - 1.0) * Ma2**2) + 1.0)

    # components of velocity
    Vr = V * np.cos(beta_rad - d)
    Vt = -(V * np.sin(beta_rad - d))

    xt = np.array([Vr, Vt])
    
    sol = solve_ivp(Taylor_Maccoll, (beta_rad, theta_rad), xt, args=(gamma,))

    # create spline functions
    Vrf = UnivariateSpline(sol.t[::-1], sol.y[0, ::-1], k=min(3, sol.t.size-1))
    Vtf = UnivariateSpline(sol.t[::-1], sol.y[1, ::-1], k=min(3, sol.t.size-1))

    return [Vrf, Vtf]

# function used in solve_ivp to define the event when tangential velocity becomes 0 
def Vt0(t, y, gamma):
    return y[1]

Vt0.terminal = True

# function used in solve ivp to find the root
def f(x, Mach, theta_rad, gamma):

    # deflection angle
    d = np.arctan(2.0 / np.tan(x) * (Mach**2 * np.sin(x)**2 - 1.0) / (Mach**2 * (gamma + np.cos(2 * x)) + 2.0))

    # calculate post shock mach number
    Ma2 = 1.0 / np.sin(x - d) * np.sqrt((1.0 + (gamma - 1.0) / 2.0 * Mach**2 * np.sin(x)**2) / (gamma * Mach**2 * np.sin(x)**2 - (gamma - 1.0) / 2.0))

    # calculate post shock velocity
    V = 1.0 / np.sqrt(2.0 / ((gamma - 1.0) * Ma2**2) + 1.0)

    # get velocity components fom V (radial and tangential)
    Vr = float(V * np.cos(x - d))
    Vt = float(-(V * np.sin(x - d)))

    # store the initial conditions
    xt = np.array([Vr, Vt])
    # Inputs:
    # (x, 0.0) --> integration interval, starting from shock angle 'x' to 0.0 (cone axis)
    # xt --> initial conditions for the TM equations
    # Vt0 --> event function to terminate the integration when tangential velocity component 'Vt' becomes zero
    sol = solve_ivp(Taylor_Maccoll, (x, 0.0), xt, events=Vt0, args=(gamma,))

    # 'sol.t[-1]' represents the angle where the integration stopped, 
    # which is where the tangential velocity 'Vt' becomes zero.
    return sol.t[-1] - theta_rad

# function which solves for the shock angle given a cone angle theta 
def shock_angle(Mach, theta_deg, gamma):

    # theta in radians
    theta=theta_deg*np.pi/180
    
    beta = fsolve(f, theta, (Mach, theta, gamma))

    # beta in deg
    return beta[0]*180/np.pi

# function which solves for a cone angle given a shock angle
def cone_angle(Mach, shock_angle_deg, gamma):

    shock_angle_rad = shock_angle_deg*np.pi/180
    
    def TM_cone(t, x):
        return Taylor_Maccoll(t, x, gamma)
    
    def cone_event(t, y):
        return Vt0(t, y, gamma)
    
    cone_event.terminal = True

    # deflection angle
    d = np.arctan(2.0 / np.tan(shock_angle_rad) * (Mach**2 * np.sin(shock_angle_rad)**2 - 1.0) / (Mach**2 * (gamma + np.cos(2 * shock_angle_rad)) + 2.0))

    # calculate post-shock Mach number
    Ma2 = 1.0 / np.sin(shock_angle_rad - d) * np.sqrt((1.0 + (gamma - 1.0) / 2.0 * Mach**2 * np.sin(shock_angle_rad)**2) / (gamma * Mach**2 * np.sin(shock_angle_rad)**2 - (gamma - 1.0) / 2.0))

    # calculate post-shock velocity
    V = 1.0 / np.sqrt(2.0 / ((gamma - 1.0) * Ma2**2) + 1.0)

    # calculate velocity components from V (radial and tangential)
    Vr = V * np.cos(shock_angle_rad - d)
    Vt = -(V * np.sin(shock_angle_rad - d))

    # store the initial conditions
    xt = np.array([Vr, Vt])

    sol = solve_ivp(TM_cone, (shock_angle_rad, 0.0), xt, events=cone_event)
    
    # extract cone angle
    cone_angle_rad = sol.t_events[0][0]
    # convert to deg
    cone_angle_deg = np.degrees(cone_angle_rad)
    
    return cone_angle_deg
