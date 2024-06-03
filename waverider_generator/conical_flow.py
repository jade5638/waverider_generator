'''
This file was taken from 
https://github.com/j-herrera/waverider
to apply conical flow locally for each osculating cone
'''
import numpy as np
from scipy.optimize import fsolve
from  scipy.integrate import solve_ivp
from scipy.interpolate import UnivariateSpline

def TM(t, x, gamma):
    A = (gamma - 1.0) / 2.0 * (1.0 - pow(x[0], 2) - pow(x[1], 2))

    dxdt = np.zeros(2)
    dxdt[0] = x[1]
    dxdt[1] = (x[1] * x[0] * x[1] - A * (2.0 * x[0] +
                                         x[1] / np.tan(t))) / (A - x[1] * x[1])
    return dxdt


def cone_field(Mach, theta, beta, gamma):
    d = np.arctan(2.0 / np.tan(beta) * (pow(Mach, 2) * pow(np.sin(beta),
                                                           2) - 1.0) / (pow(Mach, 2) * (gamma + np.cos(2 * beta)) + 2.0))
    Ma2 = 1.0 / np.sin(beta - d) * np.sqrt((1.0 + (gamma - 1.0) / 2.0 * pow(Mach, 2) * pow(
        np.sin(beta), 2)) / (gamma * pow(Mach, 2) * pow(np.sin(beta), 2) - (gamma - 1.0) / 2.0))
    V = 1.0 / np.sqrt(2.0 / ((gamma - 1.0) * pow(Ma2, 2)) + 1.0)
    Vr = V * np.cos(beta - d)
    Vt = -(V * np.sin(beta - d))

    xt = np.array([Vr, Vt])
    
    sol = solve_ivp(TM, (beta, theta), xt, args=(gamma,))
    Vrf = UnivariateSpline(sol.t[::-1], sol.y[0, ::-1], k=min(3, sol.t.size-1))
    Vtf = UnivariateSpline(sol.t[::-1], sol.y[1, ::-1], k=min(3, sol.t.size-1))
    return [Vrf, Vtf]


def Vt0(t, y, gamma):
    return y[1]
    
Vt0.terminal = True

def f(x, Mach, theta, gamma):
    d = np.arctan(2.0 / np.tan(x) * (pow(Mach, 2) * pow(np.sin(x), 2) -
                                     1.0) / (pow(Mach, 2) * (gamma + np.cos(2 * x)) + 2.0))
    Ma2 = 1.0 / np.sin(x - d) * np.sqrt((1.0 + (gamma - 1.0) / 2.0 * pow(Mach, 2) * pow(
        np.sin(x), 2)) / (gamma * pow(Mach, 2) * pow(np.sin(x), 2) - (gamma - 1.0) / 2.0))
    V = 1.0 / np.sqrt(2.0 / ((gamma - 1.0) * pow(Ma2, 2)) + 1.0)
    Vr = V * np.cos(x - d)
    Vt = -(V * np.sin(x - d))

    xt = np.array([Vr[0], Vt[0]])
    
    sol = solve_ivp(TM, (x, 0.0), xt, events=Vt0, args=(gamma,))
    
    return sol.t[-1] - theta


def shock_angle(Mach, theta, gamma):
    beta = fsolve(f, theta, (Mach, theta, gamma))
    return beta
