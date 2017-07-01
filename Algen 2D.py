from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

N0 = 10.0
j = 0.5
r = 0.07
P0 = 0.1
K = 1.0
c = 1.0
f = 0.5
g = 0.07
Q = 2.0
alpha1 = 0.04
k = 0.1
l = 0.12
phi = 47.0
m = 3.0
b = 0.9

def alpha(t):
    #return 0.1*(1 + np.sin(t/365 * 2 * np.pi - 0.8 * np.pi))
    return 0.23*(1.5 + np.sin(t/365 * 2 * np.pi - np.pi * 0.5))

def dP(P, H, t):
    return alpha(t) * P - (c * (P - P0)) / (k + P - P0) * H

def dH(P, H, t):
    return (f * c * (P - P0)) / (k + P - P0) * H - g * H

def solve(x, t):
    P, H = x
    return (dP(P, H, t), dH(P, H, t))

time = np.linspace(0,5*365,10000)
solution = odeint(solve, (0.2,0.25), time)
plt.show(plt.plot(time, solution[:,0], time, solution[:,1]))
