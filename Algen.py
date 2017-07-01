import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N0 = 10
j = 0.5
r = 0.07
P0 = 0.1
K = 1.0
c = 0.5
f = 0.5
g = 0.07
Q = 2
alpha1 = 0.04
k = 0.10
l = 0.12
phi = 47.0
m = 0.03
b = 0.9

def zeta(t):
    return 0.3*np.sin(t/365 * 2 * np.pi)

def alpha(t, M, P):
    return 50.0 * P / M * (1.5 + np.sin(t/365 * 2 * np.pi - np.pi * 0.5))

def dM(M, N, P, H, t):
    return zeta(t)

def dN(M, N, P, H, t):
    return -((alpha(t, M, P) * N) / (j + N) - r) * P + ((m + max(zeta(t), 0)) / M) * (N0 - N)

def dP(M, N, P, H, t):
    return ((alpha(t, M, P) * N) / (j + N) - r) * P - (c * (P - P0) * H) / (K + P - P0) - ((m + max(zeta(t), 0)) / M) * P

def dH(M, N, P, H, t):
    return (f * c * (P - P0) * H) / (K + P - P0) - g * H - zeta(t) / M * H

def solve(x, t):
    M, N, P, H = x
    return (dM(M, N, P, H, t),
            dN(M, N, P, H, t),
            dP(M, N, P, H, t),
            dH(M, N, P, H, t))

time = np.linspace(0,10*365,10000)
solution = odeint(solve, (30,0.5,0.12,0.18), time)
plt.show(plt.plot(time, solution[:,0]))
plt.show(plt.plot(time, solution[:,1], time, solution[:,2], time, solution[:,3]))
