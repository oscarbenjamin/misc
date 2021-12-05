#!/usr/bin/env python

from sympy import *
from random import gauss

# Paricle moving at constant velocity in 1D:

x0 = 0  # initial position
v0 = 1  # initial velocity
X0 = Matrix([[x0], [v0]])     # initial state
F = Matrix([[1, 1], [0, 1]])  # state transition matrix
H = Matrix([[1, 0]])  # Observation matrix (observe position)

# Noisy perturbations to velocity:

sigma_state = 0.01  # standard deviation of change in velocity
Q = Matrix([[0, 0], [0, sigma_state**2]])  # dynamic noise covariance
wk = lambda: Matrix([[0], [gauss(0, sigma_state)]])  # noise function

# Noisy errors in sensor:

sigma_sensor = 0.1
R = Matrix([[sigma_sensor]])
vk = lambda: Matrix([[gauss(0, sigma_sensor)]])

Xk = X0 # true initial conditions
Xhat = Matrix([[0], [0]])    # Prior estimate of state
P = Matrix([[1, 0], [0, 1]]) # Prior covariance for state

Zk = None  # No measurements yet

n = 0
while True:
    print(f'\n --------time {n} ---------------------------- ')
    print('Z0:')
    pprint(Zk)
    print('X0:')
    pprint(Xk.T)
    print('Xhat:')
    pprint(Xhat.T)
    print('P:')
    pprint(P)
    input('\nPress enter to continue...') # pause here

    # Update simulation:
    Xk = F*Xk + wk() # Update state
    Zk = H*Xk + vk() # Make sensor reading

    # Kalman predict:
    Xhat = F*Xhat    # Predict new state
    P = F*P*F.T + Q  # Estimate prior covariance

    # Kalman update:
    Yhat = Zk - H*Xhat  # prefit residual
    Sk = H*P*H.T + R    # prefit covariance
    Kk = P*H.T*Sk**-1   # optimal gain
    Xhat = Xhat + Kk*Yhat  # posterior estimate
    P = (eye(2) - Kk*H)*P  # posterior covariance
    n += 1
