import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

t = np.linspace(0,300,3001)
Tm,b = 11.16, 0.62

m,c,w,s = 0.02, 0.41, 3.58, 0.86

x = lambda t: m*np.sin(t/w-s)+0.015*np.cos(t/w/1.3)+c

def dy(t, y):
    x = 0.02 * np.sin(t/3.58) + 0.41
    if x >= y:
        b = 0.2
        b = 1
    else:
        b = 0.03
        b = 0.1
    return b * (x - y)

def dy(t, y):
    x = 0.02 * np.sin(t) + 0.4
    if x >= y:
        b = 0.2
    else:
        b = 0.03
    return b * (x - y)

def rk4(t,y,args):
    k1 =  dy(t, y, *args)
    k2 =  dy(t, y + k1/2, *args)
    k3 =  dy(t, y + k2/2, *args)
    k4 =  dy(t, y + k3, *args)
    return (k1 + 2*k2 + 2*k3 + k4)/6

def solve_rk4(dy,t_span,y0,args,max_step):
    dt = t_span[1]-t_span[0]
    sol = np.empty_like(t_span)
    sol[0] = y0
    sub_step_n = int(dt/max_step)
    for i,t in enumerate(t_span[:-1]):
        yi = sol[i]
        for j in range(sub_step_n):
            yi += max_step*rk4(t,yi,args)
        sol[i+1] = yi
    return sol

sol1 = solve_ivp(dy, t[[0,-1]], (0.42,),
                 t_eval=t).y.squeeze()
sol2 = solve_ivp(dy, t[[0,-1]], (0.4200001,),
                 t_eval=t).y.squeeze()
sol3 = solve_rk4(dy, t, 0.42, (), 0.1)
sol4 = solve_rk4(dy, t, 0.43, (), 0.1)
plt.plot(t,sol1,label='solve_ivp 0.42')
plt.plot(t,sol2,label='solve_ivp 0.4200001')
plt.plot(t,sol3,label='RK4 0.42')
plt.plot(t,sol4,label='RK4 0.43')
plt.legend()
plt.show()
