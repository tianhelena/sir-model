import numpy as np 
import math
from math import exp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
N=1000
I0 ,R0=1, 0 
S0=N-I0-R0
beta, gamma=0.5, 1/10 # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
t = np.linspace(0, 160, 160)
apha=0.9 ##effectiveness of vaccine

# The SIR model differential equations 
def deriv(y, t, N, beta, gamma, apha): 
    S, I, R =y 
    dSdt=-beta*(1-apha)*I*S/N
    dIdt=beta*(1-apha)*I*S/N-gamma*I
    dRdt=gamma*I
    return dSdt, dIdt, dRdt 
##initial condition vector 
y0=S0, I0, R0

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma,apha))
S, I, R = ret.T


reta = odeint(deriv, y0, t, args=(N, beta, gamma,0.6))
Sa, Ia, Ra = reta.T

retb = odeint(deriv, y0, t, args=(N, beta, gamma,0.9))
Sb, Ib, Rb = retb.T


##y= odeint(model, y0, t)
###model: Function name that returns derivative value at requested y and t values as dydt=model(y,t)

 





# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, Ia/1000, 'b',  lw=2, label='Infected_a')
ax.plot(t, Ib/1000, 'r', lw=2, label='Infected_b')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)

