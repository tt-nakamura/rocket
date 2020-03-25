import numpy as np
import matplotlib.pyplot as plt
from rocket import max_radius,min_time,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

title = ['Mars','Jupiter','Venus','Mercury']
rf = [1.52, 5.2, 0.723, 0.387]
tf = [192.134, 477.62, 139.315, 252.031]
dt = [np.linspace(-30,50,40),
      np.linspace(-40,70,40),
      np.linspace(-30,50,40),
      np.linspace(-30,45,40)]

plt.figure(figsize=(6.4,6.4))

for i,r in enumerate(rf):
    t = tf[i]*DAY
    s = max_radius(accel,t,min=(r<1))

    T = []
    tau = s.tau/DAY + dt[i]
    for t in tau*DAY:
        s = min_time(accel,r,t,init=s)
        T.append(s.p[0]/DAY)

    plt.subplot(2,2,i+1)
    plt.plot(tau,T)
    plt.xlim(tau[[0,-1]])
    plt.title(title[i], fontsize=14)

    if i&2: plt.xlabel(r'$\tau$ = launch time  / day', fontsize=14)
    if i&1==0: plt.ylabel(r'$T$ = flight time  / day', fontsize=14)

plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
