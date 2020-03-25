import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from rocket import max_radius,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

tf = np.linspace(20, 490, 50)
s1,s2 = None,None
r1,r2,t1,t2 = [],[],[],[]

for t in tf*DAY:
    s1 = max_radius(accel,t,init=s1)
    s2 = max_radius(accel,t,init=s2,min=True)
    r1.append(s1.y[0,-1])
    r2.append(s2.y[0,-1])
    t1.append(s1.tau/DAY)
    t2.append(s2.tau/DAY)

plt.axis([0,r1[-1],0,t2[-1]])
plt.plot(r1, t1, 'r', label='maximum')
plt.plot(r2, t2, 'b', label='minimum')

t1 = interp1d(r1,t1,3)
t2 = interp1d(r2,t2,3)

r = [1.52, 5.2, 0.723, 0.387]# mars,jupiter,venus,mercury
t = [t1(r[0]),t1(r[1]),t2(r[2]),t2(r[3])]
print(t)

plt.plot(r, t, '.k')
plt.xlabel(r'$r_1$ = final orbital radius  / AU', fontsize=14)
plt.ylabel(r'$\tau$ = launch time before conjunction / day', fontsize=13)
plt.text(r[0], t[0], 'Mars', va='top', fontsize=14)
plt.text(r[1], t[1], 'Jupiter', ha='right', fontsize=14)
plt.text(r[2], t[2], 'Venus', fontsize=14)
plt.text(r[3], t[3], 'Mercury', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
