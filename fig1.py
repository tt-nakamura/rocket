import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from rocket import max_radius,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

t = np.linspace(20, 490, 50)
s1,s2 = None,None
r1,r2 = [],[]

for tf in t*DAY:
    s1 = max_radius(accel,tf,init=s1)
    s2 = max_radius(accel,tf,init=s2,min=True)
    r1.append(s1.y[0,-1])
    r2.append(s2.y[0,-1])

plt.axis([0,t[-1],0,r1[-1]])
plt.plot(t, r1, 'r', label='maximum')
plt.plot(t, r2, 'b', label='minimum')

t1 = interp1d(r1,t,3)
t2 = interp1d(r2,t,3)

r = [1.52, 5.2, 0.723, 0.387]# mars,jupiter,venus,mercury
t = [t1(r[0]),t1(r[1]),t2(r[2]),t2(r[3])]
print(t)

plt.plot(t, r, '.k')
plt.xlabel(r'$T$ = flight time  / day', fontsize=14)
plt.ylabel(r'$r_1$ = final orbital radius  / AU', fontsize=14)
plt.text(t[0], r[0], 'Mars', va='top', fontsize=14)
plt.text(t[1], r[1], 'Jupiter', ha='right', fontsize=14)
plt.text(t[2], r[2], 'Venus', fontsize=14)
plt.text(t[3], r[3], 'Mercury', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
