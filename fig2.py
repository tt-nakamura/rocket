import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from rocket import max_radius,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

tf = 192.13*DAY

s = max_radius(accel,tf)
t = s.x/DAY
r,u,v,lr,lu,lv = s.y
th,phi = s.th,s.phi

plt.figure(figsize=(6.4, 8))
plt.subplots_adjust(left=0.12, right=0.98,
                    bottom=0.08, top=0.98, hspace=0)

plt.subplot(3,1,1)
plt.plot(t,r,label='r')
plt.plot(t,u,label='u')
plt.plot(t,v,label='v')
plt.xlim(t[[0,-1]])
plt.xticks([])
plt.ylabel('r, u, v', fontsize=14)
plt.legend(labelspacing=0.2)

plt.subplot(3,1,2)
plt.plot(t,lr,label=r'$\lambda_r$')
plt.plot(t,lu,label=r'$\lambda_u$')
plt.plot(t,lv,label=r'$\lambda_v$')
plt.xlim(t[[0,-1]])
plt.xticks([])
plt.ylabel(r'$\lambda_r$, $\lambda_u$, $\lambda_v$', fontsize=14)
plt.legend(labelspacing=0.2)

plt.subplot(3,1,3)
plt.plot(t,phi)
plt.xlim(t[[0,-1]])
plt.xticks()
plt.xlabel(r'$t - t_0$  / day', fontsize=14)
plt.ylabel(r'$\phi$  / radian', fontsize=14)

plt.savefig('fig2.eps')
plt.show()
