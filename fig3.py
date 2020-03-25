import numpy as np
import matplotlib.pyplot as plt
from rocket import max_radius,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

title = ['Mars','Jupiter','Venus','Mercury']
tf = [192.134, 477.62, 139.315, 252.031]
min = [False,False,True,True]
step = [4,4,4,2]

plt.figure(figsize=(6.4,6.4))

for i,t in enumerate(tf):
    t *= DAY
    s = max_radius(accel,t,min=min[i])
    t,r,th = s.x,s.y[0],s.th
    tau = s.tau/DAY
    a = accel(t)
    phi = np.pi/2 - s.phi + th

    x,y = r*np.cos(th),r*np.sin(th)
    a,b = a*np.cos(phi),a*np.sin(phi)

    th = th[-1] + (t - t[-1])/np.sqrt(r[-1]**3)
    th = th[th-th[0]<=2*np.pi]        
    t = t[t<=2*np.pi]

    plt.subplot(2,2,i+1)
    plt.axis('equal')
    plt.plot(r[-1]*np.cos(th),r[-1]*np.sin(th),'--g')
    plt.plot(np.cos(t),np.sin(t),'--k')
    plt.plot(x,y,'b')

    x,y = x[::step[i]],y[::step[i]]
    a,b = a[::step[i]],b[::step[i]]
    plt.quiver(x,y,a,b,width=0.006,color='r')

    ha = 'right' if i&1 else 'left'
    tx = plt.gca().transAxes
    plt.text(i&1,0.07,r'$T$=%.2f days'%tf[i],ha=ha,transform=tx)
    plt.text(i&1,0.01,r'$\tau$=%.2f days'%tau,ha=ha,transform=tx)
    plt.title(title[i], fontsize=14)

plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
