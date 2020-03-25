import numpy as np
import matplotlib.pyplot as plt
from rocket import min_time,DAY

F = 0.1405
m_dot = 0.07489
def accel(t): return F/(1 - m_dot*t)

title = ['Mars','Jupiter','Venus','Mercury']
rf = [1.52, 5.2, 0.723, 0.387]
tf = [192.134, 477.62, 139.315, 252.031]
tau = [60,140,120,230,50,130,140,210]
lr = [1,1,0,1,1,0,1,1]

plt.figure(figsize=(6.4,12.8))

for i,t in enumerate(tau):
    t *= DAY
    T = tf[i>>1]*DAY
    s = min_time(accel,rf[i>>1],t,T)
    t,r,th = s.y[:3]
    T = s.p[0]/DAY
    a = accel(t)
    phi = np.pi/2 - s.phi + th

    x,y = r*np.cos(th),r*np.sin(th)
    a,b = a*np.cos(phi),a*np.sin(phi)

    th = th[-1] + (t - t[-1])/np.sqrt(r[-1]**3)
    th = th[th-th[0]<=2*np.pi]        
    t = t[t<=2*np.pi]

    plt.subplot(4,2,i+1)
    plt.axis('equal')
    plt.plot(r[-1]*np.cos(th),r[-1]*np.sin(th),'--g')
    plt.plot(np.cos(t),np.sin(t),'--k')
    plt.plot(x,y,'b')

    x,y = x[::4],y[::4]
    a,b = a[::4],b[::4]
    plt.quiver(x,y,a,b,width=0.005,color='r')

    ha = 'right' if lr[i] else 'left'
    tx = plt.gca().transAxes
    plt.text(lr[i],0.07,r'$T$=%.2f days'%T,ha=ha,transform=tx)
    plt.text(lr[i],0.01,r'$\tau$=%d days'%tau[i],ha=ha,transform=tx)
    plt.title(title[i>>1], fontsize=14)

plt.tight_layout()
plt.savefig('fig6.eps')
plt.show()
