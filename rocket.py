import numpy as np
from scipy.integrate import solve_bvp,cumtrapz
from scipy.integrate._bvp import BVPResult

YEAR = 2*np.pi
DAY = YEAR/365.25

def max_radius(accel,tf,init=None,N=128,min=False):
    """ Optimize thrust direction so as to maximize
      radius of circular orbit in given time tf.
    Initially (t=0) and finally (t=tf), oribits are
      heliocentric and circular in ecliptic.
    Initially (r,theta)=(1,0);
    length unit = AU; time unit = year/2pi;
    gravitational const GM = 1 (M = solar mass)
    input:
      accel = thrust acceleration (function of time)
      tf = flight time (in year/2pi)
      init = initial guess as BVPResult object,
             must be output of max_radius
      N = number of grid points in time [0,tf]
          used only if init is not given
      min = True if minimizing rf (for inner planets)
    return:
      s = BVPResult object with following attributes:
      s.x = time t in [0,tf] (shape (N,))
            N is automatically adapted by solve_bvp
      s.y = solution of BVP (shape (6,N))
      s.y[0] = radial coordinate r
               s.y[0,0]=1, s.y[0,-1]=rf
      s.y[1] = radial velocity u
               s.y[1,0]=0, s.y[1,-1]=0
      s.y[2] = tangential velocity v
               s.y[2,0]=1, s.y[2,-1]=1/rf**0.5
      s.y[3:6] = costates of r,u,v
      s.th = ecliptic longitude theta (shape (N,))
             s.th[0]=0
      s.phi = thrust direction measured clockwise
              from tangential direction (shape (N,))
      s.tau = time of latest conjunction of
              earth and planet (0<tau<tf)
    reference:
      A.E.Bryson "Dynamic Optimization" section 3.4
    """
    e = 2*int(min) - 1

    def difeq(t,y):
        r,u,v,lr,lu,lv = y
        l = np.sqrt(lu**2 + lv**2)
        sin_phi = -lu/l
        cos_phi = -lv/l
        vr = v/r
        uvr = u*vr
        a = accel(t)
        return [u,
                v*vr - 1/r**2 + a*sin_phi,
                -uvr + a*cos_phi,
                lu*(vr**2 - 2/r**3) - lv*uvr/r,
                -lr + lv*vr,
                -2*lu*vr + lv*u/r]

    def bc(ya,yb):
        ra,ua,va = ya[:3]
        rb,ub,vb = yb[:3]
        lr,lu,lv = yb[3:]
        return [ra - 1,
                ua,
                va - 1,
                ub,
                vb - 1/np.sqrt(rb),
                lr - lv/2/np.sqrt(rb**3) - e]

    if isinstance(init, BVPResult):
        t = init.x/init.x[-1]*tf
        y = init.y
        s = solve_bvp(difeq, bc, t, y)
    else:
        x = np.linspace(0,1,N)
        y = [np.ones(N),
             np.zeros(N),
             np.ones(N),
             np.full(N,e),
             np.full(N,e),
             np.full(N,e)]

        t0 = 2 if min else 3
        dt = 0.5
        n = np.ceil(np.abs(tf-t0)/dt)
        for t in np.linspace(t0,tf,n+1):
            s = solve_bvp(difeq, bc, x*t, y)
            x = s.x/s.x[-1]
            y = s.y

    s.th = cumtrapz(s.y[2]/s.y[0], s.x, initial=0)
    s.phi = np.unwrap(np.arctan2(-s.y[4], -s.y[5]))
    om = 1/np.sqrt(s.y[0,-1]**3)
    s.tau = (s.th[-1] - om*tf)/(1-om)
    return s


def min_time(accel,rf,tau,tf=1,init=None,N=128):
    """ Optimize thrust direction so as to minimize
      flight time tf to circular orbit of radius rf.
    Initially (t=0) and finally (t=tf), oribits are
      heliocentric and circular in ecliptic.
    Initially (r,theta)=(1,0);
    length unit = AU; time unit = year/2pi;
    gravitational const GM = 1 (M = solar mass)
    input:
      accel = thrust acceleration (function of time)
      rf = radius of final circular orbit (in AU)
      tau = time of latest conjunction of
            earth and planet (0<tau<tf)
      tf = initial guess for flight time,
           used only if init is not given
      init = initial guess as BVPResult object,
             must be output of max_radius or min_time
      N = number of grid points in time [0,tf]
          used only if init is not given
    return:
      s = BVPResult object with following attributes:
      s.x = normalized time in [0,1] (shape (N,))
            N is automatically adapted by solve_bvp
      s.y = solution of BVP (shape (8,N))
      s.y[0] = time t in [0,tf]
      s.y[1] = radial coordinate r
               s.y[1,0]=1, s.y[1,-1]=rf
      s.y[2] = ecliptic longitude theta
               s.y[2,0]=0
      s.y[3] = radial velocity u
               s.y[3,0]=0, s.y[3,-1]=0
      s.y[4] = tangential velocity v
               s.y[4,0]=1, s.y[4,-1]=1/rf**0.5
      s.y[5:8] = costates of r,u,v
      s.p[0] = minimized flight time tf
      s.p[1] = costate of theta (scalar const)
      s.phi = thrust direction measured clockwise
              from tangential direction (shape (N,))
    reference:
      A.E.Bryson "Dynamic Optimization" section 4.5
    """
    i = not isinstance(init, BVPResult)
    if i or init.p is None:
        if i: s = max_radius(accel,tf,N=N,min=(rf<1))
        else: s = init
        s.p = [s.x[-1], 1]
        s.y = np.vstack((s.x,s.y[0],s.th,s.y[1:]))
        s.x = s.x/s.x[-1]

        dt = 0.2 if rf<1 else 0.5
        n = np.ceil(np.abs(tau - s.tau)/dt)
        for t in np.linspace(s.tau, tau, n+1):
            s = min_time(accel,rf,t,init=s)

        return s

    vf = 1/np.sqrt(rf)
    om = vf/rf
    th0 = (1-om)*tau

    def difeq(x,y,p):
        t,r,th,u,v,lr,lu,lv = y
        tf,lth = p
        l = np.sqrt(lu**2 + lv**2)
        sin_phi = -lu/l
        cos_phi = -lv/l
        vr = v/r
        uvr = u*vr
        a = accel(t)
        f = [np.ones_like(x),
             u,
             vr,
             v*vr - 1/r**2 + a*sin_phi,
             -uvr + a*cos_phi,
             lth*vr/r + lu*(vr**2 - 2/r**3) - lv*uvr/r,
             -lr + lv*vr,
             -lth/r - 2*lu*vr + lv*u/r]
        return tf*np.array(f)

    def bc(ya,yb,p):
        ta,ra,tha,ua,va = ya[:5]
        rb,thb,ub,vb = yb[1:5]
        lu,lv = yb[6:]
        tf = p[0]
        a = accel(tf)
        return [ta,
                ra - 1,
                tha,
                ua,
                va - 1,
                rb - rf,
                thb - (th0 + om*tf),
                ub,
                vb - vf,
                a*np.sqrt(lu**2 + lv**2) - 1]

    x = init.x
    y = init.y
    p = init.p
    s = solve_bvp(difeq, bc, x, y, p)
    s.phi = np.unwrap(np.arctan2(-s.y[6], -s.y[7]))
    return s
