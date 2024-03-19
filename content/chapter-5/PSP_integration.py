#PSP_integration.py
from __future__ import division

def trapzd(func, a, b, n):
    if (n == 1): 
        trapzd.s = 0.5*(b-a)*(func(a) + func(b))
        return trapzd.s
    else:
        it = 1 
        for j in range(1, n-1):
            it <<= 1
        tnm = float(it)
        del_ = (b-a)/tnm
        x = a+0.5*del_
        sum = 0.0
        for j in range(1, it+1): 
            sum += func(x)
            x += del_
        trapzd.s = 0.5*(trapzd.s + (b-a)*sum/tnm)
        return trapzd.s

def qsimp(func, a, b):
    EPS = 1.0e-6
    JMAX = 20
    ost=0.0
    os=0.0
    for j in range(1, JMAX+1):
        st = trapzd(func,a,b,j)
        s = (4.0*st-ost)/3.0
        if (j > 5):
            if (abs(s-os) < EPS*abs(os) or (s == 0.0 and os == 0.0)): 
                return s
        os=s
        ost=st