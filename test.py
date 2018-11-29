

import numpy as np
from scipy.integrate import odeint
from odesparse import odeints
from scipy.sparse import *

def deriv(y, t):
#    dy=np.zeros(y.size)
    yext=np.concatenate(([y[-1]],y,[y[0]]))
    dy=yext[2:]-2.0*y+yext[0:-2]
    return dy

N=400000

y0=np.zeros((N,))
y0[0]=1.0
y1=np.array(y0)

def l_ones(n):
    return [1 for x in xrange(n)]
JPat=lil_diags([l_ones(N-1), l_ones(N), l_ones(N-1)],[-1,0,1],(N,N))
JPat[0,N-1]=1
JPat[N-1,0]=1
#JPat=np.ones((5,5))
#JPat=csc_matrix((1000,1000))
#JPat[0,0]=1
#JPat[999,999]=1

#sol,op=odeint(deriv,y0,[0,0.1],full_output=1, rtol=1e-8, atol=1e-12)
#print sol[-1,:]

#print op

sol1,op=odeints(deriv,y1,[0,0.1],lrw=10000000,full_output=1, rtol=1e-8, atol=1e-12,JPat=JPat)
print sol1[-1,:]

print op

