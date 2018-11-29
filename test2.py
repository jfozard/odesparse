

import numpy as np
from scipy.integrate import odeint
from odesparse import odeints
from scipy.sparse import *

def deriv(y, t):
    dy=np.array([y[1],y[2],y[0]])
    return dy

N=3

y0=np.zeros((N,))
y0[0]=1.0
y1=np.array(y0)

JPat=lil_matrix((3,3))
JPat[0,1]=1
JPat[1,2]=1
JPat[2,0]=1

sol1,op=odeints(deriv,y1,[0,0.1],lrw=10000000,full_output=1, rtol=1e-8, atol=1e-12,JPat=JPat)
print sol1[-1,:]

print op

