import numpy as np
import matplotlib.pyplot as plt

def BrownBridge(n,d):
    W = np.random.normal(loc=0,scale=1,size=(n,d))
    W[0,:] = 0
    B = np.cumsum(W,axis=0)
    t = np.linspace(0, 1, n)
    return B - np.outer(t,B[n-1,:])

def BesselBridge(n,d):
    BB = BrownBridge(n,d)
    return np.sqrt(np.square(BB).mean(axis=1))

def StochasticSTF(n,r,d=2):
    n1 = int(np.floor(n*r/(1+r)))
    n2 = n-n1+1
    B1 = BesselBridge(n1+1,d)
    B2 = BesselBridge(n2+1,d)
    STF = np.convolve(B1[0:n1], B2[0:n2], mode='full')
    return STF/np.sum(STF)

n=1000
r=1.0
STF = StochasticSTF(n, r)
plt.plot(STF)
plt.show()
