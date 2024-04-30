import numpy as np
import matplotlib.pyplot as plt
from Constant import *

P = np.linspace(0, hm, n_mesh) # y
Q = np.linspace(-L, L, n_mesh) # x
p,q = np.meshgrid(P,Q)

print(L)
f = np.zeros_like(Q)
Mps = 0

for i in range(1,n_harm+1):
    '''
    if i%2 == 0:
        Mps = 0.0
    else:
        Mps = 2*Mp/(i*np.pi) * (np.cos(np.pi*i*(d-L)/L) - np.cos(np.pi*i))
    '''
    Mps = 2*Mp/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))
    #Mps = -( 2*L*Mp*np.sin(np.pi*d*i/(4*L)) * np.sin(i*(np.pi - 3*np.pi*d/(4*L))) - 0.6366/L )/i
    #Mps = L*Mp*0.3183/(L*i) * ( (np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L)) + (np.cos(i*np.pi) - np.cos(i*np.pi*(L-d)/L)))
    #Mps = 2*L*Mp*0.3183/(L*i) * ( (np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L)))
    #Mps = 2*Mp/(i*np.pi) * (np.cos(np.pi*i) - np.cos(np.pi*i*(d-L)/L) )
    print("Mps[{}] = {}".format(i, Mps))
    """
    plt.figure()
    plt.plot(Q,f)
    plt.axvline(-L+d/2, c='r')
    plt.axvline(L-d/2, c='r')
    """
    f += Mps*np.sin(i*np.pi*Q/L)


plt.figure()
plt.plot(Q,f)
plt.axvline(-L+d/2, c='r')
plt.axvline(L-d/2, c='r')

plt.show()
