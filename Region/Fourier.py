import numpy as np
import matplotlib.pyplot as plt
from Constant import *

P = np.linspace(0, hm, n_mesh) # y
Q = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p,q = np.meshgrid(P,Q)

print(L)
f = np.zeros_like(Q)
Mps = 0

for k in range(1,n_harm+1):
    #f = f + 2*L/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(1-1/2*d/L)))*np.sin(k*np.pi*x/L)
    #f = f + (-4*L/(np.pi*k))*np.sin(np.pi*d*k/(4*L))*np.sin(k*(np.pi - 3*np.pi*d/(4*L))) * np.sin(k*np.pi*x/L)
    #Mps = Mps + 2*L/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(d-L)/L))
    Mps = 2*Mp/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(d-L)/L))
    print(k)
    print("-----")
    print(Mps)
    """
    plt.figure()
    plt.plot(Q,f)
    plt.axvline(-L+d/2, c='r')
    plt.axvline(L-d/2, c='r')
    """
    f = f + Mps*np.sin(k*np.pi*Q/L)
    if (k == 2):
        print(f)


plt.figure()
plt.plot(Q,f)
plt.axvline(-L+d/2, c='r')
plt.axvline(L-d/2, c='r')

plt.show()
