import numpy as np
import matplotlib.pyplot as plt
from Constant import *



x = np.linspace(-L,L, n_mesh)

f = np.zeros(len(x))
Mps = 0

for k in range(1,n_harm+1):
    #f = f + 2*L/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(1-1/2*d/L)))*np.sin(k*np.pi*x/L)
    #f = f + (-4*L/(np.pi*k))*np.sin(np.pi*d*k/(4*L))*np.sin(k*(np.pi - 3*np.pi*d/(4*L))) * np.sin(k*np.pi*x/L)
    f = f + 2*L/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(d-L)/L))*np.sin(k*np.pi*x/L)
    Mps = Mps + 2*L/(np.pi*k)*(np.cos(k*np.pi) - np.cos(k*np.pi*(d-L)/L))


plt.plot(x,f)
plt.axvline(-L+d/2, c='r')
plt.axvline(L-d/2, c='r')
plt.show()

print(Mps)