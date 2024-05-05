import numpy as np
import matplotlib.pyplot as plt
from Constant import *

q = np.linspace(-25e-3, 25e-3, 1000)
f = np.zeros_like(q)
e = 0

for i in range(1, n_harm+1):
    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))
    f += Mps*np.sin(np.pi*i/tau_k * q)



# Plotting
plt.figure(figsize=(8, 6))
plt.plot(q, f, color='b', label='Fourier')
plt.xlabel('q [mm]')
plt.ylabel('$M_p$')
plt.title('Magnetization Profile')
plt.grid(True)
plt.legend()
plt.show()