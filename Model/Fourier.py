import numpy as np
import matplotlib.pyplot as plt
from Constant import *


f = np.zeros_like(q)


for i in range(1, n_harm+1):
    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))
    f += Mps*np.sin(np.pi*i/tau_k * q)

plt.figure()
plt.plot(q, f)
plt.show()