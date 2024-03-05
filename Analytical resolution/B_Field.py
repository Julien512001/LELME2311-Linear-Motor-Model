import numpy as np
from Region1 import *
from Region2 import *

x = np.linspace(-L, L, n_mesh)
y = np.linspace(0, h, n_mesh)

'''
Bq_1 = np.zeros((len(x), len(y)))
Bp_1 = np.zeros((len(x), len(y)))

Bq_2 = np.zeros((len(x), len(y)))
Bp_2 = np.zeros((len(x), len(y)))
'''

Bq_1 = np.zeros(len(x))
Bp_1 = np.zeros(len(x))

Bq_2 = np.zeros(len(x))
Bp_2 = np.zeros(len(x))


for i in range(1, n_harm+1):
    omega_n = np.pi*i/L
    Bq_1 += Bqs_1*np.sin(omega_n*x) + Bqc_1*np.cos(omega_n*x)
    Bp_1 += Bpc_1*np.cos(omega_n*x) + Bps_1*np.sin(omega_n*x)

    Bq_2 += Bqs_2*np.sin(omega_n*x) + Bqc_2*np.cos(omega_n*x)
    Bp_2 += Bpc_2*np.cos(omega_n*x) + Bps_2*np.sin(omega_n*x)

B1 = np.zeros_like(Bq_1)
B2 = np.zeros_like(Bq_1)

for i in range(0,len(Bq_1)):
    B1[i] = np.linalg.norm([Bq_1[i], Bp_1[i]])
    B2[i] = np.linalg.norm([Bq_2[i], Bp_2[i]])

plt.plot(x,B1)
#plt.plot(x, B2)
plt.show()