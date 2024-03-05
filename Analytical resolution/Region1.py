import numpy as np
from SystemResolution import *
from Constant import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# Compute the coefficient of B for the region 1 with the magnets

p = np.linspace(0, hm, n_mesh) # y
q = np.linspace(-L, L, n_mesh) # x
p,q = np.meshgrid(p,q)

Mps = 0

Bq0_1 = np.zeros_like(p)
Bqc_1 = np.zeros_like(p)
Bqs_1 = np.zeros_like(p)
Bpc_1 = np.zeros_like(p)
Bps_1 = np.zeros_like(p)
'''
Bq0_1 = 0
Bqc_1 = 0
Bqs_1 = 0
Bpc_1 = 0
Bps_1 = 0
'''


#Bq1 = np.zeros(len(q))
#Bp1 = np.zeros(len(p))
Bq1 = np.zeros_like(p)
Bp1 = np.zeros_like(p)

for i in range(1, n_harm+1):

    Mps = Mps + 2*L/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))
    omega_n = np.pi*i/L
    an_1 = b[0,i-1]
    bn_1 = b[1,i-1]
    cn_1 = b[2,i-1]
    dn_1 = b[3,i-1]

    Bqc_1 = an_1*np.exp(omega_n*p) - bn_1*np.exp(-omega_n*p)
    Bqs_1 = cn_1*np.exp(omega_n*p) - dn_1*np.exp(-omega_n*p)

    Bpc_1 = -cn_1*np.exp(omega_n*p) - dn_1*np.exp(-omega_n*p)
    Bps_1 = an_1*np.exp(omega_n*p) + bn_1*np.exp(-omega_n*p) + mu0*Mps

    Bq1 = Bq1 + Bqs_1*np.sin(omega_n*q) + Bqc_1*np.cos(omega_n*q)
    Bp1 = Bp1 + Bpc_1*np.sin(omega_n*q) + Bps_1*np.sin(omega_n*q)


B = np.zeros_like(p)
for i in range(len(p)):
    for j in range(len(q)):
        B[i][j] = np.sqrt(Bq1[i][j]**2 + Bp1[i][j]**2)

fig = plt.figure()
plt.contourf(q,p,B)
plt.show()