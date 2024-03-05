import numpy as np
from SystemResolution import *
from Constant import *
import matplotlib.pyplot as plt

# Compute the coefficient of B for the region 1 with the magnets

p = np.linspace(0, hm, n_mesh) # y
q = np.linspace(-L, L, n_mesh) # x
#p,q = np.meshgrid(p,q)

Mps = 0
'''
Bq0_1 = np.zeros_like(p)
Bqc_1 = np.zeros_like(p)
Bqs_1 = np.zeros_like(p)
Bpc_1 = np.zeros_like(p)
Bps_1 = np.zeros_like(p)
'''
Bq0_2 = 0
Bqc_2 = 0
Bqs_2 = 0
Bpc_2 = 0
Bps_2 = 0

p = hm

#Bq2 = np.zeros(len(q))
#Bp2 = np.zeros(len(p))

Bq2 = 0
Bp2 = 0

for i in range(1, n_harm+1):

    Mps = Mps + 2*L/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))
    omega_n = np.pi*i/L
    an_2 = b[0,i-1]
    bn_2 = b[1,i-1]
    cn_2 = b[2,i-1]
    dn_2 = b[3,i-1]

    Bqc_2 = an_2*np.exp(omega_n*p) - bn_2*np.exp(-omega_n*p)
    Bqs_2 = cn_2*np.exp(omega_n*p) - dn_2*np.exp(-omega_n*p)


    print(Bqc_2, Bqs_2)
    Bpc_2 = -cn_2*np.exp(omega_n*p) - dn_2*np.exp(-omega_n*p)
    Bps_2 = an_2*np.exp(omega_n*p) + bn_2*np.exp(-omega_n*p) + mu0*Mps

    Bq2 = Bq2 + Bqs_2*np.sin(omega_n*q) + Bqc_2*np.cos(omega_n*q)
    Bp2 = Bp2 + Bps_2*np.sin(omega_n*q) + Bpc_2*np.sin(omega_n*q)

plt.show()