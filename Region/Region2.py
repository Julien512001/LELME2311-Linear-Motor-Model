import numpy as np
from SystemResolution import *
from Constant import *
import matplotlib.pyplot as plt

# Compute the coefficient of B for the region 1 with the magnets

P = np.linspace(hm, h, n_mesh) # y
Q = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p,q = np.meshgrid(P,Q)

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


Bq2 = np.zeros_like(p)
Bp2 = np.zeros_like(p)


for i in range(1, n_harm+1):

    Mps = 2*Mp/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))
    
    omega_n = np.pi*i/L
    an_2 = b[4,i-1]
    bn_2 = b[5,i-1]
    cn_2 = b[6,i-1]
    dn_2 = b[7,i-1]

    Bqc_2 = an_2*np.exp(omega_n*p) - bn_2*np.exp(-omega_n*p)
    Bqs_2 = cn_2*np.exp(omega_n*p) - dn_2*np.exp(-omega_n*p)

    Bpc_2 = -cn_2*np.exp(omega_n*p) - dn_2*np.exp(-omega_n*p)
    Bps_2 = an_2*np.exp(omega_n*p) + bn_2*np.exp(-omega_n*p)

    Bq2 = Bq2 + Bqs_2*np.sin(omega_n*q) + Bqc_2*np.cos(omega_n*q)
    Bp2 = Bp2 + Bpc_2*np.cos(omega_n*q) + Bps_2*np.sin(omega_n*q)


B2 = np.zeros_like(p)
for i in range(len(p)):
    for j in range(len(q)):
        B2[i][j] = np.sqrt(Bq2[i][j]**2 + Bp2[i][j]**2)
