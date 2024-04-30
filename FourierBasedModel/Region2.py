import numpy as np
from SystemResolution import *
from Constant import *

# Compute the coefficient of B for the region 1 with the magnets



Mps = 0

Bq0_2 = np.zeros_like(P2)
Bqc_2 = np.zeros_like(P2)
Bqs_2 = np.zeros_like(P2)
Bpc_2 = np.zeros_like(P2)
Bps_2 = np.zeros_like(P2)


Bq2 = np.zeros_like(P2)
Bp2 = np.zeros_like(P2)


for i in range(1, n_harm+1):

    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))
    
    omega_n = np.pi*i/tau_k
    an_2 = b[4,i-1]
    bn_2 = b[5,i-1]
    cn_2 = b[6,i-1]
    dn_2 = b[7,i-1]

    Bqc_2 = an_2*np.exp(omega_n*P2) - bn_2*np.exp(-omega_n*P2)
    Bqs_2 = cn_2*np.exp(omega_n*P2) - dn_2*np.exp(-omega_n*P2)

    Bpc_2 = -cn_2*np.exp(omega_n*P2) - dn_2*np.exp(-omega_n*P2)
    Bps_2 = an_2*np.exp(omega_n*P2) + bn_2*np.exp(-omega_n*P2)

    Bq2 = Bq2 + Bqs_2*np.sin(omega_n*Q2) + Bqc_2*np.cos(omega_n*Q2)
    Bp2 = Bp2 + Bpc_2*np.cos(omega_n*Q2) + Bps_2*np.sin(omega_n*Q2)