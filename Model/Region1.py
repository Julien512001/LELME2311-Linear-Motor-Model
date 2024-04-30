import numpy as np
from SystemResolution import *
from Constant import *

# Compute the coefficient of B for the region 1 with the magnets


Mps = 0

Bq0_1 = np.zeros_like(P1)
Bqc_1 = np.zeros_like(P1)
Bqs_1 = np.zeros_like(P1)
Bpc_1 = np.zeros_like(P1)
Bps_1 = np.zeros_like(P1)

Bq1 = np.zeros_like(P1)
Bp1 = np.zeros_like(P1)
f = np.zeros_like(q)



for i in range(1, n_harm+1):

    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))

    omega_n = np.pi*i/tau_k
    an_1 = b[0,i-1]
    bn_1 = b[1,i-1]
    cn_1 = b[2,i-1]
    dn_1 = b[3,i-1]

    Bqc_1 = an_1*np.exp(omega_n*P1) - bn_1*np.exp(-omega_n*P1)
    Bqs_1 = cn_1*np.exp(omega_n*P1) - dn_1*np.exp(-omega_n*P1)

    Bpc_1 = -cn_1*np.exp(omega_n*P1) - dn_1*np.exp(-omega_n*P1)
    Bps_1 = an_1*np.exp(omega_n*P1) + bn_1*np.exp(-omega_n*P1) + mu0*Mps

    Bq1 = Bq1 + Bqs_1*np.sin(omega_n*Q1) + Bqc_1*np.cos(omega_n*Q1)
    Bp1 = Bp1 + Bpc_1*np.cos(omega_n*Q1) + Bps_1*np.sin(omega_n*Q1)